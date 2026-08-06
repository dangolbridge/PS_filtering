"""Global-assignment phase-singularity tracking and segment filtering.

This updated module improves the original ``rotor_track`` workflow by:

1. matching active tracks and current PS points with a global one-to-one
   assignment rather than greedy nearest-neighbour assignment;
2. splitting tracks first at temporal gaps and then at excessive rolling
   drift over a configurable time window;
3. evaluating and counting only valid segments, rather than counting every
   point from a track that contains at least one valid segment; and
4. reporting occupancy, time-gap, movement, and drift-split diagnostics.

The maximum allowed duration of one valid continuous segment remains:

    max_single_segment_fraction * (number_of_TCL_cycles * mean_TCL)

Track matching and rolling-drift checks use bounded shortest-path distance
along the mesh graph. Euclidean distance is retained only as a safe local
prefilter and for descriptive movement diagnostics.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

from .mesh_geodesic import MeshGeodesicGraph


MappedPSTimestep = tuple[float, np.ndarray]
SegmentKey = tuple[int, int]


@dataclass
class RotorTrack:
    """One spatially linked phase-singularity track."""

    id: int
    times: list[float] = field(default_factory=list)
    step_indices: list[int] = field(default_factory=list)
    vertices: list[int] = field(default_factory=list)
    coordinates: list[np.ndarray] = field(default_factory=list)

    def append(
        self,
        *,
        time: float,
        step_index: int,
        vertex: int,
        coordinate: np.ndarray,
    ) -> None:
        """Append one matched PS detection."""
        self.times.append(float(time))
        self.step_indices.append(int(step_index))
        self.vertices.append(int(vertex))
        self.coordinates.append(
            np.asarray(coordinate, dtype=float).copy()
        )

    @property
    def origin_coordinate(self) -> np.ndarray:
        return self.coordinates[0]

    @property
    def last_coordinate(self) -> np.ndarray:
        return self.coordinates[-1]

    @property
    def last_step_index(self) -> int:
        return self.step_indices[-1]


@dataclass(frozen=True)
class RotorTrackSegment:
    """One final track segment after time-gap and rolling-drift splitting."""

    track_id: int
    segment_index: int
    parent_temporal_segment_index: int
    rolling_piece_index: int
    started_after_rolling_drift_split: bool
    times: np.ndarray
    step_indices: np.ndarray
    vertices: np.ndarray
    coordinates: np.ndarray

    @property
    def start_time(self) -> float:
        return float(self.times[0])

    @property
    def end_time(self) -> float:
        return float(self.times[-1])

    @property
    def duration_time(self) -> float:
        # Keep the previously selected first-to-last timestamp definition.
        if self.times.size < 2:
            return 0.0
        return float(self.times[-1] - self.times[0])

    @property
    def n_hits(self) -> int:
        return int(self.times.size)


@dataclass
class RotorFilterResult:
    """Track-filtering output with valid-segment membership."""

    accepted_tracks: list[RotorTrack]
    rejected_tracks: list[RotorTrack]
    all_segments: list[RotorTrackSegment]
    valid_segments: list[RotorTrackSegment]
    invalid_segments: list[RotorTrackSegment]
    valid_segment_keys: set[SegmentKey]
    track_report: pd.DataFrame
    track_summary: pd.DataFrame
    segment_report: pd.DataFrame
    drift_split_report: pd.DataFrame


@dataclass
class RotorTrackingResult:
    """Complete per-cycle and full-range rotor-tracking output."""

    cycle_counts: dict[int, np.ndarray]
    cycle_intervals: list[tuple[float, float]]
    all_tracks: list[RotorTrack]
    accepted_tracks: list[RotorTrack]
    rejected_tracks: list[RotorTrack]
    valid_segments: list[RotorTrackSegment]
    invalid_segments: list[RotorTrackSegment]
    full_counts: np.ndarray
    normalized_counts: np.ndarray
    hit_counts: np.ndarray
    duration_map: np.ndarray
    track_report: pd.DataFrame
    track_summary: pd.DataFrame
    segment_report: pd.DataFrame
    drift_split_report: pd.DataFrame
    track_points: pd.DataFrame
    timestep_dt: float
    reference_tcl: float
    nominal_simulation_duration: float
    maximum_allowed_segment_time: float | None
    match_radius: float
    max_gap_steps: int
    gap_factor: float
    min_segment_time: float
    drift_window_time: float
    drift_radius_factor: float
    max_single_segment_fraction: float | None


def _prepare_mesh_coordinates(
    mesh_coordinates: np.ndarray,
) -> np.ndarray:
    """Validate mesh coordinates as an ``(n_vertices, 3)`` array."""
    mesh_coordinates = np.asarray(mesh_coordinates, dtype=float)

    if (
        mesh_coordinates.ndim != 2
        or mesh_coordinates.shape[1] != 3
    ):
        raise ValueError(
            "mesh_coordinates must have shape (n_vertices, 3); "
            f"received {mesh_coordinates.shape}."
        )

    if mesh_coordinates.shape[0] == 0:
        raise ValueError("The mesh contains no vertices.")

    if not np.all(np.isfinite(mesh_coordinates)):
        raise ValueError(
            "Mesh coordinates contain NaN or infinite values."
        )

    return mesh_coordinates


def _prepare_timesteps(
    timesteps: Sequence[MappedPSTimestep],
    n_vertices: int,
) -> list[MappedPSTimestep]:
    """Validate, sort, and deduplicate mapped PS timesteps."""
    prepared: list[MappedPSTimestep] = []

    for timestep_index, (time_value, vertices) in enumerate(timesteps):
        time_value = float(time_value)
        vertices = np.asarray(vertices, dtype=int).reshape(-1)

        if not np.isfinite(time_value):
            raise ValueError(
                f"Timestep {timestep_index} has a non-finite time."
            )

        if (
            np.any(vertices < 0)
            or np.any(vertices >= n_vertices)
        ):
            raise ValueError(
                f"Timestep {timestep_index} contains a vertex outside "
                f"0 to {n_vertices - 1}."
            )

        prepared.append((time_value, np.unique(vertices)))

    if not prepared:
        raise ValueError("No PS timesteps were provided.")

    prepared.sort(key=lambda item: item[0])
    times = np.asarray([item[0] for item in prepared], dtype=float)

    if times.size > 1 and np.any(np.diff(times) <= 0):
        raise ValueError(
            "PS timestep times must be strictly increasing."
        )

    return prepared


def estimate_timestep_dt(
    timesteps: Sequence[MappedPSTimestep],
) -> float:
    """Return the median positive PS timestep interval."""
    times = np.asarray(
        [float(item[0]) for item in timesteps],
        dtype=float,
    )

    if times.size < 2:
        return 1.0

    differences = np.diff(times)
    positive = differences[differences > 0]

    if positive.size == 0:
        raise ValueError(
            "Could not determine a positive PS timestep interval."
        )

    return float(np.median(positive))


def build_tracking_cycle_intervals(
    timesteps: Sequence[MappedPSTimestep],
    reference_tcl: float,
) -> list[tuple[float, float]]:
    """Divide the PS time range into TCL-length intervals."""
    reference_tcl = float(reference_tcl)

    if not np.isfinite(reference_tcl) or reference_tcl <= 0:
        raise ValueError(
            "reference_tcl must be positive and finite."
        )

    start_time = float(timesteps[0][0])
    end_time = float(timesteps[-1][0])

    if end_time <= start_time:
        return [(start_time, end_time)]

    starts = np.arange(start_time, end_time, reference_tcl)

    return [
        (
            float(start),
            float(min(start + reference_tcl, end_time)),
        )
        for start in starts
    ]


def _global_track_point_assignment(
    active_tracks: Sequence[RotorTrack],
    current_vertices: Sequence[int],
    mesh_coordinates: np.ndarray,
    geodesic_graph: MeshGeodesicGraph,
    match_radius: float,
) -> list[tuple[int, int, float]]:
    """Return globally optimal one-to-one geodesic assignments.

    Euclidean distance is used only as a safe prefilter: a mesh geodesic
    cannot be shorter than the straight-line distance. Candidate pairs that
    survive that prefilter are evaluated with bounded Dijkstra searches.

    The cost construction first maximizes the number of valid matches and,
    among assignments with equal cardinality, minimizes total geodesic
    distance.
    """
    n_tracks = len(active_tracks)
    n_points = len(current_vertices)

    if n_tracks == 0 or n_points == 0:
        return []

    current_vertex_array = np.asarray(current_vertices, dtype=int)
    point_coordinates = mesh_coordinates[current_vertex_array]
    distances = np.full((n_tracks, n_points), np.inf, dtype=float)

    for row_index, track in enumerate(active_tracks):
        euclidean = np.linalg.norm(
            point_coordinates - track.last_coordinate,
            axis=1,
        )
        candidate_indices = np.flatnonzero(euclidean <= match_radius)

        if candidate_indices.size == 0:
            continue

        target_vertices = current_vertex_array[candidate_indices]
        geodesic = geodesic_graph.distances_to_targets(
            int(track.vertices[-1]),
            target_vertices,
            max_distance=match_radius,
        )
        distances[row_index, candidate_indices] = geodesic

    cardinality_scale = max(n_tracks, n_points) + 1
    unmatched_cost = cardinality_scale * (match_radius + 1.0)
    invalid_cost = unmatched_cost * (cardinality_scale + 1.0)

    cost_matrix = np.full(
        (n_tracks, n_points + n_tracks),
        unmatched_cost,
        dtype=float,
    )
    valid = np.isfinite(distances) & (distances <= match_radius)
    cost_matrix[:, :n_points] = np.where(
        valid,
        distances,
        invalid_cost,
    )

    epsilon = np.finfo(float).eps * max(1.0, match_radius)
    if epsilon > 0:
        track_order = np.asarray(
            [track.id for track in active_tracks],
            dtype=float,
        )[:, None]
        point_order = np.arange(n_points, dtype=float)[None, :]
        cost_matrix[:, :n_points] += epsilon * (
            track_order + point_order / max(1, n_points)
        )

    row_indices, column_indices = linear_sum_assignment(cost_matrix)
    assignments: list[tuple[int, int, float]] = []

    for row_index, column_index in zip(row_indices, column_indices):
        if column_index >= n_points:
            continue

        distance = float(distances[row_index, column_index])
        if np.isfinite(distance) and distance <= match_radius:
            assignments.append(
                (int(row_index), int(column_index), distance)
            )

    assignments.sort(key=lambda item: active_tracks[item[0]].id)
    return assignments



def build_rotor_tracks(
    timesteps: Sequence[MappedPSTimestep],
    mesh_coordinates: np.ndarray,
    geodesic_graph: MeshGeodesicGraph,
    match_radius: float,
    max_gap_steps: int = 5,
) -> list[RotorTrack]:
    """Build tracks using global one-to-one geodesic assignment."""
    mesh_coordinates = _prepare_mesh_coordinates(mesh_coordinates)
    timesteps = _prepare_timesteps(
        timesteps,
        mesh_coordinates.shape[0],
    )
    match_radius = float(match_radius)
    max_gap_steps = int(max_gap_steps)

    if geodesic_graph.info.n_vertices != mesh_coordinates.shape[0]:
        raise ValueError(
            "The geodesic graph and coordinate array have different "
            "numbers of vertices."
        )

    if not np.isfinite(match_radius) or match_radius <= 0:
        raise ValueError("match_radius must be positive and finite.")
    if max_gap_steps < 1:
        raise ValueError("max_gap_steps must be at least 1.")

    active_tracks: list[RotorTrack] = []
    finished_tracks: list[RotorTrack] = []
    next_track_id = 0

    for step_index, (time_value, vertices_array) in enumerate(timesteps):
        still_matchable: list[RotorTrack] = []
        for track in active_tracks:
            if step_index - track.last_step_index > max_gap_steps:
                finished_tracks.append(track)
            else:
                still_matchable.append(track)
        active_tracks = still_matchable

        current_vertices = [int(vertex) for vertex in vertices_array]
        assignments = _global_track_point_assignment(
            active_tracks=active_tracks,
            current_vertices=current_vertices,
            mesh_coordinates=mesh_coordinates,
            geodesic_graph=geodesic_graph,
            match_radius=match_radius,
        )

        used_point_indices: set[int] = set()
        for active_index, point_index, _distance in assignments:
            track = active_tracks[active_index]
            vertex = current_vertices[point_index]
            track.append(
                time=time_value,
                step_index=step_index,
                vertex=vertex,
                coordinate=mesh_coordinates[vertex],
            )
            used_point_indices.add(point_index)

        for point_index, vertex in enumerate(current_vertices):
            if point_index in used_point_indices:
                continue
            track = RotorTrack(id=next_track_id)
            track.append(
                time=time_value,
                step_index=step_index,
                vertex=vertex,
                coordinate=mesh_coordinates[vertex],
            )
            active_tracks.append(track)
            next_track_id += 1

    finished_tracks.extend(active_tracks)
    finished_tracks.sort(key=lambda track: track.id)
    return finished_tracks



def _make_segment(
    track: RotorTrack,
    *,
    segment_index: int,
    parent_temporal_segment_index: int,
    rolling_piece_index: int,
    started_after_rolling_drift_split: bool,
    indices: np.ndarray,
) -> RotorTrackSegment:
    """Create one segment from selected positions in a track."""
    times = np.asarray(track.times, dtype=float)
    step_indices = np.asarray(track.step_indices, dtype=int)
    vertices = np.asarray(track.vertices, dtype=int)
    coordinates = np.asarray(track.coordinates, dtype=float)

    return RotorTrackSegment(
        track_id=track.id,
        segment_index=segment_index,
        parent_temporal_segment_index=parent_temporal_segment_index,
        rolling_piece_index=rolling_piece_index,
        started_after_rolling_drift_split=(
            started_after_rolling_drift_split
        ),
        times=times[indices].copy(),
        step_indices=step_indices[indices].copy(),
        vertices=vertices[indices].copy(),
        coordinates=coordinates[indices].copy(),
    )


def _temporal_segment_index_ranges(
    track: RotorTrack,
    gap_tolerance_time: float,
) -> list[np.ndarray]:
    """Return track-index arrays separated by large temporal gaps."""
    times = np.asarray(track.times, dtype=float)

    if times.size == 0:
        return []

    if np.any(np.diff(times) < 0):
        raise ValueError(
            f"Track {track.id} times are not sorted."
        )

    breaks = np.where(
        np.diff(times) > gap_tolerance_time
    )[0]
    starts = np.concatenate(([0], breaks + 1))
    ends = np.concatenate((breaks, [times.size - 1]))

    return [
        np.arange(start, end + 1, dtype=int)
        for start, end in zip(starts, ends)
    ]


def _split_indices_by_rolling_drift(
    track: RotorTrack,
    indices: np.ndarray,
    *,
    geodesic_graph: MeshGeodesicGraph,
    drift_window_time: float,
    maximum_rolling_drift_distance: float,
    parent_temporal_segment_index: int,
) -> tuple[list[np.ndarray], list[dict[str, object]]]:
    """Split a temporal segment at rolling geodesic-drift violations."""
    if indices.size == 0:
        return [], []

    times = np.asarray(track.times, dtype=float)
    vertices = np.asarray(track.vertices, dtype=int)
    pieces: list[np.ndarray] = []
    split_rows: list[dict[str, object]] = []
    piece_start_position = 0
    rolling_piece_index = 1

    while piece_start_position < indices.size:
        split_position: int | None = None

        for local_position in range(piece_start_position + 1, indices.size):
            current_global_index = int(indices[local_position])
            current_time = float(times[current_global_index])
            target_time = current_time - drift_window_time

            candidate_local_positions = np.arange(
                piece_start_position,
                local_position,
                dtype=int,
            )
            candidate_global_indices = indices[candidate_local_positions]
            candidate_times = times[candidate_global_indices]
            eligible = np.where(candidate_times <= target_time)[0]

            if eligible.size == 0:
                continue

            anchor_candidate_position = int(eligible[-1])
            anchor_local_position = int(
                candidate_local_positions[anchor_candidate_position]
            )
            anchor_global_index = int(indices[anchor_local_position])
            anchor_vertex = int(vertices[anchor_global_index])
            current_vertex = int(vertices[current_global_index])

            displacement = geodesic_graph.distance(
                anchor_vertex,
                current_vertex,
                max_distance=maximum_rolling_drift_distance,
            )

            if not np.isfinite(displacement) or (
                displacement > maximum_rolling_drift_distance
            ):
                split_position = local_position
                previous_global_index = int(indices[local_position - 1])
                split_rows.append(
                    {
                        "track_id": track.id,
                        "parent_temporal_segment_index": (
                            parent_temporal_segment_index
                        ),
                        "rolling_piece_before_split": rolling_piece_index,
                        "split_time": current_time,
                        "anchor_time": float(times[anchor_global_index]),
                        "actual_window_time": float(
                            current_time - times[anchor_global_index]
                        ),
                        "rolling_displacement": displacement,
                        "maximum_allowed_rolling_displacement": (
                            maximum_rolling_drift_distance
                        ),
                        "anchor_vertex": anchor_vertex,
                        "previous_vertex": int(vertices[previous_global_index]),
                        "current_vertex": current_vertex,
                        "distance_exceeded_or_disconnected": True,
                    }
                )
                break

        if split_position is None:
            pieces.append(indices[piece_start_position:].copy())
            break

        pieces.append(indices[piece_start_position:split_position].copy())
        piece_start_position = split_position
        rolling_piece_index += 1

    return pieces, split_rows



def split_track_into_final_segments(
    track: RotorTrack,
    *,
    geodesic_graph: MeshGeodesicGraph,
    gap_tolerance_time: float,
    drift_window_time: float,
    maximum_rolling_drift_distance: float,
) -> tuple[list[RotorTrackSegment], list[dict[str, object]]]:
    """Apply temporal-gap and rolling geodesic-drift splitting."""
    temporal_ranges = _temporal_segment_index_ranges(
        track,
        gap_tolerance_time,
    )
    final_segments: list[RotorTrackSegment] = []
    drift_split_rows: list[dict[str, object]] = []
    final_segment_index = 1

    for temporal_segment_index, indices in enumerate(
        temporal_ranges,
        start=1,
    ):
        rolling_pieces, split_rows = _split_indices_by_rolling_drift(
            track,
            indices,
            geodesic_graph=geodesic_graph,
            drift_window_time=drift_window_time,
            maximum_rolling_drift_distance=(
                maximum_rolling_drift_distance
            ),
            parent_temporal_segment_index=temporal_segment_index,
        )
        drift_split_rows.extend(split_rows)

        for rolling_piece_index, piece_indices in enumerate(
            rolling_pieces,
            start=1,
        ):
            final_segments.append(
                _make_segment(
                    track,
                    segment_index=final_segment_index,
                    parent_temporal_segment_index=temporal_segment_index,
                    rolling_piece_index=rolling_piece_index,
                    started_after_rolling_drift_split=(
                        rolling_piece_index > 1
                    ),
                    indices=piece_indices,
                )
            )
            final_segment_index += 1

    return final_segments, drift_split_rows



def _movement_statistics(
    coordinates: np.ndarray,
) -> dict[str, float]:
    """Calculate Euclidean movement diagnostics."""
    coordinates = np.asarray(coordinates, dtype=float)

    if coordinates.shape[0] == 0:
        return {
            "path_length": 0.0,
            "net_displacement": 0.0,
            "maximum_displacement_from_origin": 0.0,
            "mean_step_distance": 0.0,
            "maximum_step_distance": 0.0,
        }

    displacement = np.linalg.norm(
        coordinates - coordinates[0],
        axis=1,
    )

    if coordinates.shape[0] == 1:
        return {
            "path_length": 0.0,
            "net_displacement": 0.0,
            "maximum_displacement_from_origin": 0.0,
            "mean_step_distance": 0.0,
            "maximum_step_distance": 0.0,
        }

    steps = np.linalg.norm(
        np.diff(coordinates, axis=0),
        axis=1,
    )

    return {
        "path_length": float(steps.sum()),
        "net_displacement": float(
            np.linalg.norm(coordinates[-1] - coordinates[0])
        ),
        "maximum_displacement_from_origin": float(
            displacement.max()
        ),
        "mean_step_distance": float(steps.mean()),
        "maximum_step_distance": float(steps.max()),
    }


def _time_gap_statistics(
    times: np.ndarray,
    timestep_dt: float,
) -> dict[str, float | int]:
    """Calculate occupancy and temporal-gap diagnostics."""
    times = np.asarray(times, dtype=float)
    n_hits = int(times.size)

    if n_hits == 0:
        return {
            "sampled_span_time": 0.0,
            "observed_time": 0.0,
            "occupancy_fraction": 0.0,
            "minimum_time_gap": np.nan,
            "mean_time_gap": np.nan,
            "median_time_gap": np.nan,
            "maximum_time_gap": np.nan,
            "n_internal_missing_steps": 0,
        }

    duration_time = (
        float(times[-1] - times[0])
        if n_hits > 1
        else 0.0
    )
    sampled_span_time = duration_time + timestep_dt
    observed_time = n_hits * timestep_dt
    occupancy_fraction = min(
        1.0,
        observed_time / sampled_span_time,
    )

    if n_hits == 1:
        return {
            "sampled_span_time": sampled_span_time,
            "observed_time": observed_time,
            "occupancy_fraction": occupancy_fraction,
            "minimum_time_gap": np.nan,
            "mean_time_gap": np.nan,
            "median_time_gap": np.nan,
            "maximum_time_gap": np.nan,
            "n_internal_missing_steps": 0,
        }

    gaps = np.diff(times)
    estimated_steps = np.maximum(
        np.rint(gaps / timestep_dt).astype(int),
        1,
    )
    missing_steps = np.maximum(estimated_steps - 1, 0)

    return {
        "sampled_span_time": sampled_span_time,
        "observed_time": observed_time,
        "occupancy_fraction": occupancy_fraction,
        "minimum_time_gap": float(gaps.min()),
        "mean_time_gap": float(gaps.mean()),
        "median_time_gap": float(np.median(gaps)),
        "maximum_time_gap": float(gaps.max()),
        "n_internal_missing_steps": int(missing_steps.sum()),
    }


def _maximum_rolling_displacement(
    segment: RotorTrackSegment,
    drift_window_time: float,
    geodesic_graph: MeshGeodesicGraph,
    maximum_search_distance: float,
) -> float:
    """Return the largest bounded geodesic rolling displacement."""
    times = segment.times
    vertices = segment.vertices
    maximum = 0.0

    for current_index in range(1, times.size):
        target_time = times[current_index] - drift_window_time
        eligible = np.where(times[:current_index] <= target_time)[0]
        if eligible.size == 0:
            continue

        anchor_index = int(eligible[-1])
        distance = geodesic_graph.distance(
            int(vertices[anchor_index]),
            int(vertices[current_index]),
            max_distance=maximum_search_distance,
        )
        if not np.isfinite(distance):
            return float("inf")
        maximum = max(maximum, float(distance))

    return maximum



def filter_rotor_tracks(
    tracks: Sequence[RotorTrack],
    *,
    geodesic_graph: MeshGeodesicGraph,
    timestep_dt: float,
    nominal_simulation_duration: float,
    min_segment_time: float = 120.0,
    gap_factor: float = 3.0,
    drift_window_time: float = 120.0,
    drift_radius_factor: float = 2.0,
    base_radius: float = 2000.0,
    max_single_segment_fraction: float | None = 0.90,
) -> RotorFilterResult:
    """Split, evaluate, and classify track segments.

    A track is accepted when it contains at least one valid final segment.
    A final segment is valid when it reaches ``min_segment_time`` and does
    not exceed the simulation-duration cap. Rolling-drift violations split
    segments; they do not automatically reject either side of the split.
    """
    timestep_dt = float(timestep_dt)
    nominal_simulation_duration = float(
        nominal_simulation_duration
    )
    min_segment_time = float(min_segment_time)
    gap_factor = float(gap_factor)
    drift_window_time = float(drift_window_time)
    drift_radius_factor = float(drift_radius_factor)
    base_radius = float(base_radius)

    positive = {
        "timestep_dt": timestep_dt,
        "nominal_simulation_duration": (
            nominal_simulation_duration
        ),
        "gap_factor": gap_factor,
        "drift_window_time": drift_window_time,
        "drift_radius_factor": drift_radius_factor,
        "base_radius": base_radius,
    }

    for name, value in positive.items():
        if not np.isfinite(value) or value <= 0:
            raise ValueError(
                f"{name} must be positive and finite."
            )

    if not np.isfinite(min_segment_time) or min_segment_time < 0:
        raise ValueError(
            "min_segment_time must be non-negative and finite."
        )

    if max_single_segment_fraction is None:
        maximum_allowed_segment_time = None
    else:
        max_single_segment_fraction = float(
            max_single_segment_fraction
        )
        if (
            not np.isfinite(max_single_segment_fraction)
            or max_single_segment_fraction <= 0
        ):
            raise ValueError(
                "max_single_segment_fraction must be positive and "
                "finite or None."
            )
        maximum_allowed_segment_time = (
            max_single_segment_fraction
            * nominal_simulation_duration
        )

    gap_tolerance_time = gap_factor * timestep_dt
    maximum_rolling_drift_distance = (
        drift_radius_factor * base_radius
    )

    accepted_tracks: list[RotorTrack] = []
    rejected_tracks: list[RotorTrack] = []
    all_segments: list[RotorTrackSegment] = []
    valid_segments: list[RotorTrackSegment] = []
    invalid_segments: list[RotorTrackSegment] = []
    valid_segment_keys: set[SegmentKey] = set()
    track_rows: list[dict[str, object]] = []
    segment_rows: list[dict[str, object]] = []
    drift_split_rows: list[dict[str, object]] = []

    for track in tracks:
        segments, track_drift_splits = split_track_into_final_segments(
            track,
            geodesic_graph=geodesic_graph,
            gap_tolerance_time=gap_tolerance_time,
            drift_window_time=drift_window_time,
            maximum_rolling_drift_distance=(
                maximum_rolling_drift_distance
            ),
        )
        all_segments.extend(segments)
        drift_split_rows.extend(track_drift_splits)

        track_valid_segments: list[RotorTrackSegment] = []
        track_invalid_reasons: list[str] = []

        for segment in segments:
            movement = _movement_statistics(segment.coordinates)
            gap_stats = _time_gap_statistics(
                segment.times,
                timestep_dt,
            )
            maximum_rolling_displacement = (
                _maximum_rolling_displacement(
                    segment,
                    drift_window_time,
                    geodesic_graph,
                    maximum_rolling_drift_distance,
                )
            )

            invalid_reasons: list[str] = []

            if segment.duration_time < min_segment_time:
                invalid_reasons.append("too_short")

            exceeds_cap = (
                maximum_allowed_segment_time is not None
                and segment.duration_time
                > maximum_allowed_segment_time
            )
            if exceeds_cap:
                invalid_reasons.append(
                    "segment_exceeds_simulation_cap"
                )

            valid = len(invalid_reasons) == 0

            if valid:
                valid_segments.append(segment)
                track_valid_segments.append(segment)
                valid_segment_keys.add(
                    (track.id, segment.segment_index)
                )
                reason = "valid"
            else:
                invalid_segments.append(segment)
                reason = ";".join(invalid_reasons)
                for invalid_reason in invalid_reasons:
                    if invalid_reason not in track_invalid_reasons:
                        track_invalid_reasons.append(invalid_reason)

            segment_rows.append(
                {
                    "track_id": track.id,
                    "segment_index": segment.segment_index,
                    "parent_temporal_segment_index": (
                        segment.parent_temporal_segment_index
                    ),
                    "rolling_piece_index": (
                        segment.rolling_piece_index
                    ),
                    "started_after_rolling_drift_split": (
                        segment.started_after_rolling_drift_split
                    ),
                    "valid": valid,
                    "reason": reason,
                    "start_time": segment.start_time,
                    "end_time": segment.end_time,
                    "duration_time": segment.duration_time,
                    "n_hits": segment.n_hits,
                    "sampled_span_time": gap_stats[
                        "sampled_span_time"
                    ],
                    "observed_time": gap_stats["observed_time"],
                    "occupancy_fraction": gap_stats[
                        "occupancy_fraction"
                    ],
                    "minimum_time_gap": gap_stats[
                        "minimum_time_gap"
                    ],
                    "mean_time_gap": gap_stats["mean_time_gap"],
                    "median_time_gap": gap_stats[
                        "median_time_gap"
                    ],
                    "maximum_time_gap": gap_stats[
                        "maximum_time_gap"
                    ],
                    "n_internal_missing_steps": gap_stats[
                        "n_internal_missing_steps"
                    ],
                    "path_length": movement["path_length"],
                    "net_displacement": movement[
                        "net_displacement"
                    ],
                    "maximum_displacement_from_segment_origin": (
                        movement["maximum_displacement_from_origin"]
                    ),
                    "mean_step_distance": movement[
                        "mean_step_distance"
                    ],
                    "maximum_step_distance": movement[
                        "maximum_step_distance"
                    ],
                    "maximum_rolling_displacement": (
                        maximum_rolling_displacement
                    ),
                    "maximum_allowed_rolling_displacement": (
                        maximum_rolling_drift_distance
                    ),
                    "maximum_allowed_segment_time": (
                        maximum_allowed_segment_time
                    ),
                }
            )

        accepted = len(track_valid_segments) > 0
        if accepted:
            accepted_tracks.append(track)
            track_reason = "accepted"
        else:
            rejected_tracks.append(track)
            track_reason = (
                "no_valid_segment"
                if not track_invalid_reasons
                else "no_valid_segment:" + ";".join(
                    track_invalid_reasons
                )
            )

        times = np.asarray(track.times, dtype=float)
        coordinates = np.asarray(track.coordinates, dtype=float)
        movement = _movement_statistics(coordinates)
        valid_durations = np.asarray(
            [segment.duration_time for segment in track_valid_segments],
            dtype=float,
        )
        valid_hits = int(
            sum(segment.n_hits for segment in track_valid_segments)
        )

        track_rows.append(
            {
                "track_id": track.id,
                "accepted": accepted,
                "reason": track_reason,
                "n_hits": int(times.size),
                "n_final_segments": len(segments),
                "n_valid_segments": len(track_valid_segments),
                "n_invalid_segments": (
                    len(segments) - len(track_valid_segments)
                ),
                "n_rolling_drift_splits": len(track_drift_splits),
                "first_time": (
                    float(times[0]) if times.size else np.nan
                ),
                "last_time": (
                    float(times[-1]) if times.size else np.nan
                ),
                "lifetime_span": (
                    float(times[-1] - times[0])
                    if times.size > 1
                    else (0.0 if times.size == 1 else np.nan)
                ),
                "valid_n_hits": valid_hits,
                "valid_duration_sum": float(
                    valid_durations.sum()
                ) if valid_durations.size else 0.0,
                "longest_valid_segment_time": float(
                    valid_durations.max()
                ) if valid_durations.size else 0.0,
                "maximum_allowed_segment_time": (
                    maximum_allowed_segment_time
                ),
                "path_length": movement["path_length"],
                "net_displacement": movement[
                    "net_displacement"
                ],
                "maximum_displacement_from_origin": movement[
                    "maximum_displacement_from_origin"
                ],
                "mean_step_distance": movement[
                    "mean_step_distance"
                ],
                "maximum_step_distance": movement[
                    "maximum_step_distance"
                ],
            }
        )

    track_report = pd.DataFrame(track_rows)
    segment_report = pd.DataFrame(segment_rows)
    drift_split_report = pd.DataFrame(drift_split_rows)

    summary_rows: list[dict[str, object]] = [
        {"reason": "accepted", "count": len(accepted_tracks)},
        {"reason": "rejected", "count": len(rejected_tracks)},
        {"reason": "valid_segments", "count": len(valid_segments)},
        {"reason": "invalid_segments", "count": len(invalid_segments)},
        {
            "reason": "rolling_drift_splits",
            "count": len(drift_split_rows),
        },
    ]

    if not segment_report.empty:
        invalid_reason_counts: dict[str, int] = {}
        invalid_rows = segment_report.loc[
            ~segment_report["valid"].astype(bool)
        ]
        for reason_text in invalid_rows["reason"].astype(str):
            for reason in reason_text.split(";"):
                invalid_reason_counts[reason] = (
                    invalid_reason_counts.get(reason, 0) + 1
                )
        for reason, count in sorted(invalid_reason_counts.items()):
            summary_rows.append(
                {"reason": reason, "count": count}
            )

    track_summary = pd.DataFrame(
        summary_rows,
        columns=["reason", "count"],
    )

    return RotorFilterResult(
        accepted_tracks=accepted_tracks,
        rejected_tracks=rejected_tracks,
        all_segments=all_segments,
        valid_segments=valid_segments,
        invalid_segments=invalid_segments,
        valid_segment_keys=valid_segment_keys,
        track_report=track_report,
        track_summary=track_summary,
        segment_report=segment_report,
        drift_split_report=drift_split_report,
    )



def segments_to_vertex_counts(
    valid_segments: Sequence[RotorTrackSegment],
    n_vertices: int,
    *,
    unique_per_track: bool = True,
) -> np.ndarray:
    """Create a vertex map using valid segments only."""
    n_vertices = int(n_vertices)

    if n_vertices < 1:
        raise ValueError(
            "n_vertices must be at least 1."
        )

    counts = np.zeros(n_vertices, dtype=int)

    if unique_per_track:
        vertices_by_track: dict[int, set[int]] = {}
        for segment in valid_segments:
            vertices_by_track.setdefault(segment.track_id, set()).update(
                int(vertex) for vertex in segment.vertices
            )

        for vertices in vertices_by_track.values():
            for vertex in vertices:
                if vertex < 0 or vertex >= n_vertices:
                    raise ValueError(
                        f"Invalid vertex {vertex} in valid segment."
                    )
                counts[vertex] += 1
    else:
        for segment in valid_segments:
            for vertex in segment.vertices:
                vertex = int(vertex)
                if vertex < 0 or vertex >= n_vertices:
                    raise ValueError(
                        f"Invalid vertex {vertex} in valid segment."
                    )
                counts[vertex] += 1

    return counts


def segments_to_hit_counts(
    valid_segments: Sequence[RotorTrackSegment],
    n_vertices: int,
) -> np.ndarray:
    """Count every valid-segment detection at each vertex."""
    return segments_to_vertex_counts(
        valid_segments,
        n_vertices,
        unique_per_track=False,
    )


def build_track_points_table(
    tracks: Sequence[RotorTrack],
    all_segments: Sequence[RotorTrackSegment],
    valid_segment_keys: set[SegmentKey],
) -> pd.DataFrame:
    """Return one row per point with final-segment validity labels."""
    track_accepted_ids = {
        track_id for track_id, _ in valid_segment_keys
    }
    rows: list[dict[str, object]] = []

    for segment in all_segments:
        key = (segment.track_id, segment.segment_index)
        segment_valid = key in valid_segment_keys
        track_accepted = segment.track_id in track_accepted_ids

        for point_index, (
            step_index,
            time_value,
            vertex,
            coordinate,
        ) in enumerate(
            zip(
                segment.step_indices,
                segment.times,
                segment.vertices,
                segment.coordinates,
            )
        ):
            rows.append(
                {
                    "track_id": segment.track_id,
                    "track_accepted": track_accepted,
                    "segment_index": segment.segment_index,
                    "segment_valid": segment_valid,
                    "parent_temporal_segment_index": (
                        segment.parent_temporal_segment_index
                    ),
                    "rolling_piece_index": (
                        segment.rolling_piece_index
                    ),
                    "point_index_in_segment": point_index,
                    "step_index": int(step_index),
                    "time": float(time_value),
                    "vertex": int(vertex),
                    "x": float(coordinate[0]),
                    "y": float(coordinate[1]),
                    "z": float(coordinate[2]),
                }
            )

    return pd.DataFrame(rows)


def analyze_rotor_tracking(
    mapped_timesteps: Sequence[MappedPSTimestep],
    mesh_coordinates: np.ndarray,
    geodesic_graph: MeshGeodesicGraph,
    reference_tcl: float,
    match_radius: float,
    *,
    max_gap_steps: int = 5,
    min_segment_time: float = 120.0,
    gap_factor: float = 3.0,
    drift_window_time: float = 120.0,
    drift_radius_factor: float = 2.0,
    max_single_segment_fraction: float | None = 0.90,
    unique_per_track: bool = True,
) -> RotorTrackingResult:
    """Run cycle-by-cycle and full-range moving PS tracking."""
    mesh_coordinates = _prepare_mesh_coordinates(
        mesh_coordinates
    )
    timesteps = _prepare_timesteps(
        mapped_timesteps,
        mesh_coordinates.shape[0],
    )
    if geodesic_graph.info.n_vertices != mesh_coordinates.shape[0]:
        raise ValueError(
            "The geodesic graph and coordinate array have different "
            "numbers of vertices."
        )

    reference_tcl = float(reference_tcl)

    if not np.isfinite(reference_tcl) or reference_tcl <= 0:
        raise ValueError(
            "reference_tcl must be positive and finite."
        )

    timestep_dt = estimate_timestep_dt(timesteps)
    cycle_intervals = build_tracking_cycle_intervals(
        timesteps,
        reference_tcl,
    )
    nominal_simulation_duration = (
        len(cycle_intervals) * reference_tcl
    )
    maximum_allowed_segment_time = (
        None
        if max_single_segment_fraction is None
        else float(max_single_segment_fraction)
        * nominal_simulation_duration
    )

    cycle_counts: dict[int, np.ndarray] = {}

    for cycle_number, (cycle_start, cycle_end) in enumerate(
        cycle_intervals,
        start=1,
    ):
        is_last = cycle_number == len(cycle_intervals)

        cycle_timesteps = [
            item
            for item in timesteps
            if cycle_start <= item[0]
            and (
                item[0] <= cycle_end
                if is_last
                else item[0] < cycle_end
            )
        ]

        if not cycle_timesteps:
            continue

        cycle_tracks = build_rotor_tracks(
            cycle_timesteps,
            mesh_coordinates,
            geodesic_graph,
            match_radius,
            max_gap_steps,
        )
        cycle_filter = filter_rotor_tracks(
            cycle_tracks,
            geodesic_graph=geodesic_graph,
            timestep_dt=timestep_dt,
            nominal_simulation_duration=(
                nominal_simulation_duration
            ),
            min_segment_time=min_segment_time,
            gap_factor=gap_factor,
            drift_window_time=drift_window_time,
            drift_radius_factor=drift_radius_factor,
            base_radius=match_radius,
            max_single_segment_fraction=(
                max_single_segment_fraction
            ),
        )
        cycle_counts[cycle_number] = segments_to_vertex_counts(
            cycle_filter.valid_segments,
            mesh_coordinates.shape[0],
            unique_per_track=unique_per_track,
        )

    all_tracks = build_rotor_tracks(
        timesteps,
        mesh_coordinates,
        geodesic_graph,
        match_radius,
        max_gap_steps,
    )
    full_filter = filter_rotor_tracks(
        all_tracks,
        geodesic_graph=geodesic_graph,
        timestep_dt=timestep_dt,
        nominal_simulation_duration=(
            nominal_simulation_duration
        ),
        min_segment_time=min_segment_time,
        gap_factor=gap_factor,
        drift_window_time=drift_window_time,
        drift_radius_factor=drift_radius_factor,
        base_radius=match_radius,
        max_single_segment_fraction=(
            max_single_segment_fraction
        ),
    )

    full_counts = segments_to_vertex_counts(
        full_filter.valid_segments,
        mesh_coordinates.shape[0],
        unique_per_track=unique_per_track,
    )
    hit_counts = segments_to_hit_counts(
        full_filter.valid_segments,
        mesh_coordinates.shape[0],
    )
    duration_map = hit_counts.astype(float) * timestep_dt

    maximum_count = int(full_counts.max()) if full_counts.size else 0
    normalized_counts = (
        full_counts.astype(float) / maximum_count
        if maximum_count > 0
        else full_counts.astype(float)
    )

    track_points = build_track_points_table(
        all_tracks,
        full_filter.all_segments,
        full_filter.valid_segment_keys,
    )

    return RotorTrackingResult(
        cycle_counts=cycle_counts,
        cycle_intervals=cycle_intervals,
        all_tracks=all_tracks,
        accepted_tracks=full_filter.accepted_tracks,
        rejected_tracks=full_filter.rejected_tracks,
        valid_segments=full_filter.valid_segments,
        invalid_segments=full_filter.invalid_segments,
        full_counts=full_counts,
        normalized_counts=normalized_counts,
        hit_counts=hit_counts,
        duration_map=duration_map,
        track_report=full_filter.track_report,
        track_summary=full_filter.track_summary,
        segment_report=full_filter.segment_report,
        drift_split_report=full_filter.drift_split_report,
        track_points=track_points,
        timestep_dt=timestep_dt,
        reference_tcl=reference_tcl,
        nominal_simulation_duration=(
            nominal_simulation_duration
        ),
        maximum_allowed_segment_time=(
            maximum_allowed_segment_time
        ),
        match_radius=float(match_radius),
        max_gap_steps=int(max_gap_steps),
        gap_factor=float(gap_factor),
        min_segment_time=float(min_segment_time),
        drift_window_time=float(drift_window_time),
        drift_radius_factor=float(drift_radius_factor),
        max_single_segment_fraction=(
            None
            if max_single_segment_fraction is None
            else float(max_single_segment_fraction)
        ),
    )
