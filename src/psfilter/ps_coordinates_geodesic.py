"""Coordinate-based phase-singularity counting and filtering.

This module reproduces the intent of the ``PS_coords`` workflow:

1. Divide mapped PS timesteps into TCL-length cycles.
2. Count a PS point only when it is spatially connected to a point from a
   recent timestep.
3. Build per-cycle and full-range vertex count maps.
4. Filter vertex histories by duration and continuity criteria.

Spatial continuation is evaluated with bounded shortest-path distance along
the mesh graph. Euclidean distance is used only as a safe local prefilter:
because a mesh path cannot be shorter than the straight-line distance,
vertices farther than ``radius`` in Euclidean space cannot satisfy the
geodesic continuation criterion.

"""

from __future__ import annotations

from collections import deque
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from .mesh_geodesic import MeshGeodesicGraph


MappedPSTimestep = tuple[float, np.ndarray]


@dataclass
class PSHitResult:
    """Unfiltered spatial-continuity counts for one timestep range."""

    counts: np.ndarray
    first_hit_time: np.ndarray
    last_hit_time: np.ndarray
    hit_times: list[list[float]]


@dataclass
class PSCoordinateResult:
    """Complete coordinate-based PS analysis output."""

    cycle_counts: dict[int, np.ndarray]
    cycle_intervals: list[tuple[float, float]]
    total_counts: np.ndarray
    filtered_counts: np.ndarray
    normalized_counts: np.ndarray
    accepted_mask: np.ndarray
    segment_report: pd.DataFrame
    lifetime_report: pd.DataFrame
    timestep_dt: float
    reference_tcl: float
    simulation_duration: float
    maximum_allowed_segment_time: float


def _prepare_mesh_coordinates(mesh_coordinates: np.ndarray) -> np.ndarray:
    """Validate mesh coordinates."""
    mesh_coordinates = np.asarray(mesh_coordinates, dtype=float)

    if mesh_coordinates.ndim != 2 or mesh_coordinates.shape[1] != 3:
        raise ValueError(
            "mesh_coordinates must have shape (n_vertices, 3); "
            f"received {mesh_coordinates.shape}."
        )

    if mesh_coordinates.shape[0] == 0:
        raise ValueError("The mesh contains no vertices.")

    if not np.all(np.isfinite(mesh_coordinates)):
        raise ValueError("Mesh coordinates contain NaN or infinite values.")

    return mesh_coordinates


def _prepare_timesteps(
    timesteps: Sequence[MappedPSTimestep],
    n_vertices: int,
) -> list[MappedPSTimestep]:
    """Validate and sort mapped PS timesteps."""
    prepared: list[MappedPSTimestep] = []

    for index, (time_value, vertices) in enumerate(timesteps):
        time_value = float(time_value)
        vertices = np.asarray(vertices, dtype=int).reshape(-1)

        if not np.isfinite(time_value):
            raise ValueError(f"Timestep {index} has a non-finite time.")

        if np.any(vertices < 0) or np.any(vertices >= n_vertices):
            raise ValueError(
                f"Timestep {index} contains a vertex outside 0 to "
                f"{n_vertices - 1}."
            )

        # A coordinate file may map multiple PS coordinates to the same
        # vertex at one timestep. Count that vertex once per timestep.
        prepared.append((time_value, np.unique(vertices)))

    prepared.sort(key=lambda item: item[0])

    if not prepared:
        raise ValueError("No PS timesteps were provided.")

    times = np.asarray([time for time, _ in prepared], dtype=float)

    if times.size > 1 and np.any(np.diff(times) <= 0):
        raise ValueError("PS timestep times must be strictly increasing.")

    return prepared


def estimate_timestep_dt(
    timesteps: Sequence[MappedPSTimestep],
) -> float:
    """Estimate the median interval between PS timesteps."""
    times = np.asarray([float(time) for time, _ in timesteps], dtype=float)

    if times.size < 2:
        return 1.0

    differences = np.diff(times)
    positive = differences[differences > 0]

    if positive.size == 0:
        raise ValueError("Could not determine a positive PS timestep interval.")

    return float(np.median(positive))


def build_cycle_intervals(
    timesteps: Sequence[MappedPSTimestep],
    reference_tcl: float,
) -> list[tuple[float, float]]:
    """Create consecutive TCL-length intervals across the PS time range."""
    reference_tcl = float(reference_tcl)

    if not np.isfinite(reference_tcl) or reference_tcl <= 0:
        raise ValueError("reference_tcl must be a positive finite number.")

    start_time = float(timesteps[0][0])
    end_time = float(timesteps[-1][0])

    if end_time <= start_time:
        return [(start_time, end_time)]

    cycle_starts = np.arange(start_time, end_time, reference_tcl)

    return [
        (float(start), float(min(start + reference_tcl, end_time)))
        for start in cycle_starts
    ]


def _validate_geodesic_graph(
    mesh_coordinates: np.ndarray,
    geodesic_graph: MeshGeodesicGraph,
) -> None:
    """Ensure that the graph and coordinate array describe the same mesh."""
    if not isinstance(geodesic_graph, MeshGeodesicGraph):
        raise TypeError(
            "geodesic_graph must be a MeshGeodesicGraph instance."
        )

    n_vertices = mesh_coordinates.shape[0]

    if geodesic_graph.info.n_vertices != n_vertices:
        raise ValueError(
            "The geodesic graph and mesh coordinate array have different "
            f"vertex counts: {geodesic_graph.info.n_vertices} versus "
            f"{n_vertices}."
        )


def count_continuing_ps_hits(
    timesteps: Sequence[MappedPSTimestep],
    mesh_coordinates: np.ndarray,
    geodesic_graph: MeshGeodesicGraph,
    radius: float,
    history_steps: int = 5,
) -> PSHitResult:
    """Count PS detections connected through the mesh to recent detections.

    A current vertex is counted when its shortest-path distance along mesh
    edges to at least one vertex from the previous ``history_steps`` timesteps
    is no greater than ``radius``.

    Before running bounded Dijkstra, previous vertices farther than ``radius``
    in straight-line distance are removed. This prefilter is safe because a
    path along mesh edges cannot be shorter than the Euclidean displacement.
    """
    mesh_coordinates = _prepare_mesh_coordinates(mesh_coordinates)
    n_vertices = mesh_coordinates.shape[0]
    timesteps = _prepare_timesteps(timesteps, n_vertices)
    _validate_geodesic_graph(mesh_coordinates, geodesic_graph)

    radius = float(radius)
    history_steps = int(history_steps)

    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("radius must be a positive finite number.")

    if history_steps < 1:
        raise ValueError("history_steps must be at least 1.")

    counts = np.zeros(n_vertices, dtype=int)
    first_hit_time = np.full(n_vertices, np.nan, dtype=float)
    last_hit_time = np.full(n_vertices, np.nan, dtype=float)
    hit_times: list[list[float]] = [[] for _ in range(n_vertices)]

    history: deque[set[int]] = deque(maxlen=history_steps)

    for time_value, current_vertices in timesteps:
        previous_vertices = (
            set().union(*history)
            if history
            else set()
        )

        if previous_vertices:
            previous_array = np.fromiter(
                sorted(previous_vertices),
                dtype=int,
            )
            previous_coordinates = mesh_coordinates[previous_array]
        else:
            previous_array = np.empty(0, dtype=int)
            previous_coordinates = np.empty((0, 3), dtype=float)

        for vertex_value in current_vertices:
            vertex = int(vertex_value)

            if previous_array.size == 0:
                continue

            # Safe lower-bound prefilter before mesh-graph searches.
            euclidean_distances = np.linalg.norm(
                previous_coordinates - mesh_coordinates[vertex],
                axis=1,
            )
            candidate_vertices = previous_array[
                euclidean_distances <= radius
            ]

            if candidate_vertices.size == 0:
                continue

            geodesic_distances = geodesic_graph.distances_to_targets(
                source=vertex,
                targets=candidate_vertices,
                max_distance=radius,
            )

            if np.any(np.isfinite(geodesic_distances)):
                counts[vertex] += 1

                if np.isnan(first_hit_time[vertex]):
                    first_hit_time[vertex] = time_value

                last_hit_time[vertex] = time_value
                hit_times[vertex].append(float(time_value))

        history.append(
            set(int(vertex) for vertex in current_vertices)
        )

    return PSHitResult(
        counts=counts,
        first_hit_time=first_hit_time,
        last_hit_time=last_hit_time,
        hit_times=hit_times,
    )

def _split_hit_times(
    times: np.ndarray,
    gap_tolerance: float,
) -> list[np.ndarray]:
    """Split sorted hit times into continuous segments."""
    if times.size == 0:
        return []

    breaks = np.where(np.diff(times) > gap_tolerance)[0]
    starts = np.concatenate(([0], breaks + 1))
    ends = np.concatenate((breaks, [times.size - 1]))

    return [times[start:end + 1] for start, end in zip(starts, ends)]


def filter_vertex_ps_histories(
    hit_result: PSHitResult,
    timestep_dt: float,
    simulation_duration: float,
    *,
    min_segment_time: float = 120.0,
    gap_factor: float = 3.0,
    hard_cap_fraction: float = 0.90,
    reject_single_segment: bool = False,
) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    """Filter per-vertex PS histories.

    A vertex is rejected when any one continuous PS segment lasts longer
    than ``hard_cap_fraction * simulation_duration``.

    ``simulation_duration`` is calculated by :func:`analyze_ps_coordinates`
    as ``number_of_cycles * reference_tcl``. The cap is therefore applied to
    the duration of each individual PS segment, not to the cumulative number
    of hits across several separate segments.

    The old drift test is not applied here because each per-vertex history
    stores detections for one fixed mesh vertex. Spatial drift belongs to the
    later moving-track analysis.
    """
    timestep_dt = float(timestep_dt)
    simulation_duration = float(simulation_duration)
    min_segment_time = float(min_segment_time)
    gap_factor = float(gap_factor)
    hard_cap_fraction = float(hard_cap_fraction)

    if not np.isfinite(timestep_dt) or timestep_dt <= 0:
        raise ValueError("timestep_dt must be positive and finite.")

    if (
        not np.isfinite(simulation_duration)
        or simulation_duration <= 0
    ):
        raise ValueError(
            "simulation_duration must be positive and finite."
        )

    if not np.isfinite(min_segment_time) or min_segment_time < 0:
        raise ValueError("min_segment_time must be non-negative and finite.")

    if not np.isfinite(gap_factor) or gap_factor <= 0:
        raise ValueError("gap_factor must be positive and finite.")

    if (
        not np.isfinite(hard_cap_fraction)
        or hard_cap_fraction <= 0
        or hard_cap_fraction > 1
    ):
        raise ValueError(
            "hard_cap_fraction must be greater than 0 and no greater than 1."
        )

    n_vertices = hit_result.counts.size
    accepted_mask = np.zeros(n_vertices, dtype=bool)
    segment_rows: list[dict[str, int | float | str | bool]] = []
    lifetime_rows: list[dict[str, int | float | bool]] = []

    gap_tolerance = gap_factor * timestep_dt
    maximum_allowed_segment_time = (
        hard_cap_fraction * simulation_duration
    )

    for vertex in range(n_vertices):
        times = np.asarray(hit_result.hit_times[vertex], dtype=float)
        times.sort()

        segments = _split_hit_times(times, gap_tolerance)
        segment_durations = np.asarray(
            [
                float(segment[-1] - segment[0])
                if segment.size > 1
                else 0.0
                for segment in segments
            ],
            dtype=float,
        )
        segment_counts = np.asarray(
            [segment.size for segment in segments],
            dtype=int,
        )

        n_hits = int(times.size)
        n_segments = len(segments)
        longest_segment_time = (
            float(segment_durations.max())
            if segment_durations.size
            else 0.0
        )
        lifetime_span = (
            float(times[-1] - times[0])
            if times.size > 1
            else (0.0 if times.size == 1 else float("nan"))
        )
        estimated_hit_time = float(n_hits * timestep_dt)
        exceeds_simulation_cap = bool(
            segment_durations.size
            and np.any(
                segment_durations > maximum_allowed_segment_time
            )
        )

        rejection_reasons: list[str] = []

        if n_hits == 0:
            rejection_reasons.append("no_hits")

        if reject_single_segment and n_segments == 1:
            rejection_reasons.append("single_segment")

        if exceeds_simulation_cap:
            rejection_reasons.append(
                "segment_exceeds_simulation_cap"
            )

        if n_hits > 0 and longest_segment_time < min_segment_time:
            rejection_reasons.append("too_short")

        accepted = len(rejection_reasons) == 0
        accepted_mask[vertex] = accepted

        segment_rows.append(
            {
                "vertex": vertex,
                "n_hits": n_hits,
                "estimated_hit_time": estimated_hit_time,
                "n_segments": n_segments,
                "segment_counts": ";".join(
                    str(int(value)) for value in segment_counts
                ),
                "segment_durations": ";".join(
                    f"{float(value):.12g}" for value in segment_durations
                ),
                "longest_segment_time": longest_segment_time,
                "lifetime_span": lifetime_span,
                "simulation_duration": simulation_duration,
                "hard_cap_fraction": hard_cap_fraction,
                "maximum_allowed_segment_time": (
                    maximum_allowed_segment_time
                ),
                "exceeds_simulation_cap": exceeds_simulation_cap,
                "accepted": accepted,
                "reason": (
                    "accepted"
                    if accepted
                    else ";".join(rejection_reasons)
                ),
            }
        )

        first_hit = hit_result.first_hit_time[vertex]
        last_hit = hit_result.last_hit_time[vertex]

        lifetime_rows.append(
            {
                "vertex": vertex,
                "first_hit_t": first_hit,
                "last_hit_t": last_hit,
                "lifetime": (
                    float(last_hit - first_hit)
                    if np.isfinite(first_hit) and np.isfinite(last_hit)
                    else float("nan")
                ),
                "n_hits": n_hits,
                "accepted": accepted,
            }
        )

    return (
        accepted_mask,
        pd.DataFrame(segment_rows),
        pd.DataFrame(lifetime_rows),
    )

def analyze_ps_coordinates(
    mapped_timesteps: Sequence[MappedPSTimestep],
    mesh_coordinates: np.ndarray,
    geodesic_graph: MeshGeodesicGraph,
    reference_tcl: float,
    radius: float,
    *,
    history_steps: int = 5,
    min_segment_time: float = 120.0,
    gap_factor: float = 3.0,
    hard_cap_fraction: float = 0.90,
    reject_single_segment: bool = False,
) -> PSCoordinateResult:
    """Run coordinate-based PS analysis using mesh-geodesic continuity."""
    mesh_coordinates = _prepare_mesh_coordinates(mesh_coordinates)
    _validate_geodesic_graph(mesh_coordinates, geodesic_graph)
    timesteps = _prepare_timesteps(
        mapped_timesteps,
        mesh_coordinates.shape[0],
    )
    timestep_dt = estimate_timestep_dt(timesteps)
    cycle_intervals = build_cycle_intervals(timesteps, reference_tcl)

    # The filtering cap is defined relative to the nominal simulation
    # duration: number of TCL cycles multiplied by the mean TCL.
    simulation_duration = len(cycle_intervals) * float(reference_tcl)
    maximum_allowed_segment_time = (
        hard_cap_fraction * simulation_duration
    )

    cycle_counts: dict[int, np.ndarray] = {}

    for cycle_number, (cycle_start, cycle_end) in enumerate(
        cycle_intervals,
        start=1,
    ):
        is_last_cycle = cycle_number == len(cycle_intervals)

        if is_last_cycle:
            cycle_timesteps = [
                item
                for item in timesteps
                if cycle_start <= item[0] <= cycle_end
            ]
        else:
            cycle_timesteps = [
                item
                for item in timesteps
                if cycle_start <= item[0] < cycle_end
            ]

        if not cycle_timesteps:
            continue

        cycle_hit_result = count_continuing_ps_hits(
            timesteps=cycle_timesteps,
            mesh_coordinates=mesh_coordinates,
            geodesic_graph=geodesic_graph,
            radius=radius,
            history_steps=history_steps,
        )
        cycle_counts[cycle_number] = cycle_hit_result.counts

    total_hit_result = count_continuing_ps_hits(
        timesteps=timesteps,
        mesh_coordinates=mesh_coordinates,
        geodesic_graph=geodesic_graph,
        radius=radius,
        history_steps=history_steps,
    )

    accepted_mask, segment_report, lifetime_report = (
        filter_vertex_ps_histories(
            hit_result=total_hit_result,
            timestep_dt=timestep_dt,
            simulation_duration=simulation_duration,
            min_segment_time=min_segment_time,
            gap_factor=gap_factor,
            hard_cap_fraction=hard_cap_fraction,
            reject_single_segment=reject_single_segment,
        )
    )

    filtered_counts = total_hit_result.counts.copy()
    filtered_counts[~accepted_mask] = 0

    maximum = int(filtered_counts.max()) if filtered_counts.size else 0
    normalized_counts = (
        filtered_counts.astype(float) / maximum
        if maximum > 0
        else filtered_counts.astype(float)
    )

    return PSCoordinateResult(
        cycle_counts=cycle_counts,
        cycle_intervals=cycle_intervals,
        total_counts=total_hit_result.counts,
        filtered_counts=filtered_counts,
        normalized_counts=normalized_counts,
        accepted_mask=accepted_mask,
        segment_report=segment_report,
        lifetime_report=lifetime_report,
        timestep_dt=timestep_dt,
        reference_tcl=float(reference_tcl),
        simulation_duration=simulation_duration,
        maximum_allowed_segment_time=maximum_allowed_segment_time,
    )
