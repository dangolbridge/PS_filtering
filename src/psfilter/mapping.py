"""Spatial mapping utilities for phase-singularity coordinates."""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Sequence

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


PSTimestep = tuple[float, np.ndarray]
MappedPSTimestep = tuple[float, np.ndarray]


@dataclass
class PSMappingResult:
    """Result of mapping PS coordinates to nearest mesh vertices."""

    timesteps: list[MappedPSTimestep]
    report: pd.DataFrame
    n_input_points: int
    n_mapped_points: int
    n_dropped_points: int
    maximum_distance: float


def _prepare_mesh_coordinates(mesh_coordinates: np.ndarray) -> np.ndarray:
    """Validate mesh coordinates as an ``(n_vertices, 3)`` array."""
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


def map_ps_coordinates_to_vertices(
    ps_timesteps: Sequence[PSTimestep],
    mesh_coordinates: np.ndarray,
    max_distance: float | None = None,
) -> PSMappingResult:
    """Map every PS coordinate to its nearest mesh vertex.

    Parameters
    ----------
    ps_timesteps
        Sequence of ``(time, coordinates)`` items returned by
        :func:`psfilter.io.read_pts_t_file`. Each coordinate array must have
        shape ``(n_points, 3)``.
    mesh_coordinates
        Mesh vertex coordinates with shape ``(n_vertices, 3)``.
    max_distance
        Optional maximum accepted nearest-vertex distance. Points farther
        away are excluded. When omitted, every point is mapped.

    Returns
    -------
    PSMappingResult
        Mapped timesteps and a row-by-row mapping report.
    """
    mesh_coordinates = _prepare_mesh_coordinates(mesh_coordinates)

    if max_distance is not None:
        max_distance = float(max_distance)
        if not np.isfinite(max_distance) or max_distance < 0:
            raise ValueError(
                "max_distance must be a non-negative finite number or None."
            )

    tree = cKDTree(mesh_coordinates)
    mapped_timesteps: list[MappedPSTimestep] = []
    report_rows: list[dict[str, int | float | bool]] = []

    n_input = 0
    n_mapped = 0
    maximum_distance = 0.0

    for timestep_index, (time_value, point_coordinates) in enumerate(ps_timesteps):
        time_value = float(time_value)
        point_coordinates = np.asarray(point_coordinates, dtype=float)

        if point_coordinates.size == 0:
            point_coordinates = np.empty((0, 3), dtype=float)

        if point_coordinates.ndim != 2 or point_coordinates.shape[1] != 3:
            raise ValueError(
                "Each PS coordinate array must have shape (n_points, 3); "
                f"timestep {timestep_index} has shape {point_coordinates.shape}."
            )

        if not np.all(np.isfinite(point_coordinates)):
            raise ValueError(
                f"PS coordinates at timestep {timestep_index} contain "
                "NaN or infinite values."
            )

        mapped_vertices: list[int] = []

        if point_coordinates.shape[0] > 0:
            distances, vertex_indices = tree.query(point_coordinates, k=1)
            distances = np.asarray(distances, dtype=float).reshape(-1)
            vertex_indices = np.asarray(vertex_indices, dtype=int).reshape(-1)

            for point_index, (distance, vertex_index) in enumerate(
                zip(distances, vertex_indices)
            ):
                n_input += 1
                maximum_distance = max(maximum_distance, float(distance))
                accepted = (
                    max_distance is None
                    or float(distance) <= max_distance
                )

                if accepted:
                    mapped_vertices.append(int(vertex_index))
                    n_mapped += 1

                report_rows.append(
                    {
                        "timestep_index": timestep_index,
                        "time": time_value,
                        "point_index": point_index,
                        "vertex": int(vertex_index),
                        "distance": float(distance),
                        "accepted": bool(accepted),
                    }
                )

        mapped_timesteps.append(
            (
                time_value,
                np.asarray(mapped_vertices, dtype=int),
            )
        )

    report = pd.DataFrame(
        report_rows,
        columns=[
            "timestep_index",
            "time",
            "point_index",
            "vertex",
            "distance",
            "accepted",
        ],
    )

    return PSMappingResult(
        timesteps=mapped_timesteps,
        report=report,
        n_input_points=n_input,
        n_mapped_points=n_mapped,
        n_dropped_points=n_input - n_mapped,
        maximum_distance=float(maximum_distance),
    )
