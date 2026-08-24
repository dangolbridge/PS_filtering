"""Command-line interface for the psfilter analysis package.

The CLI exposes the existing numerical modules through five subcommands:

``tcl``
    Calculate TCL and, optionally, activation-pattern outputs.
``ps-detect``
    Create a phase-singularity ``.pts_t`` file with igbhead/igbfilament.
``ps-coords``
    Run coordinate-based PS continuity and filtering using either backend.
``rotor-track``
    Run globally assigned rotor tracking using either backend.
``all``
    Execute the complete workflow while reusing mesh coordinates and mapped
    PS timesteps, and build one geodesic graph only when requested.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import pandas as pd

from .activation import analyze_activation_patterns, plot_cycle_activation_figures
from .io import (
    read_mean_tcl,
    read_pts_file,
    read_pts_t_file,
    read_transmembrane_file,
    write_array,
    write_dataframe,
    write_tcl_summary,
)
from .mapping import PSMappingResult, map_ps_coordinates_to_vertices
from .mesh_geodesic import MeshGeodesicGraph
from .phase_singularity import ensure_phase_singularity_file
from .preprocessing import ensure_transmembrane_file
from .ps_coordinates import (
    PSCoordinateResult as EuclideanPSCoordinateResult,
    analyze_ps_coordinates as analyze_ps_coordinates_euclidean,
)
from .ps_coordinates_geodesic import (
    PSCoordinateResult as GeodesicPSCoordinateResult,
    analyze_ps_coordinates as analyze_ps_coordinates_geodesic,
)
from .rotor_tracking import (
    RotorTrackingResult as EuclideanRotorTrackingResult,
    analyze_rotor_tracking as analyze_rotor_tracking_euclidean,
)
from .rotor_tracking_geodesic import (
    RotorTrackingResult as GeodesicRotorTrackingResult,
    analyze_rotor_tracking as analyze_rotor_tracking_geodesic,
)
from .tcl import TCLResult, calculate_tcl



PSCoordinateResultType = EuclideanPSCoordinateResult | GeodesicPSCoordinateResult
RotorTrackingResultType = (
    EuclideanRotorTrackingResult | GeodesicRotorTrackingResult
)
DISTANCE_MODES = ("geodesic", "euclidean")

TCL_OUTPUTS = [
    "tcl_peaks_AP.csv",
    "tcl_per_signal_AP.csv",
    "global_tcl_AP.txt",
    "tcl_summary.csv",
]
ACTIVATION_OUTPUTS = [
    "cycle_activations.csv",
    "cycle_summary.csv",
    "global_activation_order.csv",
    "activation_sequences.csv",
    "sequence_counts.csv",
]
PS_COORD_OUTPUTS = [
    "rotor_counts_ptst.dat",
    "rotor_counts_normalized_ptst.dat",
    "rotor_segments_filtered.csv",
    "rotor_lifetimes.csv",
    "ps_coords_summary.csv",
]
ROTOR_TRACK_OUTPUTS = [
    "rotor_track_counts_ptst.dat",
    "rotor_track_counts_normalized_ptst.dat",
    "rotor_tracks_counts_report.csv",
    "rotor_summary.csv",
    "rotor_tracking_run_summary.csv",
]


@dataclass(frozen=True)
class PipelinePaths:
    """Resolved paths used by the ``all`` subcommand."""

    work_dir: Path
    transmembrane_file: Path
    vm_igb: Path
    mesh: Path
    mesh_points: Path
    mesh_elements: Path
    points_time: Path
    tcl_output_dir: Path
    ps_output_dir: Path
    transformed_points: Path | None
    query_points: Path | None
    node_indices: Path | None


def _nonempty_file(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def _outputs_ready(directory: Path, filenames: Sequence[str]) -> bool:
    return all(_nonempty_file(directory / filename) for filename in filenames)


def _read_labels_file(path: Path | None) -> list[str] | None:
    if path is None:
        return None
    if not path.is_file():
        raise FileNotFoundError(f"Signal-label file not found: {path}")
    labels = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not labels:
        raise ValueError(f"Signal-label file is empty: {path}")
    return labels


def _resolve_relative(path: Path | None, base: Path) -> Path | None:
    if path is None:
        return None
    return path if path.is_absolute() else base / path


def _resolve_pipeline_paths(args: argparse.Namespace) -> PipelinePaths:
    work_dir = Path(args.work_dir)
    mesh = _resolve_relative(Path(args.mesh), work_dir)
    assert mesh is not None

    tcl_output_dir = _resolve_relative(args.tcl_output_dir, work_dir)
    if tcl_output_dir is None:
        tcl_output_dir = work_dir / "tcl_results"

    ps_output_dir = _resolve_relative(args.ps_output_dir, work_dir)
    if ps_output_dir is None:
        ps_output_dir = work_dir / "PS_results"

    transmembrane_file = _resolve_relative(args.input, work_dir)
    if transmembrane_file is None:
        transmembrane_file = work_dir / "transmembrane_v.dat"

    vm_igb = _resolve_relative(args.vm_igb, work_dir)
    if vm_igb is None:
        vm_igb = work_dir / "vm.igb"

    mesh_points = _resolve_relative(args.mesh_points, work_dir)
    if mesh_points is None:
        mesh_points = Path(f"{mesh}.pts")

    mesh_elements = _resolve_relative(args.mesh_elements, work_dir)
    if mesh_elements is None:
        mesh_elements = Path(f"{mesh}.elem")

    points_time = _resolve_relative(args.points_time, work_dir)
    if points_time is None:
        points_time = ps_output_dir / f"{mesh.name}.pts_t"

    return PipelinePaths(
        work_dir=work_dir,
        transmembrane_file=transmembrane_file,
        vm_igb=vm_igb,
        mesh=mesh,
        mesh_points=mesh_points,
        mesh_elements=mesh_elements,
        points_time=points_time,
        tcl_output_dir=tcl_output_dir,
        ps_output_dir=ps_output_dir,
        transformed_points=_resolve_relative(args.transformed_points, work_dir),
        query_points=_resolve_relative(args.query_points, work_dir),
        node_indices=_resolve_relative(args.node_indices, work_dir),
    )


def _require_values(**values: object) -> None:
    missing = [name for name, value in values.items() if value is None]
    if missing:
        joined = ", ".join(f"--{name.replace('_', '-')}" for name in missing)
        raise ValueError(f"Missing required arguments for this operation: {joined}")


def _ensure_transmembrane(
    *,
    input_path: Path,
    vm_igb: Path | None,
    transformed_points: Path | None,
    mesh: Path | None,
    query_points: Path | None,
    node_indices: Path | None,
    meshtool: Path | str | None,
    igbextract: Path | str | None,
    search_radius: float,
    overwrite: bool,
) -> Path:
    # ``--overwrite`` controls numerical outputs, not the extracted input.
    # Reuse an existing nonempty transmembrane file to avoid an unnecessary
    # meshtool/igbextract pass.
    if _nonempty_file(input_path):
        return input_path

    _require_values(
        vm_igb=vm_igb,
        transformed_points=transformed_points,
        mesh=mesh,
        query_points=query_points,
        node_indices=node_indices,
        meshtool=meshtool,
        igbextract=igbextract,
    )

    return ensure_transmembrane_file(
        output_path=input_path,
        vm_igb_path=vm_igb,
        transformed_points_path=transformed_points,
        mesh_path=mesh,
        query_points_path=query_points,
        node_indices_path=node_indices,
        meshtool_executable=meshtool,
        igbextract_executable=igbextract,
        search_radius=search_radius,
        overwrite=False,
    )


def _write_tcl_outputs(
    result: TCLResult,
    output_dir: Path,
    *,
    run_activation: bool,
    fraction_window: float,
    fraction_cluster: float,
    max_groups: int | None,
    stop_on_max_groups: bool,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    write_dataframe(output_dir / "tcl_peaks_AP.csv", result.peaks_table)
    write_dataframe(
        output_dir / "tcl_per_signal_AP.csv",
        result.per_signal_statistics,
    )
    write_tcl_summary(
        output_dir / "global_tcl_AP.txt",
        mean_rr=result.reference_mean_tcl,
        std_rr=result.reference_std_tcl,
    )

    summary = pd.DataFrame(
        [
            {
                "reference_index": result.reference_index,
                "reference_label": result.signal_labels[result.reference_index],
                "reference_mean_tcl": result.reference_mean_tcl,
                "reference_std_tcl": result.reference_std_tcl,
                "reference_n_intervals": len(result.reference_intervals),
                "pooled_mean_tcl": result.pooled_mean_tcl,
                "pooled_std_tcl": result.pooled_std_tcl,
                "pooled_n_intervals": len(result.pooled_intervals),
                "dt": result.dt,
                "derivative_threshold": result.derivative_threshold,
                "min_peak_distance": result.min_peak_distance,
                "derivative_per_time": result.derivative_per_time,
            }
        ]
    )
    write_dataframe(output_dir / "tcl_summary.csv", summary)

    if not run_activation:
        return

    activation = analyze_activation_patterns(
        tcl_result=result,
        fraction_window=fraction_window,
        fraction_cluster=fraction_cluster,
        max_groups=max_groups,
        stop_on_max_groups=stop_on_max_groups,
    )
    write_dataframe(
        output_dir / "cycle_activations.csv",
        activation.cycle_activations,
    )
    write_dataframe(output_dir / "cycle_summary.csv", activation.cycle_summary)
    write_dataframe(
        output_dir / "global_activation_order.csv",
        activation.global_activation_order,
    )
    write_dataframe(
        output_dir / "activation_sequences.csv",
        activation.activation_sequences,
    )
    write_dataframe(output_dir / "sequence_counts.csv", activation.sequence_counts)

    plot_cycle_activation_figures(
    activation,
    output_dir=output_dir / "activation_cycle_figures",
    time_mode="delay",
)
    
def _run_tcl_stage(
    *,
    input_path: Path,
    output_dir: Path,
    dt: float,
    reference_index: int,
    derivative_threshold: float,
    min_peak_distance: float,
    labels_file: Path | None,
    derivative_per_time: bool,
    run_activation: bool,
    fraction_window: float,
    fraction_cluster: float,
    max_groups: int | None,
    stop_on_max_groups: bool,
    overwrite: bool,
    vm_igb: Path | None = None,
    transformed_points: Path | None = None,
    mesh: Path | None = None,
    query_points: Path | None = None,
    node_indices: Path | None = None,
    meshtool: Path | str | None = None,
    igbextract: Path | str | None = None,
    search_radius: float = 50.0,
) -> tuple[str, TCLResult | None]:
    required = TCL_OUTPUTS + (ACTIVATION_OUTPUTS if run_activation else [])
    if _outputs_ready(output_dir, required) and not overwrite:
        print(f"[tcl] Reusing outputs in {output_dir}")
        return "reused", None

    transmembrane_path = _ensure_transmembrane(
        input_path=input_path,
        vm_igb=vm_igb,
        transformed_points=transformed_points,
        mesh=mesh,
        query_points=query_points,
        node_indices=node_indices,
        meshtool=meshtool,
        igbextract=igbextract,
        search_radius=search_radius,
        overwrite=overwrite,
    )

    signals = read_transmembrane_file(transmembrane_path)
    labels = _read_labels_file(labels_file)
    print(
        f"[tcl] Loaded {signals.shape[0]} samples and "
        f"{signals.shape[1]} signals from {transmembrane_path}"
    )

    result = calculate_tcl(
        signals=signals,
        dt=dt,
        reference_index=reference_index,
        derivative_threshold=derivative_threshold,
        min_peak_distance=min_peak_distance,
        signal_labels=labels,
        derivative_per_time=derivative_per_time,
    )
    _write_tcl_outputs(
        result,
        output_dir,
        run_activation=run_activation,
        fraction_window=fraction_window,
        fraction_cluster=fraction_cluster,
        max_groups=max_groups,
        stop_on_max_groups=stop_on_max_groups,
    )

    print(
        f"[tcl] Reference TCL: {result.reference_mean_tcl:.6g} "
        f"± {result.reference_std_tcl:.6g}"
    )
    return "created", result


def _run_ps_detection_stage(
    *,
    points_time: Path,
    vm_igb: Path,
    mesh: Path,
    igbhead: Path | str,
    igbfilament: Path | str,
    cleaned_igb: Path | None,
    threshold: float,
    filament_dt: float,
    overwrite: bool,
    keep_cleaned_igb: bool,
    dry_run: bool,
) -> str:
    result = ensure_phase_singularity_file(
        points_time_path=points_time,
        vm_igb_path=vm_igb,
        mesh_path=mesh,
        igbhead_executable=igbhead,
        igbfilament_executable=igbfilament,
        cleaned_igb_path=cleaned_igb,
        threshold=threshold,
        filament_dt=filament_dt,
        overwrite=overwrite,
        keep_cleaned_igb=keep_cleaned_igb,
        dry_run=dry_run,
    )
    if result.reused_existing:
        print(f"[ps-detect] Reused {result.points_time_path}")
        return "reused"
    if dry_run:
        print(f"[ps-detect] Dry run; expected output: {result.points_time_path}")
        return "dry-run"
    print(f"[ps-detect] Created {result.points_time_path}")
    return "created"


def _build_geodesic_graph(
    mesh_coordinates,
    mesh_elements: Path,
    *,
    cache_size: int,
    edge_chunk_size: int,
) -> MeshGeodesicGraph:
    print(f"[geodesic] Building mesh graph from {mesh_elements}")
    graph = MeshGeodesicGraph.from_elem_file(
        mesh_coordinates,
        mesh_elements,
        cache_size=cache_size,
        edge_chunk_size=edge_chunk_size,
    )
    stats = graph.statistics()
    print(
        f"[geodesic] {stats['n_undirected_edges']} edges; "
        f"build {stats['graph_build_seconds']:.3f} s; "
        f"compact memory {stats['compact_memory_bytes'] / (1024 ** 2):.1f} MiB"
    )
    return graph


def _graph_summary_fields(
    graph: MeshGeodesicGraph | None,
) -> dict[str, int | float | None]:
    """Return stable graph-summary fields for either distance backend."""
    if graph is None:
        return {
            "n_mesh_edges": None,
            "n_mesh_elements": None,
            "graph_build_seconds": None,
            "graph_compact_memory_bytes": None,
            "geodesic_cache_size_limit": None,
            "geodesic_cache_entries": None,
            "geodesic_distance_requests": None,
            "geodesic_dijkstra_runs": None,
            "geodesic_cache_hits": None,
            "geodesic_cache_misses": None,
            "geodesic_visited_vertices": None,
            "geodesic_relaxed_edges": None,
        }

    stats = graph.statistics()
    return {
        "n_mesh_edges": stats["n_undirected_edges"],
        "n_mesh_elements": stats["element_count"],
        "graph_build_seconds": stats["graph_build_seconds"],
        "graph_compact_memory_bytes": stats["compact_memory_bytes"],
        "geodesic_cache_size_limit": stats["cache_size_limit"],
        "geodesic_cache_entries": stats["cache_entries"],
        "geodesic_distance_requests": stats["distance_requests"],
        "geodesic_dijkstra_runs": stats["dijkstra_runs"],
        "geodesic_cache_hits": stats["cache_hits"],
        "geodesic_cache_misses": stats["cache_misses"],
        "geodesic_visited_vertices": stats["visited_vertices"],
        "geodesic_relaxed_edges": stats["relaxed_edges"],
    }


def _write_ps_coordinate_outputs(
    result: PSCoordinateResultType,
    mapping: PSMappingResult,
    graph: MeshGeodesicGraph | None,
    output_dir: Path,
    *,
    distance_mode: str,
    mesh_vertex_count: int,
    radius: float,
    history_steps: int,
    min_segment_time: float,
    gap_factor: float,
    hard_cap_fraction: float,
    reject_single_segment: bool,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for cycle_number, counts in result.cycle_counts.items():
        write_array(
            output_dir / f"rotor_counts_cycle_{cycle_number}.dat",
            counts,
            fmt="%d",
        )
    write_array(output_dir / "rotor_counts_ptst.dat", result.filtered_counts, fmt="%d")
    write_array(
        output_dir / "rotor_counts_normalized_ptst.dat",
        result.normalized_counts,
        fmt="%.6f",
    )
    write_dataframe(output_dir / "rotor_segments_filtered.csv", result.segment_report)
    write_dataframe(output_dir / "rotor_lifetimes.csv", result.lifetime_report)
    write_dataframe(output_dir / "ps_mapping_report.csv", mapping.report)

    summary_row = {
        "distance_method": (
            "bounded_mesh_geodesic"
            if distance_mode == "geodesic"
            else "euclidean"
        ),
        "n_mesh_vertices": mesh_vertex_count,
        "n_ps_timesteps": len(mapping.timesteps),
        "n_input_ps_points": mapping.n_input_points,
        "n_mapped_ps_points": mapping.n_mapped_points,
        "n_dropped_ps_points": mapping.n_dropped_points,
        "maximum_mapping_distance": mapping.maximum_distance,
        "reference_tcl": result.reference_tcl,
        "ps_timestep_dt": result.timestep_dt,
        "n_cycles": len(result.cycle_counts),
        "n_vertices_with_raw_hits": int((result.total_counts > 0).sum()),
        "n_accepted_vertices": int(result.accepted_mask.sum()),
        "radius": radius,
        "history_steps": history_steps,
        "min_segment_time": min_segment_time,
        "gap_factor": gap_factor,
        "hard_cap_fraction": hard_cap_fraction,
        "simulation_duration": result.simulation_duration,
        "maximum_allowed_segment_time": result.maximum_allowed_segment_time,
        "reject_single_segment": reject_single_segment,
    }
    summary_row.update(_graph_summary_fields(graph))
    write_dataframe(
        output_dir / "ps_coords_summary.csv",
        pd.DataFrame([summary_row]),
    )

def _run_ps_coordinates_stage(
    *,
    mesh_coordinates,
    mapping: PSMappingResult,
    graph: MeshGeodesicGraph | None,
    distance_mode: str,
    reference_tcl: float,
    output_dir: Path,
    radius: float,
    history_steps: int,
    min_segment_time: float,
    gap_factor: float,
    hard_cap_fraction: float,
    reject_single_segment: bool,
    overwrite: bool,
) -> tuple[str, PSCoordinateResultType | None]:
    if _outputs_ready(output_dir, PS_COORD_OUTPUTS) and not overwrite:
        print(f"[ps-coords:{distance_mode}] Reusing outputs in {output_dir}")
        return "reused", None

    common_arguments = {
        "mapped_timesteps": mapping.timesteps,
        "mesh_coordinates": mesh_coordinates,
        "reference_tcl": reference_tcl,
        "radius": radius,
        "history_steps": history_steps,
        "min_segment_time": min_segment_time,
        "gap_factor": gap_factor,
        "hard_cap_fraction": hard_cap_fraction,
        "reject_single_segment": reject_single_segment,
    }

    if distance_mode == "geodesic":
        if graph is None:
            raise RuntimeError(
                "A geodesic graph is required for geodesic PS_coords."
            )
        result = analyze_ps_coordinates_geodesic(
            geodesic_graph=graph,
            **common_arguments,
        )
    elif distance_mode == "euclidean":
        result = analyze_ps_coordinates_euclidean(**common_arguments)
    else:
        raise ValueError(f"Unsupported distance mode: {distance_mode}")

    _write_ps_coordinate_outputs(
        result,
        mapping,
        graph,
        output_dir,
        distance_mode=distance_mode,
        mesh_vertex_count=mesh_coordinates.shape[0],
        radius=radius,
        history_steps=history_steps,
        min_segment_time=min_segment_time,
        gap_factor=gap_factor,
        hard_cap_fraction=hard_cap_fraction,
        reject_single_segment=reject_single_segment,
    )
    print(
        f"[ps-coords:{distance_mode}] Accepted vertices: "
        f"{int(result.accepted_mask.sum())}"
    )
    return "created", result

def _write_rotor_outputs(
    result: RotorTrackingResultType,
    mapping: PSMappingResult,
    graph: MeshGeodesicGraph | None,
    output_dir: Path,
    *,
    distance_mode: str,
    mesh_vertex_count: int,
    count_every_hit: bool,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for cycle_number, counts in result.cycle_counts.items():
        write_array(
            output_dir / f"rotor_track_counts_cycle_{cycle_number}.dat",
            counts,
            fmt="%d",
        )
    write_array(output_dir / "rotor_track_counts_ptst.dat", result.full_counts, fmt="%d")
    write_array(
        output_dir / "rotor_track_counts_normalized_ptst.dat",
        result.normalized_counts,
        fmt="%.6f",
    )
    write_array(
        output_dir / "rotor_track_hit_counts_ptst.dat",
        result.hit_counts,
        fmt="%d",
    )
    write_array(
        output_dir / "rotor_track_duration_ptst.dat",
        result.duration_map,
        fmt="%.6f",
    )
    write_dataframe(
        output_dir / "rotor_tracks_counts_report.csv",
        result.track_report,
    )
    write_dataframe(output_dir / "rotor_summary.csv", result.track_summary)
    write_dataframe(output_dir / "rotor_track_segments.csv", result.segment_report)
    write_dataframe(
        output_dir / "rotor_track_drift_splits.csv",
        result.drift_split_report,
    )
    write_dataframe(output_dir / "rotor_track_points.csv", result.track_points)
    write_dataframe(output_dir / "rotor_track_mapping_report.csv", mapping.report)

    summary_row = {
        "matching_method": "global_linear_assignment",
        "distance_method": (
            "bounded_mesh_geodesic"
            if distance_mode == "geodesic"
            else "euclidean"
        ),
        "n_mesh_vertices": mesh_vertex_count,
        "n_ps_timesteps": len(mapping.timesteps),
        "n_input_ps_points": mapping.n_input_points,
        "n_mapped_ps_points": mapping.n_mapped_points,
        "n_dropped_ps_points": mapping.n_dropped_points,
        "maximum_mapping_distance": mapping.maximum_distance,
        "reference_tcl": result.reference_tcl,
        "ps_timestep_dt": result.timestep_dt,
        "n_cycles": len(result.cycle_intervals),
        "nominal_simulation_duration": result.nominal_simulation_duration,
        "maximum_allowed_segment_time": result.maximum_allowed_segment_time,
        "n_all_tracks": len(result.all_tracks),
        "n_accepted_tracks": len(result.accepted_tracks),
        "n_rejected_tracks": len(result.rejected_tracks),
        "n_valid_segments": len(result.valid_segments),
        "n_invalid_segments": len(result.invalid_segments),
        "n_rolling_drift_splits": len(result.drift_split_report),
        "n_vertices_with_track_counts": int((result.full_counts > 0).sum()),
        "n_vertices_with_hit_counts": int((result.hit_counts > 0).sum()),
        "match_radius": result.match_radius,
        "max_gap_steps": result.max_gap_steps,
        "gap_factor": result.gap_factor,
        "min_segment_time": result.min_segment_time,
        "drift_window_time": result.drift_window_time,
        "drift_radius_factor": result.drift_radius_factor,
        "maximum_rolling_drift_distance": (
            result.drift_radius_factor * result.match_radius
        ),
        "hard_cap_fraction": result.max_single_segment_fraction,
        "unique_per_track": not count_every_hit,
    }
    summary_row.update(_graph_summary_fields(graph))
    write_dataframe(
        output_dir / "rotor_tracking_run_summary.csv",
        pd.DataFrame([summary_row]),
    )

def _run_rotor_stage(
    *,
    mesh_coordinates,
    mapping: PSMappingResult,
    graph: MeshGeodesicGraph | None,
    distance_mode: str,
    reference_tcl: float,
    output_dir: Path,
    match_radius: float,
    max_gap_steps: int,
    min_segment_time: float,
    gap_factor: float,
    drift_window_time: float,
    drift_radius_factor: float,
    hard_cap_fraction: float | None,
    count_every_hit: bool,
    overwrite: bool,
) -> tuple[str, RotorTrackingResultType | None]:
    if _outputs_ready(output_dir, ROTOR_TRACK_OUTPUTS) and not overwrite:
        print(f"[rotor-track:{distance_mode}] Reusing outputs in {output_dir}")
        return "reused", None

    common_arguments = {
        "mapped_timesteps": mapping.timesteps,
        "mesh_coordinates": mesh_coordinates,
        "reference_tcl": reference_tcl,
        "match_radius": match_radius,
        "max_gap_steps": max_gap_steps,
        "min_segment_time": min_segment_time,
        "gap_factor": gap_factor,
        "drift_window_time": drift_window_time,
        "drift_radius_factor": drift_radius_factor,
        "max_single_segment_fraction": hard_cap_fraction,
        "unique_per_track": not count_every_hit,
    }

    if distance_mode == "geodesic":
        if graph is None:
            raise RuntimeError(
                "A geodesic graph is required for geodesic rotor tracking."
            )
        result = analyze_rotor_tracking_geodesic(
            geodesic_graph=graph,
            **common_arguments,
        )
    elif distance_mode == "euclidean":
        result = analyze_rotor_tracking_euclidean(**common_arguments)
    else:
        raise ValueError(f"Unsupported distance mode: {distance_mode}")

    _write_rotor_outputs(
        result,
        mapping,
        graph,
        output_dir,
        distance_mode=distance_mode,
        mesh_vertex_count=mesh_coordinates.shape[0],
        count_every_hit=count_every_hit,
    )
    print(
        f"[rotor-track:{distance_mode}] "
        f"Accepted tracks: {len(result.accepted_tracks)}; "
        f"valid segments: {len(result.valid_segments)}"
    )
    return "created", result

def _add_tcl_arguments(parser: argparse.ArgumentParser, *, include_paths: bool = True) -> None:
    if include_paths:
        parser.add_argument("--ascii-vm-input", dest="input",type=Path, required=True, help=("transmembrane_v.dat. " 
                            "The calculated transmembrane voltage ascii file path for your simulation."))
        parser.add_argument("--tcl-output-dir", dest="output_dir",type=Path, default=Path("tcl_results"), help=("Output directory for TCL calculations. "
        "Default: ./tcl_results."
    ),)
    parser.add_argument("--dt", type=float, required=True, help=("Time per signal sample "
                                                                 "Time per signal sample. (Equal to space_dt in CARP... still needs improvement.)"))
    parser.add_argument("--reference-index", type=int, required=True, help=("The signal reference index for mean TCL report"))
    parser.add_argument("--derivative-threshold", type=float, default=5.0)
    parser.add_argument("--min-peak-distance", type=float, default=10.0)
    parser.add_argument("--labels-file", type=Path, default=None, help=("The signal names file path. "
                                                                        "This file path contain the name of the signals for TCL calculation"))
    parser.add_argument("--derivative-per-time", action="store_true")
    parser.add_argument("--fraction-window", type=float, default=0.9, help=("The structured window for signal sequences. "
                                                                            "The TCL ratio used for making the window for signal sequencing."))
    parser.add_argument("--fraction-cluster", type=float, default=0.1, help=("The variation in signal sequences. "
                        "The maximum accepted TCL ratio in signals sequences to be called as the same sequence."))
    parser.add_argument("--max-groups", type=int, default=20, help=("The max number of sequences. "
                                                                    "The maximum number of sequences. After this number arrhythmia will be called AF."))
    parser.add_argument("--no-stop-on-max-groups", action="store_true")
    parser.add_argument("--skip-activation", action="store_true")
    parser.add_argument("--overwrite", action="store_true")


def _add_extraction_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--vm-igb", type=Path)
    parser.add_argument("--transformed-points", type=Path)
    parser.add_argument("--mesh", type=Path)
    parser.add_argument("--query-points", type=Path)
    parser.add_argument("--node-indices", type=Path)
    parser.add_argument("--meshtool", default=None)
    parser.add_argument("--igbextract", default=None)
    parser.add_argument("--search-radius", type=float, default=50.0)


def _add_graph_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--geodesic-cache-size", type=int, default=200000)
    parser.add_argument("--edge-chunk-size", type=int, default=500000)


def _add_distance_mode_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--distance-mode",
        choices=DISTANCE_MODES,
        default="geodesic",
        help=(
            "Distance backend used by the analysis. Default: geodesic."
        ),
    )


def _add_mapping_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--max-mapping-distance", type=float, default=None)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="psfilter",
        description="TCL, phase-singularity, and rotor-tracking analysis.",
    )
    parser.add_argument("--version", action="version", version="psfilter 0.1.0")
    subparsers = parser.add_subparsers(dest="command", required=True)

    tcl_parser = subparsers.add_parser("tcl", help="Calculate TCL and activation patterns.",
                                       description=(
                                        "Calculate tachycardia cycle length and activation patterns "
                                        "from extracted transmembrane-voltage signals."), 
                                       epilog=("If the transmembrane ascii file is missing you can use this part as a way to calculate it. "
                                         "You either need a nodes number file (--node-indices) containing the element number in each line [first line is the number of nodes], "
                                         "or the pts file of the nodes (--transformed-points) containg the X Y Z of each node[first line is the number of nodes]. "
                                         "from this pts file a file is made (--query-points) that is used with meshtool to find the closest indices in a specific radius (--search-radius)."),
                                         )
    _add_tcl_arguments(tcl_parser)
    _add_extraction_arguments(tcl_parser)
    tcl_parser.set_defaults(handler=_command_tcl)

    detect_parser = subparsers.add_parser("ps-detect", help="Create a .pts_t file.",
                                          description=(
                                        "Calculate PSs using igbfilament "
                                        "from transmembrane-voltage igb file and the mesh."), 
                                       epilog=("If the pts_t file is missing you can make it using this mode. "
                                         "You need installed igbfilament path (--igbfilament) and the mesh path. "
                                         "You can change the voltage threshold (threshold) and the timestep (filament-dt) for calculation. "
                                         "In the prcoess a clean igb file is used by igbhead -j and it's used for the calculation which would be deleted after the calculations."),
                                         )
    detect_parser.add_argument("--vm-igb", type=Path, required=True)
    detect_parser.add_argument("--mesh", type=Path, required=True)
    detect_parser.add_argument("--output", type=Path, required=True)
    detect_parser.add_argument("--igbhead", default="igbhead")
    detect_parser.add_argument("--igbfilament", default="igbfilament")
    detect_parser.add_argument("--cleaned-igb", type=Path, default=None)
    detect_parser.add_argument("--threshold", type=float, default=-50.0)
    detect_parser.add_argument("--filament-dt", type=float, default=8.0)
    detect_parser.add_argument("--overwrite", action="store_true")
    detect_parser.add_argument("--keep-cleaned-igb", action="store_true")
    detect_parser.add_argument("--dry-run", action="store_true")
    detect_parser.set_defaults(handler=_command_ps_detect)

    coords_parser = subparsers.add_parser(
        "ps-coords",
        help="Run coordinate-based PS filtering.",
        description=(
    "In this stage, stable rotors are identified. "
    "For large meshes, geodesic distance can be computationally expensive. "
    "Existing nonempty results are reused by default. "
    "Use --overwrite only when you want to recalculate them."
),
        epilog=(
            "We need to filter the wrongly calculated PSs, spatially and temporally. "
            "For the stable PSs we suggest using higher hard cap of the maximum time of the stable rotor (--hard-cap-fraction). "
            "Furthermore, minimum time of the segment (--min-segment-time) and the gap between the calculated ps (--gap-factor) can be increased too. "
            "Also we suggest turn on the single segment long ones.(--allow-single-segment) "
        )
    )
    coords_parser.add_argument("--mesh-points", type=Path, required=True)
    coords_parser.add_argument("--mesh-elements", type=Path, default=None)
    coords_parser.add_argument("--points-time", type=Path, required=True)
    coords_parser.add_argument("--tcl-summary", type=Path, required=True)
    coords_parser.add_argument("--output-dir", type=Path, required=True)
    coords_parser.add_argument("--stable-ps-radius", type=float, default=2000.0)
    coords_parser.add_argument("--history-steps", type=int, default=5)
    coords_parser.add_argument("--min-segment-time", type=float, default=120.0)
    coords_parser.add_argument("--gap-factor", type=float, default=3.0)
    coords_parser.add_argument("--hard-cap-fraction", type=float, default=0.90)
    coords_parser.add_argument("--allow-single-segment", action="store_true")
    coords_parser.add_argument("--overwrite", action="store_true")
    _add_mapping_argument(coords_parser)
    _add_distance_mode_argument(coords_parser)
    _add_graph_arguments(coords_parser)
    coords_parser.set_defaults(handler=_command_ps_coords)

    rotor_parser = subparsers.add_parser(
        "rotor-track",
        help="Run globally assigned rotor tracking.",
        description=(
    "In this stage, rotor trajectories and drift are identified. "
    "For large meshes, geodesic distance can be computationally expensive. "
    "Existing nonempty results are reused by default. "
    "Use --overwrite only when you want to recalculate them."
),
        epilog=(
            "We need to filter the wrongly calculated PSs, spatially and temporally but follow the drifts. "
            "For the drift calculation we suggest using lower hard cap of the maximum time of the stable rotor (--hard-cap-fraction). "
            "However, minimum time of the segment (--min-segment-time) should be less. "
            "To make track linking more permissive, increase --track-radius " 
            "or --drift-window-time. "
            "Also you can control spatial movement with --drift-radius-factor."
        )
    )
    rotor_parser.add_argument("--mesh-points", type=Path, required=True)
    rotor_parser.add_argument("--mesh-elements", type=Path, default=None)
    rotor_parser.add_argument("--points-time", type=Path, required=True)
    rotor_parser.add_argument("--tcl-summary", type=Path, required=True)
    rotor_parser.add_argument("--output-dir", type=Path, required=True)
    rotor_parser.add_argument("--track-radius", type=float, default=2000.0)
    rotor_parser.add_argument("--max-gap-steps", type=int, default=5)
    rotor_parser.add_argument("--min-segment-time", type=float, default=120.0)
    rotor_parser.add_argument("--gap-factor", type=float, default=3.0)
    rotor_parser.add_argument(
        "--drift-window-time",
        "--drift-min-time",
        dest="drift_window_time",
        type=float,
        default=120.0,
    )
    rotor_parser.add_argument("--drift-radius-factor", type=float, default=2.0)
    rotor_parser.add_argument("--hard-cap-fraction", type=float, default=0.90)
    rotor_parser.add_argument("--count-every-hit", action="store_true")
    rotor_parser.add_argument("--overwrite", action="store_true")
    _add_mapping_argument(rotor_parser)
    _add_distance_mode_argument(rotor_parser)
    _add_graph_arguments(rotor_parser)
    rotor_parser.set_defaults(handler=_command_rotor_track)

    all_parser = subparsers.add_parser(
        "all",
        help="Run the complete analysis workflow. We don't suggest using this at first use.",
    )
    all_parser.add_argument("--work-dir", type=Path, default=Path("."))
    all_parser.add_argument(
        "--mesh",
        type=Path,
        required=True,
        help="Mesh basename, normally without .pts or .elem.",
    )
    all_parser.add_argument("--ascii-vm-input", dest="input", type=Path, help=("transmembrane_v.dat. " 
                            "The calculated transmembrane voltage ascii file path for your simulation."))
    all_parser.add_argument("--vm-igb", type=Path, default=None)
    all_parser.add_argument("--mesh-points", type=Path, default=None)
    all_parser.add_argument("--mesh-elements", type=Path, default=None)
    all_parser.add_argument("--points-time", type=Path, default=None)
    all_parser.add_argument("--tcl-output-dir", type=Path, default=None, help=(
        "Output directory for TCL calculations. "
        "Default: <work-dir>/tcl_results."
    ),)
    all_parser.add_argument("--ps-output-dir", type=Path, default=None)
    all_parser.add_argument("--transformed-points", type=Path, default=None)
    all_parser.add_argument("--query-points", type=Path, default=None)
    all_parser.add_argument("--node-indices", type=Path, default=None)
    all_parser.add_argument("--meshtool", default=None)
    all_parser.add_argument("--igbextract", default=None)
    all_parser.add_argument("--search-radius", type=float, default=50.0)

    all_parser.add_argument("--dt", type=float, required=True, help=("Time per signal sample. "
                                                                 "Time per signal sample. (Equal to space_dt in CARP... still needs improvement.)"))
    all_parser.add_argument("--reference-index", type=int, required=True, help=("The signal reference index for mean TCL report"))
    all_parser.add_argument("--derivative-threshold", type=float, default=5.0)
    all_parser.add_argument("--min-peak-distance", type=float, default=10.0)
    all_parser.add_argument("--labels-file", type=Path, default=None, help=("The signal names file path. "
                            "This file path contain the name of the signals for TCL calculation"))
    all_parser.add_argument("--derivative-per-time", action="store_true")
    all_parser.add_argument("--fraction-window", type=float, default=0.9, help=("The structured window for signal sequences. "
                                                                                "The TCL ratio used for making the window for signal sequencing."))
    all_parser.add_argument("--fraction-cluster", type=float, default=0.1, help=("The variation in signal sequences. "
                            "The maximum accepted TCL ratio in signals sequences to be called as the same sequence."))
    all_parser.add_argument("--max-groups", type=int, default=20, help=("The max number of sequences. "
                            "The maximum number of sequences. After this number arrhythmia will be called AF."))
    all_parser.add_argument("--no-stop-on-max-groups", action="store_true")
    all_parser.add_argument("--skip-activation", action="store_true")

    all_parser.add_argument("--igbhead", default="igbhead")
    all_parser.add_argument("--igbfilament", default="igbfilament")
    all_parser.add_argument("--cleaned-igb", type=Path, default=None)
    all_parser.add_argument("--threshold", type=float, default=-50.0)
    all_parser.add_argument("--filament-dt", type=float, default=8.0)
    all_parser.add_argument("--keep-cleaned-igb", action="store_true")
    all_parser.add_argument("--skip-ps-detection", action="store_true")

    all_parser.add_argument("--stable-ps-radius", type=float, default=2000.0)
    all_parser.add_argument("--history-steps", type=int, default=5)
    all_parser.add_argument("--ps-min-segment-time", type=float, default=120.0)
    all_parser.add_argument("--ps-gap-factor", type=float, default=3.0)
    all_parser.add_argument("--ps-hard-cap-fraction", type=float, default=0.90)
    all_parser.add_argument("--allow-single-segment", action="store_true")
    all_parser.add_argument("--skip-ps-coords", action="store_true")
    all_parser.add_argument(
        "--ps-distance-mode",
        choices=DISTANCE_MODES,
        default="geodesic",
        help="Distance backend for PS_coords. Default: geodesic.",
    )
    all_parser.add_argument(
        "--ps-coords-output-dir",
        type=Path,
        default=None,
        help=(
            "Output directory for PS_coords. By default, a mode-specific "
            "subdirectory is created under --ps-output-dir."
        ),
    )

    all_parser.add_argument("--track-radius", type=float, default=2000.0)
    all_parser.add_argument("--max-gap-steps", type=int, default=5)
    all_parser.add_argument("--track-min-segment-time", type=float, default=120.0)
    all_parser.add_argument("--track-gap-factor", type=float, default=3.0)
    all_parser.add_argument("--drift-window-time", type=float, default=120.0)
    all_parser.add_argument("--drift-radius-factor", type=float, default=2.0)
    all_parser.add_argument("--track-hard-cap-fraction", type=float, default=0.90)
    all_parser.add_argument("--count-every-hit", action="store_true")
    all_parser.add_argument("--skip-rotor-track", action="store_true")
    all_parser.add_argument(
        "--track-distance-mode",
        choices=DISTANCE_MODES,
        default="geodesic",
        help="Distance backend for rotor tracking. Default: geodesic.",
    )
    all_parser.add_argument(
        "--rotor-output-dir",
        type=Path,
        default=None,
        help=(
            "Output directory for rotor tracking. By default, a mode-specific "
            "subdirectory is created under --ps-output-dir."
        ),
    )

    all_parser.add_argument("--max-mapping-distance", type=float, default=None)
    _add_graph_arguments(all_parser)
    all_parser.add_argument("--overwrite", action="store_true")
    all_parser.set_defaults(handler=_command_all)

    return parser


def _command_tcl(args: argparse.Namespace) -> int:
    _run_tcl_stage(
        input_path=args.input,
        output_dir=args.output_dir,
        dt=args.dt,
        reference_index=args.reference_index,
        derivative_threshold=args.derivative_threshold,
        min_peak_distance=args.min_peak_distance,
        labels_file=args.labels_file,
        derivative_per_time=args.derivative_per_time,
        run_activation=not args.skip_activation,
        fraction_window=args.fraction_window,
        fraction_cluster=args.fraction_cluster,
        max_groups=args.max_groups,
        stop_on_max_groups=not args.no_stop_on_max_groups,
        overwrite=args.overwrite,
        vm_igb=args.vm_igb,
        transformed_points=args.transformed_points,
        mesh=args.mesh,
        query_points=args.query_points,
        node_indices=args.node_indices,
        meshtool=args.meshtool,
        igbextract=args.igbextract,
        search_radius=args.search_radius,
    )
    return 0


def _command_ps_detect(args: argparse.Namespace) -> int:
    _run_ps_detection_stage(
        points_time=args.output,
        vm_igb=args.vm_igb,
        mesh=args.mesh,
        igbhead=args.igbhead,
        igbfilament=args.igbfilament,
        cleaned_igb=args.cleaned_igb,
        threshold=args.threshold,
        filament_dt=args.filament_dt,
        overwrite=args.overwrite,
        keep_cleaned_igb=args.keep_cleaned_igb,
        dry_run=args.dry_run,
    )
    return 0


def _load_ps_inputs(
    *,
    mesh_points: Path,
    mesh_elements: Path | None,
    points_time: Path,
    tcl_summary: Path,
    max_mapping_distance: float | None,
    geodesic_cache_size: int,
    edge_chunk_size: int,
    need_geodesic: bool,
):
    mesh_coordinates = read_pts_file(mesh_points)

    graph = None
    if need_geodesic:
        if mesh_elements is None:
            mesh_elements = mesh_points.with_suffix(".elem")
        graph = _build_geodesic_graph(
            mesh_coordinates,
            mesh_elements,
            cache_size=geodesic_cache_size,
            edge_chunk_size=edge_chunk_size,
        )

    ps_timesteps = read_pts_t_file(points_time)
    mapping = map_ps_coordinates_to_vertices(
        ps_timesteps=ps_timesteps,
        mesh_coordinates=mesh_coordinates,
        max_distance=max_mapping_distance,
    )
    reference_tcl = read_mean_tcl(tcl_summary)
    return mesh_coordinates, graph, mapping, reference_tcl

def _command_ps_coords(args: argparse.Namespace) -> int:
    mesh, graph, mapping, reference_tcl = _load_ps_inputs(
        mesh_points=args.mesh_points,
        mesh_elements=args.mesh_elements,
        points_time=args.points_time,
        tcl_summary=args.tcl_summary,
        max_mapping_distance=args.max_mapping_distance,
        geodesic_cache_size=args.geodesic_cache_size,
        edge_chunk_size=args.edge_chunk_size,
        need_geodesic=args.distance_mode == "geodesic",
    )
    _run_ps_coordinates_stage(
        mesh_coordinates=mesh,
        mapping=mapping,
        graph=graph,
        distance_mode=args.distance_mode,
        reference_tcl=reference_tcl,
        output_dir=args.output_dir,
        radius=args.stable_ps_radius,
        history_steps=args.history_steps,
        min_segment_time=args.min_segment_time,
        gap_factor=args.gap_factor,
        hard_cap_fraction=args.hard_cap_fraction,
        reject_single_segment=not args.allow_single_segment,
        overwrite=args.overwrite,
    )
    return 0

def _command_rotor_track(args: argparse.Namespace) -> int:
    mesh, graph, mapping, reference_tcl = _load_ps_inputs(
        mesh_points=args.mesh_points,
        mesh_elements=args.mesh_elements,
        points_time=args.points_time,
        tcl_summary=args.tcl_summary,
        max_mapping_distance=args.max_mapping_distance,
        geodesic_cache_size=args.geodesic_cache_size,
        edge_chunk_size=args.edge_chunk_size,
        need_geodesic=args.distance_mode == "geodesic",
    )
    hard_cap = None if args.hard_cap_fraction < 0 else args.hard_cap_fraction
    _run_rotor_stage(
        mesh_coordinates=mesh,
        mapping=mapping,
        graph=graph,
        distance_mode=args.distance_mode,
        reference_tcl=reference_tcl,
        output_dir=args.output_dir,
        match_radius=args.track_radius,
        max_gap_steps=args.max_gap_steps,
        min_segment_time=args.min_segment_time,
        gap_factor=args.gap_factor,
        drift_window_time=args.drift_window_time,
        drift_radius_factor=args.drift_radius_factor,
        hard_cap_fraction=hard_cap,
        count_every_hit=args.count_every_hit,
        overwrite=args.overwrite,
    )
    return 0

def _command_all(args: argparse.Namespace) -> int:
    paths = _resolve_pipeline_paths(args)
    paths.tcl_output_dir.mkdir(parents=True, exist_ok=True)
    paths.ps_output_dir.mkdir(parents=True, exist_ok=True)

    ps_coords_output_dir = _resolve_relative(
        args.ps_coords_output_dir,
        paths.work_dir,
    )
    if ps_coords_output_dir is None:
        ps_coords_output_dir = (
            paths.ps_output_dir / f"ps_coords_{args.ps_distance_mode}"
        )

    rotor_output_dir = _resolve_relative(
        args.rotor_output_dir,
        paths.work_dir,
    )
    if rotor_output_dir is None:
        rotor_output_dir = (
            paths.ps_output_dir / f"rotor_track_{args.track_distance_mode}"
        )

    stage_rows: list[dict[str, str]] = []

    tcl_status, tcl_result = _run_tcl_stage(
        input_path=paths.transmembrane_file,
        output_dir=paths.tcl_output_dir,
        dt=args.dt,
        reference_index=args.reference_index,
        derivative_threshold=args.derivative_threshold,
        min_peak_distance=args.min_peak_distance,
        labels_file=_resolve_relative(args.labels_file, paths.work_dir),
        derivative_per_time=args.derivative_per_time,
        run_activation=not args.skip_activation,
        fraction_window=args.fraction_window,
        fraction_cluster=args.fraction_cluster,
        max_groups=args.max_groups,
        stop_on_max_groups=not args.no_stop_on_max_groups,
        overwrite=args.overwrite,
        vm_igb=paths.vm_igb,
        transformed_points=paths.transformed_points,
        mesh=paths.mesh,
        query_points=paths.query_points,
        node_indices=paths.node_indices,
        meshtool=args.meshtool,
        igbextract=args.igbextract,
        search_radius=args.search_radius,
    )
    stage_rows.append({"stage": "tcl", "status": tcl_status})

    tcl_summary_path = paths.tcl_output_dir / "global_tcl_AP.txt"
    reference_tcl = (
        tcl_result.reference_mean_tcl
        if tcl_result is not None
        else read_mean_tcl(tcl_summary_path)
    )

    if args.skip_ps_detection:
        if not _nonempty_file(paths.points_time):
            raise FileNotFoundError(
                "--skip-ps-detection was used, but the .pts_t file does not exist: "
                f"{paths.points_time}"
            )
        detection_status = "skipped-existing"
    else:
        detection_status = _run_ps_detection_stage(
            points_time=paths.points_time,
            vm_igb=paths.vm_igb,
            mesh=paths.mesh,
            igbhead=args.igbhead,
            igbfilament=args.igbfilament,
            cleaned_igb=_resolve_relative(args.cleaned_igb, paths.work_dir),
            threshold=args.threshold,
            filament_dt=args.filament_dt,
            overwrite=args.overwrite,
            keep_cleaned_igb=args.keep_cleaned_igb,
            dry_run=False,
        )
    stage_rows.append({"stage": "ps-detect", "status": detection_status})

    need_ps_coords = not args.skip_ps_coords and (
        args.overwrite or not _outputs_ready(
            ps_coords_output_dir,
            PS_COORD_OUTPUTS,
        )
    )
    need_rotor = not args.skip_rotor_track and (
        args.overwrite or not _outputs_ready(
            rotor_output_dir,
            ROTOR_TRACK_OUTPUTS,
        )
    )

    if need_ps_coords or need_rotor:
        mesh_coordinates = read_pts_file(paths.mesh_points)
        ps_timesteps = read_pts_t_file(paths.points_time)
        mapping = map_ps_coordinates_to_vertices(
            ps_timesteps=ps_timesteps,
            mesh_coordinates=mesh_coordinates,
            max_distance=args.max_mapping_distance,
        )

        need_geodesic_graph = (
            need_ps_coords and args.ps_distance_mode == "geodesic"
        ) or (
            need_rotor and args.track_distance_mode == "geodesic"
        )

        graph = (
            _build_geodesic_graph(
                mesh_coordinates,
                paths.mesh_elements,
                cache_size=args.geodesic_cache_size,
                edge_chunk_size=args.edge_chunk_size,
            )
            if need_geodesic_graph
            else None
        )
    else:
        mesh_coordinates = None
        graph = None
        mapping = None

    if args.skip_ps_coords:
        ps_coords_status = "skipped"
    elif need_ps_coords:
        assert mesh_coordinates is not None and mapping is not None
        ps_coords_status, _ = _run_ps_coordinates_stage(
            mesh_coordinates=mesh_coordinates,
            mapping=mapping,
            graph=graph,
            distance_mode=args.ps_distance_mode,
            reference_tcl=reference_tcl,
            output_dir=ps_coords_output_dir,
            radius=args.stable_ps_radius,
            history_steps=args.history_steps,
            min_segment_time=args.ps_min_segment_time,
            gap_factor=args.ps_gap_factor,
            hard_cap_fraction=args.ps_hard_cap_fraction,
            reject_single_segment=not args.allow_single_segment,
            overwrite=args.overwrite,
        )
    else:
        ps_coords_status = "reused"
    stage_rows.append(
        {
            "stage": f"ps-coords-{args.ps_distance_mode}",
            "status": ps_coords_status,
        }
    )

    if args.skip_rotor_track:
        rotor_status = "skipped"
    elif need_rotor:
        assert mesh_coordinates is not None and mapping is not None
        track_cap = (
            None
            if args.track_hard_cap_fraction < 0
            else args.track_hard_cap_fraction
        )
        rotor_status, _ = _run_rotor_stage(
            mesh_coordinates=mesh_coordinates,
            mapping=mapping,
            graph=graph,
            distance_mode=args.track_distance_mode,
            reference_tcl=reference_tcl,
            output_dir=rotor_output_dir,
            match_radius=args.track_radius,
            max_gap_steps=args.max_gap_steps,
            min_segment_time=args.track_min_segment_time,
            gap_factor=args.track_gap_factor,
            drift_window_time=args.drift_window_time,
            drift_radius_factor=args.drift_radius_factor,
            hard_cap_fraction=track_cap,
            count_every_hit=args.count_every_hit,
            overwrite=args.overwrite,
        )
    else:
        rotor_status = "reused"
    stage_rows.append(
        {
            "stage": f"rotor-track-{args.track_distance_mode}",
            "status": rotor_status,
        }
    )

    pipeline_summary = pd.DataFrame(stage_rows)
    pipeline_summary["work_dir"] = str(paths.work_dir)
    pipeline_summary["transmembrane_file"] = str(paths.transmembrane_file)
    pipeline_summary["points_time"] = str(paths.points_time)
    pipeline_summary["tcl_output_dir"] = str(paths.tcl_output_dir)
    pipeline_summary["ps_output_dir"] = str(paths.ps_output_dir)
    pipeline_summary["ps_coords_output_dir"] = str(ps_coords_output_dir)
    pipeline_summary["rotor_output_dir"] = str(rotor_output_dir)
    write_dataframe(paths.work_dir / "psfilter_pipeline_summary.csv", pipeline_summary)

    print("\nPipeline complete:")
    for row in stage_rows:
        print(f"  {row['stage']}: {row['status']}")
    return 0

def main(argv: Sequence[str] | None = None) -> int:
    """Run the psfilter command-line interface."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except (FileNotFoundError, ValueError, RuntimeError, TypeError) as error:
        parser.error(str(error))
    return 2
