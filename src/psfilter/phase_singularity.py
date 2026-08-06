"""Phase-singularity detection workflow.

This module coordinates the external ``igbhead`` and ``igbfilament`` tools.
It reproduces the ``PS_run`` stage without changing the working directory
and without mixing command execution into the CLI.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .external_tools import run_igbfilament, run_igbhead
from .preprocessing import mesh_input_exists


PathLike = str | Path


@dataclass(frozen=True)
class PhaseSingularityDetectionResult:
    """Paths and status returned by phase-singularity detection."""

    points_time_path: Path
    output_prefix: Path
    cleaned_igb_path: Path
    reused_existing: bool
    cleaned_igb_kept: bool


def _is_nonempty_file(path: Path) -> bool:
    """Return True when ``path`` is an existing nonempty regular file."""
    return path.is_file() and path.stat().st_size > 0


def _output_prefix_from_pts_t(points_time_path: Path) -> Path:
    """Derive the igbfilament ``-a`` prefix from a ``.pts_t`` filename."""
    if points_time_path.suffix != ".pts_t":
        raise ValueError(
            "points_time_path must end in '.pts_t'; "
            f"received: {points_time_path}"
        )

    return points_time_path.with_suffix("")


def ensure_phase_singularity_file(
    *,
    points_time_path: PathLike,
    vm_igb_path: PathLike,
    mesh_path: PathLike,
    igbhead_executable: PathLike = "igbhead",
    igbfilament_executable: PathLike = "igbfilament",
    cleaned_igb_path: PathLike | None = None,
    threshold: float = -50.0,
    filament_dt: float = 8.0,
    overwrite: bool = False,
    keep_cleaned_igb: bool = False,
    dry_run: bool = False,
) -> PhaseSingularityDetectionResult:
    """Ensure that an igbfilament ``.pts_t`` output exists.

    The workflow is:

    1. Run ``igbhead -j`` on ``vm.igb`` to create a cleaned temporary IGB.
    2. Run ``igbfilament`` on the cleaned IGB and the mesh.
    3. Verify that the requested ``.pts_t`` file was created.
    4. Remove the cleaned temporary IGB unless requested otherwise.

    Parameters
    ----------
    points_time_path
        Exact expected ``.pts_t`` output path, for example
        ``PS_results/Reentry_surface_iac.pts_t``.
    vm_igb_path
        Original simulation ``vm.igb``.
    mesh_path
        Mesh basename, normally without ``.pts`` or ``.elem``.
    igbhead_executable
        Path or command name for ``igbhead``.
    igbfilament_executable
        Path or command name for ``igbfilament``.
    cleaned_igb_path
        Optional temporary cleaned IGB path. By default, ``clean.igb`` is
        placed beside ``vm.igb``.
    threshold
        Voltage threshold passed to igbfilament.
    filament_dt
        Temporal interval passed through igbfilament ``-d``.
    overwrite
        Recreate the output even when a nonempty ``.pts_t`` already exists.
    keep_cleaned_igb
        Retain the temporary cleaned IGB after successful processing.
    dry_run
        Print commands without executing them. Output existence is not
        verified in dry-run mode.

    Returns
    -------
    PhaseSingularityDetectionResult
        Output paths and whether an existing result was reused.
    """
    points_time_path = Path(points_time_path)
    vm_igb_path = Path(vm_igb_path)
    mesh_path = Path(mesh_path)

    if cleaned_igb_path is None:
        cleaned_igb_path = vm_igb_path.parent / "clean.igb"
    else:
        cleaned_igb_path = Path(cleaned_igb_path)

    output_prefix = _output_prefix_from_pts_t(points_time_path)

    if _is_nonempty_file(points_time_path) and not overwrite:
        return PhaseSingularityDetectionResult(
            points_time_path=points_time_path,
            output_prefix=output_prefix,
            cleaned_igb_path=cleaned_igb_path,
            reused_existing=True,
            cleaned_igb_kept=_is_nonempty_file(cleaned_igb_path),
        )

    if not vm_igb_path.is_file():
        raise FileNotFoundError(f"Input IGB file not found: {vm_igb_path}")

    if not mesh_input_exists(mesh_path):
        raise FileNotFoundError(
            f"Mesh or mesh basename not found: {mesh_path}"
        )

    threshold = float(threshold)
    filament_dt = float(filament_dt)

    if filament_dt <= 0:
        raise ValueError(
            f"filament_dt must be positive; received {filament_dt}."
        )

    points_time_path.parent.mkdir(parents=True, exist_ok=True)
    cleaned_igb_path.parent.mkdir(parents=True, exist_ok=True)

    if overwrite and points_time_path.exists():
        points_time_path.unlink()

    # Never allow an old clean.igb to be mistaken for a newly generated one.
    if cleaned_igb_path.exists():
        cleaned_igb_path.unlink()

    try:
        run_igbhead(
            igbhead_executable=igbhead_executable,
            input_igb=vm_igb_path,
            output_path=cleaned_igb_path,
            dry_run=dry_run,
        )

        if not dry_run and not _is_nonempty_file(cleaned_igb_path):
            raise RuntimeError(
                "igbhead completed, but the cleaned IGB was not created "
                f"or is empty: {cleaned_igb_path}"
            )

        run_igbfilament(
            igbfilament_executable=igbfilament_executable,
            input_igb=cleaned_igb_path,
            input_mesh=mesh_path,
            output_prefix=output_prefix,
            dt_val=filament_dt,
            threshold_val=threshold,
            dry_run=dry_run,
        )

        if not dry_run and not _is_nonempty_file(points_time_path):
            raise RuntimeError(
                "igbfilament completed, but the expected phase-singularity "
                f"file was not created or is empty: {points_time_path}"
            )

    finally:
        if not keep_cleaned_igb and cleaned_igb_path.exists():
            cleaned_igb_path.unlink()

    return PhaseSingularityDetectionResult(
        points_time_path=points_time_path,
        output_prefix=output_prefix,
        cleaned_igb_path=cleaned_igb_path,
        reused_existing=False,
        cleaned_igb_kept=keep_cleaned_igb and cleaned_igb_path.exists(),
    )
