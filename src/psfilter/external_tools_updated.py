"""Wrappers around external openCARP and meshtool command-line tools."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from .commands import execute_command


PathLike = str | Path


def run_meshtool_idxlist(
    meshtool_executable: PathLike,
    mesh_path: PathLike,
    coordinate_file: PathLike,
    *,
    dry_run: bool = False,
) -> None:
    """Use meshtool to map coordinates to mesh-node indices."""
    command = [
        meshtool_executable,
        "query",
        "idxlist",
        "-msh",
        mesh_path,
        "-coord",
        coordinate_file,
    ]

    execute_command(command, dry_run=dry_run)


def run_igbextract(
    igbextract_executable: PathLike,
    input_igb: PathLike,
    node_indices: Sequence[int],
    output_path: PathLike,
    *,
    dry_run: bool = False,
) -> None:
    """Extract selected node signals from an IGB file."""
    node_text = ",".join(str(index) for index in node_indices)

    if not node_text:
        raise ValueError("At least one node index is required for igbextract.")

    command = [
        igbextract_executable,
        "-l",
        node_text,
        "-O",
        output_path,
        input_igb,
    ]

    execute_command(command, dry_run=dry_run)


def run_igbhead(
    igbhead_executable: PathLike,
    input_igb: PathLike,
    output_path: PathLike,
    *,
    dry_run: bool = False,
) -> None:
    """Run igbhead jump correction and write a cleaned IGB file."""
    command = [
        igbhead_executable,
        "-j",
        "-f",
        output_path,
        input_igb,
    ]

    execute_command(command, dry_run=dry_run)


def run_igbfilament(
    igbfilament_executable: PathLike,
    input_igb: PathLike,
    input_mesh: PathLike,
    output_prefix: PathLike,
    *,
    dt_val: float = 8.0,
    threshold_val: float = -50.0,
    dry_run: bool = False,
) -> None:
    """Run igbfilament and write phase-singularity outputs.

    Parameters
    ----------
    igbfilament_executable
        Path or command name for ``igbfilament``.
    input_igb
        Cleaned transmembrane-voltage IGB file.
    input_mesh
        Mesh basename passed to igbfilament.
    output_prefix
        Output prefix passed through ``-a``. The expected coordinate-time
        output is normally ``<output_prefix>.pts_t``.
    dt_val
        Temporal sampling interval used by igbfilament.
    threshold_val
        Voltage threshold used by igbfilament.
    dry_run
        Print the command without executing it.
    """
    command = [
        igbfilament_executable,
        "-t",
        threshold_val,
        "-d",
        dt_val,
        "-a",
        output_prefix,
        input_mesh,
        input_igb,
    ]

    execute_command(command, dry_run=dry_run)
