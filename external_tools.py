from pathlib import Path
from collections.abc import Sequence

from .commands import execute_command


PathLike = str | Path


def run_meshtool_idxlist(
    meshtool_executable: PathLike,
    mesh_path: PathLike,
    coordinate_file: PathLike,
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

    execute_command(command)


def run_igbextract(
    igbextract_executable: PathLike,
    input_igb: PathLike,
    node_indices: Sequence[int],
    output_path: PathLike,
) -> None:
    """Extract selected node signals from an IGB file."""
    node_text = ",".join(str(index) for index in node_indices)

    command = [
        igbextract_executable,
        "-l",
        node_text,
        "-O",
        output_path,
        input_igb,
    ]

    execute_command(command)

def igbhead_run(
    igbhead_executable: PathLike,
    input_igb: PathLike,
    output_path: PathLike,
) -> None:
    """run Jive to clean the IGB file."""

    command = [
        igbhead_executable,
        "-j",
        '-f', output_path,
        input_igb,
    ]

    execute_command(command)

def run_igbfilament(
    igbfilament_executable: PathLike,
    input_igb: PathLike,
    input_mesh: PathLike,
    dt_val: Sequence[int]=8,
    threshold_val: Sequence[int]=-50,    
    output_path: PathLike,
) -> None:
    """Run igbfilament over a mesh and igbfile made from it"""
    
    command = [igbfilament_executable,
           '-t', threshold_val,
           '-d', dt_val,
           '-a', output_path,
           input_mesh,
           input_igb]
    execute_command(command)

