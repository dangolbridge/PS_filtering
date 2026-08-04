from pathlib import Path

from .external_tools import (
    run_igbextract,
    run_meshtool_idxlist,
)
from .io import (
    read_index_list,
    write_tcl_query_points,
)


PathLike = str | Path


def ensure_transmembrane_file(
    *,
    output_path: PathLike,
    vm_igb_path: PathLike,
    transformed_points_path: PathLike,
    mesh_path: PathLike,
    query_points_path: PathLike,
    node_indices_path: PathLike,
    meshtool_executable: PathLike,
    igbextract_executable: PathLike,
    search_radius: float = 50.0,
    overwrite: bool = False,
) -> Path:
    """
    Ensure that the extracted transmembrane-voltage file exists.

    If necessary, this function:

    1. Creates the coordinate-query file.
    2. Runs meshtool to determine mesh-node indices.
    3. Runs igbextract to create the voltage text file.
    """
    output_path = Path(output_path)
    vm_igb_path = Path(vm_igb_path)
    transformed_points_path = Path(transformed_points_path)
    mesh_path = Path(mesh_path)
    query_points_path = Path(query_points_path)
    node_indices_path = Path(node_indices_path)

    # The requested output already exists.
    if output_path.is_file() and not overwrite:
        return output_path

    if not vm_igb_path.is_file():
        raise FileNotFoundError(
            f"Input IGB file not found: {vm_igb_path}"
        )

    if not mesh_path.exists():
        raise FileNotFoundError(
            f"Mesh not found: {mesh_path}"
        )

    if not node_indices_path.is_file():
        if not transformed_points_path.is_file():
            raise FileNotFoundError(
                "Points used for TCL extraction were not found: "
                f"{transformed_points_path}"
            )

        if not query_points_path.is_file():
            write_tcl_query_points(
                mesh_points_path=transformed_points_path,
                output_path=query_points_path,
                search_radius=search_radius,
            )

        run_meshtool_idxlist(
            meshtool_executable=meshtool_executable,
            mesh_path=mesh_path,
            coordinate_file=query_points_path,
        )

        if not node_indices_path.is_file():
            raise RuntimeError(
                "Meshtool completed, but the expected index-list "
                f"file was not created: {node_indices_path}"
            )

    node_indices = read_index_list(node_indices_path)

    if node_indices.size == 0:
        raise ValueError(
            f"No node indices were found in {node_indices_path}."
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    run_igbextract(
        igbextract_executable=igbextract_executable,
        input_igb=vm_igb_path,
        node_indices=node_indices,
        output_path=output_path,
    )

    if not output_path.is_file():
        raise RuntimeError(
            "igbextract completed, but the expected output "
            f"was not created: {output_path}"
        )

    return output_path
