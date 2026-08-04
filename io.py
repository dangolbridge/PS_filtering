"""Input and output utilities for phase-singularity analysis."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


PathLike = str | Path
PSTimestep = tuple[float, np.ndarray]


def read_pts_file(path: PathLike) -> np.ndarray:
    """
    Read an openCARP/CARP mesh ``.pts`` file.

    Parameters
    ----------
    path
        Path to the ``.pts`` file.

    Returns
    -------
    numpy.ndarray
        Mesh coordinates with shape ``(n_points, 3)``.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the header or coordinate data is invalid.
    """
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(f"Mesh points file not found: {path}")

    with path.open("r", encoding="utf-8") as file:
        header = file.readline().strip()

        try:
            expected_points = int(header)
        except ValueError as error:
            raise ValueError(
                f"The first line of {path} must contain the number of points."
            ) from error

        coordinates = np.loadtxt(file, dtype=float)

    # np.loadtxt returns shape (3,) when the file contains one point.
    coordinates = np.atleast_2d(coordinates)

    if coordinates.ndim != 2 or coordinates.shape[1] < 3:
        raise ValueError(
            f"Expected three coordinate columns in {path}; "
            f"received shape {coordinates.shape}."
        )

    # Ignore any additional columns, if present.
    coordinates = coordinates[:, :3]

    if coordinates.shape[0] != expected_points:
        raise ValueError(
            f"Expected {expected_points} points in {path}, "
            f"but found {coordinates.shape[0]}."
        )

    return coordinates


def read_pts_t_file(path: PathLike) -> list[PSTimestep]:
    """
    Read phase-singularity coordinates from a ``.pts_t`` file.

    The expected structure is::

        #time
        number_of_points
        x y z
        x y z
        ...

    Parameters
    ----------
    path
        Path to the ``.pts_t`` file.

    Returns
    -------
    list of tuple
        Each item contains:

        ``(time, coordinates)``

        where ``coordinates`` has shape ``(n_points, 3)``.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file structure is invalid.
    """
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(
            f"Phase-singularity coordinate file not found: {path}"
        )

    with path.open("r", encoding="utf-8") as file:
        lines = file.readlines()

    timesteps: list[PSTimestep] = []
    line_index = 0

    while line_index < len(lines):
        line = lines[line_index].strip()

        if not line.startswith("#"):
            line_index += 1
            continue

        time_text = line[1:].strip()

        try:
            time_value = float(time_text)
        except ValueError:
            # Ignore non-time comment lines.
            line_index += 1
            continue

        line_index += 1

        if line_index >= len(lines):
            raise ValueError(
                f"Missing point count after time {time_value} in {path}."
            )

        try:
            number_of_points = int(lines[line_index].strip())
        except ValueError as error:
            raise ValueError(
                f"Invalid point count after time {time_value} in {path}."
            ) from error

        if number_of_points < 0:
            raise ValueError(
                f"Negative point count after time {time_value} in {path}."
            )

        line_index += 1
        points: list[list[float]] = []

        for _ in range(number_of_points):
            if line_index >= len(lines):
                raise ValueError(
                    f"Unexpected end of {path} while reading "
                    f"points at time {time_value}."
                )

            fields = lines[line_index].split()

            if len(fields) < 3:
                raise ValueError(
                    f"Expected three coordinates at line "
                    f"{line_index + 1} of {path}."
                )

            try:
                point = [
                    float(fields[0]),
                    float(fields[1]),
                    float(fields[2]),
                ]
            except ValueError as error:
                raise ValueError(
                    f"Invalid coordinates at line "
                    f"{line_index + 1} of {path}."
                ) from error

            points.append(point)
            line_index += 1

        coordinates = np.asarray(points, dtype=float).reshape(-1, 3)
        timesteps.append((time_value, coordinates))

    return timesteps


def read_mean_tcl(
    path: PathLike,
    key: str = "global_mean_rr",
) -> float:
    """
    Read the mean TCL from a key-value text file.

    Expected example::

        global_mean_rr: 257.4
        global_std_rr: 12.8
    """
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(f"TCL summary file not found: {path}")

    with path.open("r", encoding="utf-8") as file:
        for line in file:
            name, separator, value = line.partition(":")

            if separator and name.strip() == key:
                try:
                    return float(value.strip())
                except ValueError as error:
                    raise ValueError(
                        f"Invalid value for {key!r} in {path}: "
                        f"{value.strip()!r}"
                    ) from error

    raise ValueError(f"Could not find {key!r} in {path}.")


def write_tcl_summary(
    path: PathLike,
    mean_rr: float,
    std_rr: float,
) -> None:
    """Write global TCL statistics to a text file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as file:
        file.write(f"global_mean_rr: {mean_rr}\n")
        file.write(f"global_std_rr: {std_rr}\n")


def write_array(
    path: PathLike,
    values: np.ndarray,
    fmt: str = "%.6f",
) -> None:
    """Write a NumPy array to a text file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    np.savetxt(path, values, fmt=fmt)


def write_dataframe(
    path: PathLike,
    dataframe: pd.DataFrame,
) -> None:
    """Write a pandas DataFrame to a CSV file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    dataframe.to_csv(path, index=False)


def required_outputs_exist(
    output_directory: PathLike,
    filenames: Sequence[str],
) -> bool:
    """Return True when all required output files exist."""
    output_directory = Path(output_directory)

    return all(
        (output_directory / filename).is_file()
        for filename in filenames
    )


def read_transmembrane_file(path: PathLike) -> np.ndarray:
    """
    Read an extracted transmembrane-voltage text file.

    The returned array has shape:

        (number_of_time_samples, number_of_signals)

    Parameters
    ----------
    path
        Path to the extracted transmembrane-voltage file, normally
        ``transmembrane_v.dat``.

    Returns
    -------
    numpy.ndarray
        Transmembrane-voltage signals as a two-dimensional array.

    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    ValueError
        If the file is empty or cannot be interpreted as numerical data.
    """
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(
            f"Transmembrane-voltage file not found: {path}"
        )

    try:
        data = np.loadtxt(path, dtype=float, ndmin=2)
    except ValueError as error:
        raise ValueError(
            f"Could not read numerical data from {path}."
        ) from error

    if data.size == 0:
        raise ValueError(
            f"Transmembrane-voltage file is empty: {path}"
        )

    return data
