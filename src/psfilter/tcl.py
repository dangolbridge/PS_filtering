"""Tachycardia cycle length (TCL) calculation utilities.

This module detects activation markers from transmembrane-voltage signals,
calculates cycle-length intervals, and returns structured TCL results.

The expected signal-array orientation is:

    (number_of_time_samples, number_of_signals)

File reading and result writing should remain in ``io.py``. This module performs
only numerical analysis.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.signal import find_peaks


__all__ = [
    "TCLResult",
    "calculate_reference_tcl",
    "calculate_rr_intervals",
    "calculate_tcl",
    "calculate_tcl_statistics",
    "detect_activation_peaks",
    "detect_all_activation_peaks",
]


@dataclass
class TCLResult:
    """Container for the complete output of a TCL calculation.

    Attributes
    ----------
    peaks_per_signal
        Peak sample indices for every signal.
    peaks_table
        One row per detected activation peak.
    per_signal_statistics
        Per-signal TCL statistics.
    pooled_intervals
        All positive TCL intervals pooled across signals.
    pooled_mean_tcl
        Mean of all pooled intervals.
    pooled_std_tcl
        Sample standard deviation of all pooled intervals.
    reference_index
        Index of the signal used as the reference channel.
    reference_intervals
        TCL intervals from the reference signal.
    reference_mean_tcl
        Mean TCL of the reference signal.
    reference_std_tcl
        Sample standard deviation of reference-signal TCL.
    signal_labels
        Label assigned to each signal.
    dt
        Time represented by one sample.
    derivative_threshold
        Minimum derivative peak height used for activation detection.
    min_peak_distance
        Minimum time allowed between detected derivative peaks.
    derivative_per_time
        Whether the derivative was calculated per unit time.
    """

    peaks_per_signal: list[np.ndarray]
    peaks_table: pd.DataFrame
    per_signal_statistics: pd.DataFrame

    pooled_intervals: np.ndarray
    pooled_mean_tcl: float
    pooled_std_tcl: float

    reference_index: int
    reference_intervals: np.ndarray
    reference_mean_tcl: float
    reference_std_tcl: float

    signal_labels: list[str]
    dt: float
    derivative_threshold: float
    min_peak_distance: float
    derivative_per_time: bool


def _validate_positive_finite(value: float, name: str) -> float:
    """Return ``value`` as float after validating that it is positive."""
    value = float(value)

    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a positive finite number; got {value!r}.")

    return value


def _prepare_signal(signal: np.ndarray) -> np.ndarray:
    """Validate and return one signal as a one-dimensional float array."""
    signal = np.asarray(signal, dtype=float)

    if signal.ndim != 1:
        raise ValueError(
            f"A single signal must be one-dimensional; received shape {signal.shape}."
        )

    if signal.size < 3:
        raise ValueError(
            "At least three samples are required to calculate a numerical derivative."
        )

    if not np.all(np.isfinite(signal)):
        raise ValueError("The signal contains NaN or infinite values.")

    return signal


def _prepare_signal_matrix(signals: np.ndarray) -> np.ndarray:
    """Validate and return a samples-by-signals float array."""
    signals = np.asarray(signals, dtype=float)

    if signals.ndim == 1:
        signals = signals.reshape(-1, 1)

    if signals.ndim != 2:
        raise ValueError(
            "Signals must be a two-dimensional array with shape "
            f"(n_samples, n_signals); received shape {signals.shape}."
        )

    n_samples, n_signals = signals.shape

    if n_samples < 3:
        raise ValueError(
            "At least three time samples are required for TCL calculation."
        )

    if n_signals < 1:
        raise ValueError("At least one signal column is required.")

    if not np.all(np.isfinite(signals)):
        raise ValueError("The signal array contains NaN or infinite values.")

    return signals


def _resolve_signal_labels(
    n_signals: int,
    signal_labels: Sequence[str] | None,
) -> list[str]:
    """Return validated labels for all signal columns."""
    if signal_labels is None:
        return [f"Sig{index}" for index in range(n_signals)]

    labels = [str(label) for label in signal_labels]

    if len(labels) != n_signals:
        raise ValueError(
            f"Expected {n_signals} signal labels, but received {len(labels)}."
        )

    return labels


def _sample_mean_std(values: np.ndarray) -> tuple[float, float]:
    """Return mean and sample standard deviation using the old script rules."""
    values = np.asarray(values, dtype=float)

    if values.size == 0:
        return float("nan"), float("nan")

    mean_value = float(np.mean(values))
    std_value = float(np.std(values, ddof=1)) if values.size > 1 else 0.0

    return mean_value, std_value


def detect_activation_peaks(
    signal: np.ndarray,
    dt: float,
    derivative_threshold: float = 5.0,
    min_peak_distance: float = 10.0,
    derivative_per_time: bool = False,
) -> np.ndarray:
    """Detect activation peaks in one transmembrane-voltage signal.

    Activation markers are detected as positive peaks in the numerical
    derivative, reproducing the method used in the original ``TCL_AP`` code.

    Parameters
    ----------
    signal
        One transmembrane-voltage signal.
    dt
        Time represented by one sample.
    derivative_threshold
        Minimum derivative peak height. When ``derivative_per_time=False``,
        the threshold is expressed per sample, matching the original script.
        When ``True``, it is expressed per unit time.
    min_peak_distance
        Minimum time between detected peaks, in the same time unit as ``dt``.
        It is converted internally to a number of samples.
    derivative_per_time
        If ``False``, calculate ``np.gradient(signal)`` to reproduce the old
        behavior. If ``True``, calculate ``np.gradient(signal, dt)``.

    Returns
    -------
    numpy.ndarray
        Detected peak sample indices as integers.
    """
    signal = _prepare_signal(signal)
    dt = _validate_positive_finite(dt, "dt")
    min_peak_distance = _validate_positive_finite(
        min_peak_distance,
        "min_peak_distance",
    )
    derivative_threshold = float(derivative_threshold)

    if not np.isfinite(derivative_threshold):
        raise ValueError(
            "derivative_threshold must be finite; "
            f"got {derivative_threshold!r}."
        )

    distance_samples = max(1, int(round(min_peak_distance / dt)))

    if derivative_per_time:
        derivative = np.gradient(signal, dt)
    else:
        derivative = np.gradient(signal)

    peaks, _ = find_peaks(
        derivative,
        height=derivative_threshold,
        distance=distance_samples,
    )

    return peaks.astype(int, copy=False)


def detect_all_activation_peaks(
    signals: np.ndarray,
    dt: float,
    derivative_threshold: float = 5.0,
    min_peak_distance: float = 10.0,
    derivative_per_time: bool = False,
) -> list[np.ndarray]:
    """Detect activation peaks independently in every signal column."""
    signals = _prepare_signal_matrix(signals)

    return [
        detect_activation_peaks(
            signal=signals[:, signal_index],
            dt=dt,
            derivative_threshold=derivative_threshold,
            min_peak_distance=min_peak_distance,
            derivative_per_time=derivative_per_time,
        )
        for signal_index in range(signals.shape[1])
    ]


def calculate_rr_intervals(
    peak_indices: np.ndarray,
    dt: float,
) -> np.ndarray:
    """Convert consecutive peak-sample differences into TCL intervals.

    Parameters
    ----------
    peak_indices
        Detected activation sample indices.
    dt
        Time represented by one sample.

    Returns
    -------
    numpy.ndarray
        Positive TCL intervals in the same time unit as ``dt``.
    """
    dt = _validate_positive_finite(dt, "dt")
    peak_indices = np.asarray(peak_indices, dtype=int).reshape(-1)

    if peak_indices.size < 2:
        return np.empty(0, dtype=float)

    if np.any(peak_indices < 0):
        raise ValueError("Peak indices cannot be negative.")

    peak_indices = np.sort(peak_indices)
    intervals = np.diff(peak_indices).astype(float) * dt

    return intervals[intervals > 0]


def _build_peaks_table(
    peaks_per_signal: Sequence[np.ndarray],
    dt: float,
    signal_labels: Sequence[str],
) -> pd.DataFrame:
    """Create one table row for every detected activation peak."""
    rows: list[dict[str, int | float | str]] = []

    for signal_index, peaks in enumerate(peaks_per_signal):
        for peak_sample in np.asarray(peaks, dtype=int):
            rows.append(
                {
                    "signal_index": signal_index,
                    "signal_label": signal_labels[signal_index],
                    "peak_sample": int(peak_sample),
                    "peak_time": float(peak_sample * dt),
                }
            )

    return pd.DataFrame(
        rows,
        columns=[
            "signal_index",
            "signal_label",
            "peak_sample",
            "peak_time",
        ],
    )


def calculate_tcl_statistics(
    peaks_per_signal: Sequence[np.ndarray],
    dt: float,
    signal_labels: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Calculate TCL statistics separately for every signal.

    The returned DataFrame contains:

    - ``signal_index``
    - ``signal_label``
    - ``n_peaks``
    - ``n_intervals``
    - ``mean_tcl``
    - ``std_tcl``
    - ``sem_tcl``
    """
    dt = _validate_positive_finite(dt, "dt")
    peaks_per_signal = [
        np.asarray(peaks, dtype=int).reshape(-1)
        for peaks in peaks_per_signal
    ]
    labels = _resolve_signal_labels(len(peaks_per_signal), signal_labels)

    rows: list[dict[str, int | float | str]] = []

    for signal_index, peaks in enumerate(peaks_per_signal):
        intervals = calculate_rr_intervals(peaks, dt)
        mean_tcl, std_tcl = _sample_mean_std(intervals)

        if intervals.size > 1:
            sem_tcl = float(std_tcl / np.sqrt(intervals.size))
        elif intervals.size == 1:
            sem_tcl = 0.0
        else:
            sem_tcl = float("nan")

        rows.append(
            {
                "signal_index": signal_index,
                "signal_label": labels[signal_index],
                "n_peaks": int(peaks.size),
                "n_intervals": int(intervals.size),
                "mean_tcl": mean_tcl,
                "std_tcl": std_tcl,
                "sem_tcl": sem_tcl,
            }
        )

    return pd.DataFrame(rows)


def calculate_reference_tcl(
    peaks_per_signal: Sequence[np.ndarray],
    reference_index: int,
    dt: float,
) -> tuple[float, float, np.ndarray]:
    """Calculate TCL statistics for one reference signal.

    Returns
    -------
    tuple
        ``(mean_tcl, std_tcl, intervals)``

    Raises
    ------
    IndexError
        If ``reference_index`` is outside the signal range.
    ValueError
        If fewer than two peaks exist in the reference signal.
    """
    dt = _validate_positive_finite(dt, "dt")
    reference_index = int(reference_index)

    if not 0 <= reference_index < len(peaks_per_signal):
        raise IndexError(
            f"reference_index={reference_index} is outside the valid range "
            f"0 to {len(peaks_per_signal) - 1}."
        )

    intervals = calculate_rr_intervals(
        peaks_per_signal[reference_index],
        dt,
    )

    if intervals.size == 0:
        n_peaks = np.asarray(peaks_per_signal[reference_index]).size
        raise ValueError(
            "The reference signal must contain at least two distinct peaks; "
            f"signal {reference_index} contains {n_peaks}."
        )

    mean_tcl, std_tcl = _sample_mean_std(intervals)

    return mean_tcl, std_tcl, intervals


def calculate_tcl(
    signals: np.ndarray,
    dt: float,
    reference_index: int,
    derivative_threshold: float = 5.0,
    min_peak_distance: float = 10.0,
    signal_labels: Sequence[str] | None = None,
    derivative_per_time: bool = False,
) -> TCLResult:
    """Run the complete transmembrane-signal TCL calculation.

    This function does not read or write files. Use ``read_transmembrane_file``
    from ``io.py`` before calling it, and use the I/O helpers to save the
    returned tables and statistics.

    Parameters
    ----------
    signals
        Transmembrane-voltage array with shape ``(n_samples, n_signals)``.
    dt
        Time represented by one sample.
    reference_index
        Zero-based index of the reference signal.
    derivative_threshold
        Minimum positive derivative peak height.
    min_peak_distance
        Minimum time between peaks, in the same time unit as ``dt``.
    signal_labels
        Optional label for each signal column.
    derivative_per_time
        Whether to calculate the derivative per unit time.

    Returns
    -------
    TCLResult
        Structured peak and TCL results.
    """
    signals = _prepare_signal_matrix(signals)
    dt = _validate_positive_finite(dt, "dt")
    min_peak_distance = _validate_positive_finite(
        min_peak_distance,
        "min_peak_distance",
    )

    n_signals = signals.shape[1]
    labels = _resolve_signal_labels(n_signals, signal_labels)

    peaks_per_signal = detect_all_activation_peaks(
        signals=signals,
        dt=dt,
        derivative_threshold=derivative_threshold,
        min_peak_distance=min_peak_distance,
        derivative_per_time=derivative_per_time,
    )

    peaks_table = _build_peaks_table(
        peaks_per_signal=peaks_per_signal,
        dt=dt,
        signal_labels=labels,
    )

    per_signal_statistics = calculate_tcl_statistics(
        peaks_per_signal=peaks_per_signal,
        dt=dt,
        signal_labels=labels,
    )

    interval_arrays = [
        calculate_rr_intervals(peaks, dt)
        for peaks in peaks_per_signal
    ]
    nonempty_intervals = [
        intervals
        for intervals in interval_arrays
        if intervals.size > 0
    ]

    if nonempty_intervals:
        pooled_intervals = np.concatenate(nonempty_intervals)
    else:
        pooled_intervals = np.empty(0, dtype=float)

    pooled_mean_tcl, pooled_std_tcl = _sample_mean_std(pooled_intervals)

    (
        reference_mean_tcl,
        reference_std_tcl,
        reference_intervals,
    ) = calculate_reference_tcl(
        peaks_per_signal=peaks_per_signal,
        reference_index=reference_index,
        dt=dt,
    )

    return TCLResult(
        peaks_per_signal=peaks_per_signal,
        peaks_table=peaks_table,
        per_signal_statistics=per_signal_statistics,
        pooled_intervals=pooled_intervals,
        pooled_mean_tcl=pooled_mean_tcl,
        pooled_std_tcl=pooled_std_tcl,
        reference_index=int(reference_index),
        reference_intervals=reference_intervals,
        reference_mean_tcl=reference_mean_tcl,
        reference_std_tcl=reference_std_tcl,
        signal_labels=labels,
        dt=dt,
        derivative_threshold=float(derivative_threshold),
        min_peak_distance=min_peak_distance,
        derivative_per_time=bool(derivative_per_time),
    )
