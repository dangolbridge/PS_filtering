"""Cycle activation-order and recurring-sequence analysis.

This module uses activation peaks already detected by :mod:`psfilter.tcl`.
It does not read or write files. File output should be handled with the
utilities in ``io.py``.

The analysis reproduces the main rules of the original script:

1. Each reference activation starts a cycle window.
2. The first activation from every signal inside that window is selected.
3. Activations are ordered by their delay from the reference activation.
4. Consecutive delays separated by at most a fraction of the reference TCL
   are assigned to the same group.
5. Signals inside one group are treated as unordered when activation
   sequences are compared between cycles.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from .tcl import TCLResult


__all__ = [
    "ActivationResult",
    "analyze_activation_patterns",
    "build_activation_sequences",
    "build_global_activation_order",
    "cluster_ordered_delays",
    "extract_cycle_activation_orders",
    "summarize_activation_sequences",
]


_CYCLE_ACTIVATION_COLUMNS = [
    "cycle_index",
    "cycle_number",
    "cycle_start_sample",
    "cycle_start_time",
    "window_end_time",
    "signal_index",
    "signal_label",
    "activation_sample",
    "activation_time",
    "delay",
    "group",
]

_CYCLE_SUMMARY_COLUMNS = [
    "cycle_index",
    "cycle_number",
    "cycle_start_sample",
    "cycle_start_time",
    "window_end_time",
    "n_signals",
    "n_groups",
    "exceeded_max_groups",
]

_SEQUENCE_COLUMNS = [
    "cycle_index",
    "cycle_number",
    "sequence_name",
    "sequence",
    "sequence_key",
    "n_groups",
    "n_signals",
]

_SEQUENCE_COUNT_COLUMNS = [
    "sequence_name",
    "sequence",
    "n_cycles",
    "fraction_of_cycles",
    "first_cycle_number",
]


@dataclass
class ActivationResult:
    """Container for cycle activation-order and sequence results.

    Attributes
    ----------
    cycle_activations
        Long-form table containing one row per selected signal activation
        in each reference-defined cycle.
    cycle_summary
        One row per analysed cycle.
    global_activation_order
        All cycle activations sorted and clustered by delay.
    activation_sequences
        Sequence assigned to each non-empty cycle.
    sequence_counts
        Number and fraction of cycles assigned to each sequence.
    reference_index
        Zero-based index of the reference signal.
    reference_mean_tcl
        Mean TCL used to define cycle windows and delay clusters.
    window_length
        Cycle-window length.
    fraction_window
        Window length as a fraction of reference TCL.
    fraction_cluster
        Maximum consecutive delay gap, as a fraction of reference TCL,
        for membership in the same group.
    max_groups
        Maximum expected number of groups per cycle, or ``None``.
    n_global_groups
        Number of groups in the global delay ordering.
    stopped_early
        Whether cycle analysis stopped after exceeding ``max_groups``.
    """

    cycle_activations: pd.DataFrame
    cycle_summary: pd.DataFrame
    global_activation_order: pd.DataFrame
    activation_sequences: pd.DataFrame
    sequence_counts: pd.DataFrame

    reference_index: int
    reference_mean_tcl: float
    window_length: float
    fraction_window: float
    fraction_cluster: float
    max_groups: int | None
    n_global_groups: int
    stopped_early: bool


def _validate_positive_finite(value: float, name: str) -> float:
    """Return ``value`` as a positive finite float."""
    value = float(value)

    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a positive finite number; got {value!r}.")

    return value


def _validate_nonnegative_finite(value: float, name: str) -> float:
    """Return ``value`` as a non-negative finite float."""
    value = float(value)

    if not np.isfinite(value) or value < 0:
        raise ValueError(
            f"{name} must be a non-negative finite number; got {value!r}."
        )

    return value


def _validate_max_groups(max_groups: int | None) -> int | None:
    """Validate the optional maximum number of activation groups."""
    if max_groups is None:
        return None

    max_groups = int(max_groups)

    if max_groups < 1:
        raise ValueError(f"max_groups must be at least 1; got {max_groups}.")

    return max_groups


def _prepare_peaks_per_signal(
    peaks_per_signal: Sequence[np.ndarray],
) -> list[np.ndarray]:
    """Validate, sort, and return peak indices for all signals."""
    prepared: list[np.ndarray] = []

    for signal_index, peaks in enumerate(peaks_per_signal):
        raw_peaks = np.asarray(peaks).reshape(-1)

        if raw_peaks.size == 0:
            prepared.append(np.empty(0, dtype=int))
            continue

        if not np.issubdtype(raw_peaks.dtype, np.number):
            raise ValueError(
                f"Peak indices for signal {signal_index} must be numerical."
            )

        numeric_peaks = raw_peaks.astype(float)

        if not np.all(np.isfinite(numeric_peaks)):
            raise ValueError(
                f"Peak indices for signal {signal_index} contain NaN or infinity."
            )

        if not np.all(numeric_peaks == np.floor(numeric_peaks)):
            raise ValueError(
                f"Peak indices for signal {signal_index} must be integers."
            )

        integer_peaks = numeric_peaks.astype(int)

        if np.any(integer_peaks < 0):
            raise ValueError(
                f"Peak indices for signal {signal_index} cannot be negative."
            )

        prepared.append(np.unique(integer_peaks))

    if not prepared:
        raise ValueError("At least one signal is required.")

    return prepared


def _resolve_signal_labels(
    n_signals: int,
    signal_labels: Sequence[str] | None,
) -> list[str]:
    """Return one label for each signal."""
    if signal_labels is None:
        return [f"Sig{index}" for index in range(n_signals)]

    labels = [str(label) for label in signal_labels]

    if len(labels) != n_signals:
        raise ValueError(
            f"Expected {n_signals} signal labels, but received {len(labels)}."
        )

    return labels


def _empty_cycle_activations() -> pd.DataFrame:
    """Return an empty cycle-activation table with stable columns."""
    return pd.DataFrame(columns=_CYCLE_ACTIVATION_COLUMNS)


def _empty_cycle_summary() -> pd.DataFrame:
    """Return an empty cycle-summary table with stable columns."""
    return pd.DataFrame(columns=_CYCLE_SUMMARY_COLUMNS)


def _empty_sequences() -> pd.DataFrame:
    """Return an empty activation-sequence table with stable columns."""
    return pd.DataFrame(columns=_SEQUENCE_COLUMNS)


def _empty_sequence_counts() -> pd.DataFrame:
    """Return an empty sequence-count table with stable columns."""
    return pd.DataFrame(columns=_SEQUENCE_COUNT_COLUMNS)


def cluster_ordered_delays(
    sorted_delays: Sequence[float] | np.ndarray,
    tcl: float,
    fraction: float = 0.1,
    max_groups: int | None = 20,
) -> tuple[np.ndarray, int, bool]:
    """Cluster activations according to consecutive delay gaps.

    The input delays must already be sorted in ascending order. Two
    consecutive activations belong to the same group when their delay gap is
    less than or equal to ``fraction * tcl``.

    Unlike the original implementation, all group labels are calculated even
    when the number of groups exceeds ``max_groups``. The function then flags
    the condition instead of returning a shortened group array.

    Parameters
    ----------
    sorted_delays
        Activation delays sorted in ascending order.
    tcl
        Reference tachycardia cycle length.
    fraction
        Grouping threshold as a fraction of ``tcl``.
    max_groups
        Maximum expected number of groups. Use ``None`` to disable the check.

    Returns
    -------
    tuple
        ``(groups, n_groups, exceeded_max_groups)``.
    """
    tcl = _validate_positive_finite(tcl, "tcl")
    fraction = _validate_nonnegative_finite(fraction, "fraction")
    max_groups = _validate_max_groups(max_groups)

    delays = np.asarray(sorted_delays, dtype=float).reshape(-1)

    if delays.size == 0:
        return np.empty(0, dtype=int), 0, False

    if not np.all(np.isfinite(delays)):
        raise ValueError("Activation delays contain NaN or infinite values.")

    if np.any(np.diff(delays) < 0):
        raise ValueError("sorted_delays must be in ascending order.")

    threshold = fraction * tcl
    groups = np.ones(delays.size, dtype=int)

    current_group = 1
    for index in range(1, delays.size):
        if delays[index] - delays[index - 1] > threshold:
            current_group += 1

        groups[index] = current_group

    n_groups = int(current_group)
    exceeded = max_groups is not None and n_groups > max_groups

    return groups, n_groups, exceeded


def extract_cycle_activation_orders(
    peaks_per_signal: Sequence[np.ndarray],
    dt: float,
    reference_index: int,
    reference_mean_tcl: float | None = None,
    signal_labels: Sequence[str] | None = None,
    fraction_window: float = 0.9,
    fraction_cluster: float = 0.1,
    max_groups: int | None = 20,
    stop_on_max_groups: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, bool]:
    """Extract ordered activations for reference-defined cycles.

    For each reference peak, the function selects the first peak from every
    signal in the half-open interval:

    ``[cycle_start, cycle_start + fraction_window * reference_mean_tcl)``

    The selected activations are sorted by delay and assigned to delay groups.

    Parameters
    ----------
    peaks_per_signal
        Peak sample indices for all signals.
    dt
        Time represented by one sample.
    reference_index
        Zero-based reference-signal index.
    reference_mean_tcl
        Mean TCL used to define the window. When omitted, it is calculated
        from the reference peaks.
    signal_labels
        Optional signal labels.
    fraction_window
        Window length as a fraction of reference TCL.
    fraction_cluster
        Delay-group threshold as a fraction of reference TCL.
    max_groups
        Maximum expected groups in one cycle, or ``None``.
    stop_on_max_groups
        Stop after the first cycle that exceeds ``max_groups``, reproducing
        the stopping intention of the original code.

    Returns
    -------
    tuple
        ``(cycle_activations, cycle_summary, stopped_early)``.
    """
    peaks = _prepare_peaks_per_signal(peaks_per_signal)
    dt = _validate_positive_finite(dt, "dt")
    fraction_window = _validate_positive_finite(
        fraction_window,
        "fraction_window",
    )
    fraction_cluster = _validate_nonnegative_finite(
        fraction_cluster,
        "fraction_cluster",
    )
    max_groups = _validate_max_groups(max_groups)

    reference_index = int(reference_index)

    if not 0 <= reference_index < len(peaks):
        raise IndexError(
            f"reference_index={reference_index} is outside the valid range "
            f"0 to {len(peaks) - 1}."
        )

    reference_peaks = peaks[reference_index]

    if reference_peaks.size < 2:
        raise ValueError(
            "The reference signal must contain at least two peaks for "
            "cycle activation analysis."
        )

    if reference_mean_tcl is None:
        reference_intervals = np.diff(reference_peaks).astype(float) * dt
        reference_mean_tcl = float(np.mean(reference_intervals))
    else:
        reference_mean_tcl = _validate_positive_finite(
            reference_mean_tcl,
            "reference_mean_tcl",
        )

    labels = _resolve_signal_labels(len(peaks), signal_labels)
    window_length = fraction_window * reference_mean_tcl

    activation_rows: list[dict[str, int | float | str]] = []
    summary_rows: list[dict[str, int | float | bool]] = []
    stopped_early = False

    for cycle_index, reference_sample in enumerate(reference_peaks):
        cycle_number = cycle_index + 1
        cycle_start_time = float(reference_sample * dt)
        window_end_time = float(cycle_start_time + window_length)

        selected_rows: list[dict[str, int | float | str]] = []

        for signal_index, signal_peaks in enumerate(peaks):
            if signal_peaks.size == 0:
                continue

            candidate_position = int(
                np.searchsorted(signal_peaks, reference_sample, side="left")
            )

            if candidate_position >= signal_peaks.size:
                continue

            activation_sample = int(signal_peaks[candidate_position])
            activation_time = float(activation_sample * dt)

            if activation_time >= window_end_time:
                continue

            selected_rows.append(
                {
                    "signal_index": signal_index,
                    "signal_label": labels[signal_index],
                    "activation_sample": activation_sample,
                    "activation_time": activation_time,
                    "delay": float(activation_time - cycle_start_time),
                }
            )

        selected_rows.sort(
            key=lambda row: (
                float(row["delay"]),
                int(row["signal_index"]),
            )
        )

        if selected_rows:
            delays = np.asarray(
                [float(row["delay"]) for row in selected_rows],
                dtype=float,
            )
            groups, n_groups, exceeded = cluster_ordered_delays(
                sorted_delays=delays,
                tcl=reference_mean_tcl,
                fraction=fraction_cluster,
                max_groups=max_groups,
            )

            for row, group in zip(selected_rows, groups):
                activation_rows.append(
                    {
                        "cycle_index": cycle_index,
                        "cycle_number": cycle_number,
                        "cycle_start_sample": int(reference_sample),
                        "cycle_start_time": cycle_start_time,
                        "window_end_time": window_end_time,
                        **row,
                        "group": int(group),
                    }
                )
        else:
            n_groups = 0
            exceeded = False

        summary_rows.append(
            {
                "cycle_index": cycle_index,
                "cycle_number": cycle_number,
                "cycle_start_sample": int(reference_sample),
                "cycle_start_time": cycle_start_time,
                "window_end_time": window_end_time,
                "n_signals": len(selected_rows),
                "n_groups": n_groups,
                "exceeded_max_groups": bool(exceeded),
            }
        )

        if exceeded and stop_on_max_groups:
            stopped_early = True
            break

    if activation_rows:
        cycle_activations = pd.DataFrame(
            activation_rows,
            columns=_CYCLE_ACTIVATION_COLUMNS,
        )
    else:
        cycle_activations = _empty_cycle_activations()

    if summary_rows:
        cycle_summary = pd.DataFrame(
            summary_rows,
            columns=_CYCLE_SUMMARY_COLUMNS,
        )
    else:
        cycle_summary = _empty_cycle_summary()

    return cycle_activations, cycle_summary, stopped_early


def build_global_activation_order(
    cycle_activations: pd.DataFrame,
    tcl: float,
    fraction_cluster: float = 0.1,
) -> tuple[pd.DataFrame, int]:
    """Build the legacy global activation-delay ordering.

    All selected activations from all cycles are sorted together by delay and
    clustered using the same consecutive-gap rule. This reproduces the global
    summary produced by the original script while keeping it separate from
    the per-cycle analysis.
    """
    tcl = _validate_positive_finite(tcl, "tcl")
    fraction_cluster = _validate_nonnegative_finite(
        fraction_cluster,
        "fraction_cluster",
    )

    required_columns = {
        "cycle_index",
        "cycle_number",
        "signal_index",
        "signal_label",
        "activation_time",
        "delay",
    }
    missing_columns = required_columns.difference(cycle_activations.columns)

    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(
            "cycle_activations is missing required columns: "
            f"{missing_text}."
        )

    if cycle_activations.empty:
        result = cycle_activations.copy()
        result["global_group"] = pd.Series(dtype=int)
        return result, 0

    global_order = cycle_activations.sort_values(
        by=["delay", "cycle_index", "signal_index"],
        kind="stable",
    ).reset_index(drop=True)

    groups, n_groups, _ = cluster_ordered_delays(
        sorted_delays=global_order["delay"].to_numpy(dtype=float),
        tcl=tcl,
        fraction=fraction_cluster,
        max_groups=None,
    )
    global_order["global_group"] = groups

    return global_order, n_groups


def _display_labels_by_index(
    cycle_activations: pd.DataFrame,
) -> dict[int, str]:
    """Create unambiguous display labels, even when names are duplicated."""
    label_by_index = (
        cycle_activations[["signal_index", "signal_label"]]
        .drop_duplicates(subset=["signal_index"])
        .set_index("signal_index")["signal_label"]
        .astype(str)
        .to_dict()
    )

    label_counts = Counter(label_by_index.values())

    return {
        int(signal_index): (
            label
            if label_counts[label] == 1
            else f"{label}[{int(signal_index)}]"
        )
        for signal_index, label in label_by_index.items()
    }


def build_activation_sequences(
    cycle_activations: pd.DataFrame,
) -> pd.DataFrame:
    """Identify recurring grouped activation sequences across cycles.

    Signal order inside one delay group is ignored. Sequence identity is based
    on signal indices, not labels, so duplicate anatomical labels cannot merge
    distinct signals accidentally.
    """
    required_columns = {
        "cycle_index",
        "cycle_number",
        "signal_index",
        "signal_label",
        "delay",
        "group",
    }
    missing_columns = required_columns.difference(cycle_activations.columns)

    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(
            "cycle_activations is missing required columns: "
            f"{missing_text}."
        )

    if cycle_activations.empty:
        return _empty_sequences()

    display_labels = _display_labels_by_index(cycle_activations)
    sequence_names: dict[tuple[tuple[int, ...], ...], str] = {}
    sequence_rows: list[dict[str, int | str]] = []

    cycle_ids = (
        cycle_activations[["cycle_index", "cycle_number"]]
        .drop_duplicates()
        .sort_values("cycle_index")
    )

    for cycle_record in cycle_ids.itertuples(index=False):
        cycle_index = int(cycle_record.cycle_index)
        cycle_number = int(cycle_record.cycle_number)

        cycle_table = cycle_activations[
            cycle_activations["cycle_index"] == cycle_index
        ].sort_values(
            by=["group", "delay", "signal_index"],
            kind="stable",
        )

        grouped_indices: list[tuple[int, ...]] = []
        grouped_display: list[str] = []

        for _, group_table in cycle_table.groupby("group", sort=True):
            signal_indices = tuple(
                sorted(group_table["signal_index"].astype(int).tolist())
            )
            grouped_indices.append(signal_indices)

            group_labels = sorted(
                display_labels[signal_index]
                for signal_index in signal_indices
            )
            grouped_display.append(",".join(group_labels))

        sequence_key = tuple(grouped_indices)

        if sequence_key not in sequence_names:
            sequence_names[sequence_key] = f"Seq_{len(sequence_names) + 1}"

        sequence_rows.append(
            {
                "cycle_index": cycle_index,
                "cycle_number": cycle_number,
                "sequence_name": sequence_names[sequence_key],
                "sequence": " | ".join(grouped_display),
                "sequence_key": repr(sequence_key),
                "n_groups": len(grouped_indices),
                "n_signals": int(len(cycle_table)),
            }
        )

    return pd.DataFrame(sequence_rows, columns=_SEQUENCE_COLUMNS)


def summarize_activation_sequences(
    activation_sequences: pd.DataFrame,
) -> pd.DataFrame:
    """Count how often each activation sequence occurs."""
    required_columns = {
        "cycle_number",
        "sequence_name",
        "sequence",
    }
    missing_columns = required_columns.difference(activation_sequences.columns)

    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(
            "activation_sequences is missing required columns: "
            f"{missing_text}."
        )

    if activation_sequences.empty:
        return _empty_sequence_counts()

    n_cycles = len(activation_sequences)

    summary = (
        activation_sequences.groupby(
            ["sequence_name", "sequence"],
            sort=False,
            as_index=False,
        )
        .agg(
            n_cycles=("cycle_number", "size"),
            first_cycle_number=("cycle_number", "min"),
        )
    )

    summary["fraction_of_cycles"] = summary["n_cycles"] / n_cycles

    return summary[
        [
            "sequence_name",
            "sequence",
            "n_cycles",
            "fraction_of_cycles",
            "first_cycle_number",
        ]
    ]


def analyze_activation_patterns(
    tcl_result: "TCLResult",
    fraction_window: float = 0.9,
    fraction_cluster: float = 0.1,
    max_groups: int | None = 20,
    stop_on_max_groups: bool = True,
) -> ActivationResult:
    """Run cycle ordering, clustering, and sequence identification.

    Parameters
    ----------
    tcl_result
        Result returned by :func:`psfilter.tcl.calculate_tcl`.
    fraction_window
        Cycle-window length as a fraction of reference TCL.
    fraction_cluster
        Consecutive delay-gap threshold as a fraction of reference TCL.
    max_groups
        Maximum expected number of groups in one cycle.
    stop_on_max_groups
        Stop after a cycle exceeds ``max_groups``.

    Returns
    -------
    ActivationResult
        Structured activation-order and sequence results.
    """
    required_attributes = [
        "peaks_per_signal",
        "dt",
        "reference_index",
        "reference_mean_tcl",
        "signal_labels",
    ]
    missing_attributes = [
        name
        for name in required_attributes
        if not hasattr(tcl_result, name)
    ]

    if missing_attributes:
        missing_text = ", ".join(missing_attributes)
        raise TypeError(
            "tcl_result does not provide the required attributes: "
            f"{missing_text}."
        )

    reference_mean_tcl = _validate_positive_finite(
        tcl_result.reference_mean_tcl,
        "tcl_result.reference_mean_tcl",
    )

    cycle_activations, cycle_summary, stopped_early = (
        extract_cycle_activation_orders(
            peaks_per_signal=tcl_result.peaks_per_signal,
            dt=tcl_result.dt,
            reference_index=tcl_result.reference_index,
            reference_mean_tcl=reference_mean_tcl,
            signal_labels=tcl_result.signal_labels,
            fraction_window=fraction_window,
            fraction_cluster=fraction_cluster,
            max_groups=max_groups,
            stop_on_max_groups=stop_on_max_groups,
        )
    )

    global_activation_order, n_global_groups = (
        build_global_activation_order(
            cycle_activations=cycle_activations,
            tcl=reference_mean_tcl,
            fraction_cluster=fraction_cluster,
        )
    )

    activation_sequences = build_activation_sequences(cycle_activations)
    sequence_counts = summarize_activation_sequences(activation_sequences)

    return ActivationResult(
        cycle_activations=cycle_activations,
        cycle_summary=cycle_summary,
        global_activation_order=global_activation_order,
        activation_sequences=activation_sequences,
        sequence_counts=sequence_counts,
        reference_index=int(tcl_result.reference_index),
        reference_mean_tcl=reference_mean_tcl,
        window_length=float(fraction_window * reference_mean_tcl),
        fraction_window=float(fraction_window),
        fraction_cluster=float(fraction_cluster),
        max_groups=_validate_max_groups(max_groups),
        n_global_groups=n_global_groups,
        stopped_early=stopped_early,
    )
