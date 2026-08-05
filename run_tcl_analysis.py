#!/usr/bin/env python3
"""Run TCL and activation-pattern analysis from transmembrane_v.dat."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from psfilter.activation import analyze_activation_patterns
from psfilter.io import (
    read_transmembrane_file,
    write_dataframe,
    write_tcl_summary,
)
from psfilter.tcl import calculate_tcl
from psfilter.preprocessing import ensure_transmembrane_file

def read_labels_file(path: Path | None) -> list[str] | None:
    """Read one signal label per nonempty line."""
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate tachycardia cycle length and activation patterns "
            "from an extracted transmembrane-voltage text file."
        )
    )

    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to transmembrane_v.dat.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("tcl_results"),
        help="Directory for analysis outputs.",
    )
    parser.add_argument(
        "--dt",
        type=float,
        required=True,
        help="Time represented by one row/sample, normally in ms.",
    )
    parser.add_argument(
        "--reference-index",
        type=int,
        required=True,
        help="Zero-based index of the reference signal.",
    )
    parser.add_argument(
        "--derivative-threshold",
        type=float,
        default=5.0,
        help="Minimum positive derivative peak height. Default: 5.",
    )
    parser.add_argument(
        "--min-peak-distance",
        type=float,
        default=10.0,
        help=(
            "Minimum time between detected peaks, in the same unit as dt. "
            "Default: 10."
        ),
    )
    parser.add_argument(
        "--labels-file",
        type=Path,
        default=None,
        help="Optional text file containing one signal label per line.",
    )
    parser.add_argument(
        "--derivative-per-time",
        action="store_true",
        help=(
            "Calculate dV/dt using dt. Leave unset to reproduce the original "
            "np.gradient(signal) behavior."
        ),
    )
    parser.add_argument(
        "--fraction-window",
        type=float,
        default=0.9,
        help="Activation window as a fraction of reference TCL. Default: 0.9.",
    )
    parser.add_argument(
        "--fraction-cluster",
        type=float,
        default=0.1,
        help="Delay-group threshold as a fraction of TCL. Default: 0.1.",
    )
    parser.add_argument(
        "--max-groups",
        type=int,
        default=20,
        help="Maximum expected activation groups per cycle. Default: 20.",
    )
    parser.add_argument(
        "--no-stop-on-max-groups",
        action="store_true",
        help="Continue analysing cycles when max-groups is exceeded.",
    )
    parser.add_argument(
        "--skip-activation",
        action="store_true",
        help="Calculate TCL only and skip activation-pattern analysis.",
    )
    parser.add_argument("--vm-igb", type=Path)
    parser.add_argument("--transformed-points", type=Path)
    parser.add_argument("--mesh", type=Path)
    parser.add_argument("--query-points", type=Path)
    parser.add_argument("--node-indices", type=Path)
    parser.add_argument("--meshtool", type=Path)
    parser.add_argument("--igbextract", type=Path)

    return parser


def main() -> None:
    args = build_parser().parse_args()

    if not args.input.is_file():
        transmembrane_path = ensure_transmembrane_file(
            output_path=args.input,
            vm_igb_path=args.vm_igb,
            transformed_points_path=args.transformed_points,
            mesh_path=args.mesh,
            query_points_path=args.query_points,
            node_indices_path=args.node_indices,
            meshtool_executable=args.meshtool,
            igbextract_executable=args.igbextract,
        )
    else:
        transmembrane_path = args.input

    signals = read_transmembrane_file(transmembrane_path)
    
    labels = read_labels_file(args.labels_file)

    print(
        f"Loaded {signals.shape[0]} samples and "
        f"{signals.shape[1]} signals from {args.input}"
    )

    tcl_result = calculate_tcl(
        signals=signals,
        dt=args.dt,
        reference_index=args.reference_index,
        derivative_threshold=args.derivative_threshold,
        min_peak_distance=args.min_peak_distance,
        signal_labels=labels,
        derivative_per_time=args.derivative_per_time,
    )

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    write_dataframe(
        output_dir / "tcl_peaks_AP.csv",
        tcl_result.peaks_table,
    )
    write_dataframe(
        output_dir / "tcl_per_signal_AP.csv",
        tcl_result.per_signal_statistics,
    )

    # Preserve compatibility with the old downstream TCL reader.
    write_tcl_summary(
        output_dir / "global_tcl_AP.txt",
        mean_rr=tcl_result.reference_mean_tcl,
        std_rr=tcl_result.reference_std_tcl,
    )

    tcl_summary = pd.DataFrame(
        [
            {
                "reference_index": tcl_result.reference_index,
                "reference_label": (
                    tcl_result.signal_labels[tcl_result.reference_index]
                ),
                "reference_mean_tcl": tcl_result.reference_mean_tcl,
                "reference_std_tcl": tcl_result.reference_std_tcl,
                "reference_n_intervals": len(
                    tcl_result.reference_intervals
                ),
                "pooled_mean_tcl": tcl_result.pooled_mean_tcl,
                "pooled_std_tcl": tcl_result.pooled_std_tcl,
                "pooled_n_intervals": len(tcl_result.pooled_intervals),
                "dt": tcl_result.dt,
                "derivative_threshold": (
                    tcl_result.derivative_threshold
                ),
                "min_peak_distance": tcl_result.min_peak_distance,
                "derivative_per_time": (
                    tcl_result.derivative_per_time
                ),
            }
        ]
    )
    write_dataframe(output_dir / "tcl_summary.csv", tcl_summary)

    print(
        "Reference TCL: "
        f"{tcl_result.reference_mean_tcl:.6g} "
        f"± {tcl_result.reference_std_tcl:.6g}"
    )
    print(
        "Pooled TCL: "
        f"{tcl_result.pooled_mean_tcl:.6g} "
        f"± {tcl_result.pooled_std_tcl:.6g}"
    )

    if args.skip_activation:
        print(f"TCL outputs written to {output_dir}")
        return

    activation_result = analyze_activation_patterns(
        tcl_result=tcl_result,
        fraction_window=args.fraction_window,
        fraction_cluster=args.fraction_cluster,
        max_groups=args.max_groups,
        stop_on_max_groups=not args.no_stop_on_max_groups,
    )

    write_dataframe(
        output_dir / "cycle_activations.csv",
        activation_result.cycle_activations,
    )
    write_dataframe(
        output_dir / "cycle_summary.csv",
        activation_result.cycle_summary,
    )
    write_dataframe(
        output_dir / "global_activation_order.csv",
        activation_result.global_activation_order,
    )
    write_dataframe(
        output_dir / "activation_sequences.csv",
        activation_result.activation_sequences,
    )
    write_dataframe(
        output_dir / "sequence_counts.csv",
        activation_result.sequence_counts,
    )

    print(
        f"Analysed {len(activation_result.cycle_summary)} cycles; "
        f"found {len(activation_result.sequence_counts)} unique sequences."
    )
    print(f"All outputs written to {output_dir}")


if __name__ == "__main__":
    main()
