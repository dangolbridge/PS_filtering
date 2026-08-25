# PSFilter example workflow

This example demonstrates the complete PSFilter workflow from simulation data
to TCL analysis, activation-sequence analysis, phase-singularity (PS)
detection, stable-PS filtering, and rotor-drift tracking.

The example is organized so that small files can be stored directly in the
GitHub repository, while the large `vm.igb` file can be downloaded separately.

The commands below demonstrate every PSFilter CLI protocol:

```text
psfilter tcl
psfilter ps-detect
psfilter ps-coords
psfilter rotor-track
psfilter all
```

The PS-analysis examples also show both Euclidean and geodesic distance modes.

---

## 1. Example directory structure

A suggested directory structure is:

```text
examples/
├── README.md
├── data/
│   ├── example_mesh.pts
│   ├── example_mesh.elem
│   ├── transmembrane_v.dat
│   ├── signal_labels.txt              # optional
│   └── vm.igb                         # downloaded separately
│
├── results/
│   ├── tcl_results/
│   └── PS_results/
│
└── figures/
    ├── 01_input_mesh.png
    ├── 02_tcl_results.png
    ├── 03_signal_sequences.png
    ├── 04_ps_detection.png
    ├── 05_stable_ps_euclidean.png
    ├── 06_stable_ps_geodesic.png
    ├── 07_rotor_tracking_euclidean.png
    ├── 08_rotor_tracking_geodesic.png
    └── 09_complete_pipeline.png
```

Small example files can be committed directly to GitHub:

```text
example_mesh.pts
example_mesh.elem
transmembrane_v.dat
signal_labels.txt
```

The large original `vm.igb` file should be downloaded separately.

Download it from:

```text
[ADD DOWNLOAD LINK FOR vm.igb HERE]
```

and save it as:

```text
examples/data/vm.igb
```

For example:

```bash
cd examples/data
wget "[ADD DOWNLOAD LINK HERE]" -O vm.igb
```

> **Note**
>
> The `vm.igb` file is required for reproducing the PS-detection stage from the
> original simulation voltage data. The already extracted
> `transmembrane_v.dat` file is included so that TCL and activation analysis
> can be tested without downloading the large IGB file.


## 3. Example parameters

The commands below use these example values:

```text
Mesh basename:          example_mesh
TCL signal timestep:    1.0
Reference signal index: 7

PS voltage threshold:   -50 mV
PS timestep:            8 ms

Stable-PS radius:       2000
Track radius:           2000

Minimum segment time:   120 ms
Gap factor:             3
Maximum gap steps:      5

Drift window time:      120 ms
Drift radius factor:    2

Hard-cap fraction:      0.90
```

These values correspond to the current PSFilter defaults where applicable.

The default PS-filtering parameters were mainly tuned for anatomically scaled
cardiac meshes and arrhythmias with TCL values around 200 ms. They may need to
be adjusted for other meshes, sampling intervals, or arrhythmia types.

> **Units**
>
> Spatial parameters such as `--stable-ps-radius` and `--track-radius` use the
> same spatial units as the mesh coordinates.
>
> Time parameters must use units consistent with the input data.

---

# Protocol 1 — TCL and activation-sequence analysis

The `tcl` protocol calculates tachycardia cycle length from the extracted
transmembrane-voltage signals.

Activation-sequence analysis is performed by default as part of the same
command.

Create the output directory:

```bash
mkdir -p examples/results/tcl_results
```

Run:

```bash
psfilter tcl \
    --ascii-vm-input examples/data/transmembrane_v.dat \
    --tcl-output-dir examples/results/tcl_results \
    --dt 1.0 \
    --reference-index 7
```

Main TCL outputs:

```text
examples/results/tcl_results/
├── tcl_peaks_AP.csv
├── tcl_per_signal_AP.csv
├── global_tcl_AP.txt
└── tcl_summary.csv
```

Activation-analysis outputs:

```text
cycle_activations.csv
cycle_summary.csv
global_activation_order.csv
activation_sequences.csv
sequence_counts.csv
```

The file used by later PS analyses is:

```text
examples/results/tcl_results/global_tcl_AP.txt
```

## Optional signal labels

```bash
psfilter tcl \
    --ascii-vm-input examples/data/transmembrane_v.dat \
    --tcl-output-dir examples/results/tcl_results \
    --dt 1.0 \
    --reference-index 7 \
    --labels-file examples/data/signal_labels.txt
```

## TCL only

To calculate TCL without activation-sequence analysis:

```bash
psfilter tcl \
    --ascii-vm-input examples/data/transmembrane_v.dat \
    --tcl-output-dir examples/results/tcl_results \
    --dt 1.0 \
    --reference-index 7 \
    --skip-activation
```

For this complete example, do **not** use `--skip-activation`, because the
signal-sequence analysis is part of the demonstrated workflow.

### Figure placeholder — TCL calculation

> **ADD FIGURE HERE**
>
> Suggested file: `figures/02_tcl_results.png`
>
> Suggested content: example transmembrane-voltage signals, detected
> activations/peaks, and/or calculated TCL.

<!-- Replace with:
![Example TCL calculation](figures/02_tcl_results.png)
-->

### Figure placeholder — signal sequences

> **ADD FIGURE HERE**
>
> Suggested file: `figures/03_signal_sequences.png`
>
> Suggested content: activation-order sequences or sequence-count results.

<!-- Replace with:
![Example activation-sequence analysis](figures/03_signal_sequences.png)
-->

---

# Protocol 2 — Phase-singularity detection

The `ps-detect` protocol calculates PS coordinates from `vm.igb` using
`igbhead` and `igbfilament`.

Create the output directory:

```bash
mkdir -p examples/results/PS_results
```

Run:

```bash
psfilter ps-detect \
    --vm-igb examples/data/vm.igb \
    --mesh examples/data/example_mesh \
    --output examples/results/PS_results/example_mesh.pts_t \
    --threshold -50 \
    --filament-dt 8
```

If the external tools are not in `PATH`:

```bash
psfilter ps-detect \
    --vm-igb examples/data/vm.igb \
    --mesh examples/data/example_mesh \
    --output examples/results/PS_results/example_mesh.pts_t \
    --igbhead /path/to/igbhead \
    --igbfilament /path/to/igbfilament \
    --threshold -50 \
    --filament-dt 8
```

The principal output is:

```text
examples/results/PS_results/example_mesh.pts_t
```

This file contains the raw PS coordinates through time and is used by both
`ps-coords` and `rotor-track`.

## Dry run

```bash
psfilter ps-detect \
    --vm-igb examples/data/vm.igb \
    --mesh examples/data/example_mesh \
    --output examples/results/PS_results/example_mesh.pts_t \
    --threshold -50 \
    --filament-dt 8 \
    --dry-run
```

## Keep the cleaned IGB file

The cleaned intermediate IGB is normally temporary.

To retain it:

```bash
psfilter ps-detect \
    --vm-igb examples/data/vm.igb \
    --mesh examples/data/example_mesh \
    --output examples/results/PS_results/example_mesh.pts_t \
    --threshold -50 \
    --filament-dt 8 \
    --keep-cleaned-igb
```

### Figure placeholder — raw PS detection

> **ADD FIGURE HERE**
>
> Suggested file: `figures/04_ps_detection.png`
>
> Suggested content: raw PS locations calculated by `igbfilament`.

<!-- Replace with:
![Raw PS detection](figures/04_ps_detection.png)
-->

---

# Protocol 3 — Stable PS filtering

The `ps-coords` protocol filters raw PS detections spatially and temporally to
identify persistent/stable PS locations.

Two distance modes are available:

```text
euclidean
geodesic
```

For large meshes, Euclidean distance is generally much faster.

---

## 3A. Stable PS filtering — Euclidean

```bash
mkdir -p examples/results/PS_results/ps_coords_euclidean

psfilter ps-coords \
    --mesh-points examples/data/example_mesh.pts \
    --points-time examples/results/PS_results/example_mesh.pts_t \
    --tcl-summary examples/results/tcl_results/global_tcl_AP.txt \
    --output-dir examples/results/PS_results/ps_coords_euclidean \
    --distance-mode euclidean \
    --stable-ps-radius 2000 \
    --history-steps 5 \
    --min-segment-time 120 \
    --gap-factor 3 \
    --hard-cap-fraction 0.90
```

Main outputs include:

```text
rotor_counts_cycle_<N>.dat
rotor_counts_ptst.dat
rotor_counts_normalized_ptst.dat
rotor_segments_filtered.csv
rotor_lifetimes.csv
ps_mapping_report.csv
ps_coords_summary.csv
```

### Figure placeholder — stable PS, Euclidean

> **ADD FIGURE HERE**
>
> Suggested file: `figures/05_stable_ps_euclidean.png`
>
> Suggested content: stable-PS count map or filtered PS locations using
> Euclidean distance.

<!-- Replace with:
![Stable PS filtering using Euclidean distance](figures/05_stable_ps_euclidean.png)
-->

---

## 3B. Stable PS filtering — geodesic

```bash
mkdir -p examples/results/PS_results/ps_coords_geodesic

psfilter ps-coords \
    --mesh-points examples/data/example_mesh.pts \
    --mesh-elements examples/data/example_mesh.elem \
    --points-time examples/results/PS_results/example_mesh.pts_t \
    --tcl-summary examples/results/tcl_results/global_tcl_AP.txt \
    --output-dir examples/results/PS_results/ps_coords_geodesic \
    --distance-mode geodesic \
    --stable-ps-radius 2000 \
    --history-steps 5 \
    --min-segment-time 120 \
    --gap-factor 3 \
    --hard-cap-fraction 0.90
```

If `--mesh-elements` is omitted, PSFilter attempts to infer the `.elem` file
from the `--mesh-points` basename.

> **Performance note**
>
> Geodesic stable-PS filtering can be computationally expensive on large
> meshes because many graph-distance calculations may be required.

### Figure placeholder — stable PS, geodesic

> **ADD FIGURE HERE**
>
> Suggested file: `figures/06_stable_ps_geodesic.png`
>
> Suggested content: stable-PS count map using geodesic distance.

<!-- Replace with:
![Stable PS filtering using geodesic distance](figures/06_stable_ps_geodesic.png)
-->

---

## Main stable-PS parameters

### `--stable-ps-radius`

Maximum spatial distance used when testing PS continuity.

```text
default = 2000
```

### `--history-steps`

Number of previous PS timesteps considered during continuity analysis.

```text
default = 5
```

### `--min-segment-time`

Minimum duration required for an accepted continuous PS segment.

```text
default = 120
```

### `--gap-factor`

Controls the tolerated time gap inside a continuous PS segment.

```text
gap tolerance = gap factor × estimated PS timestep
```

Default:

```text
3
```

### `--hard-cap-fraction`

Maximum allowed single-segment duration relative to nominal simulation
duration.

```text
default = 0.90
```

### `--allow-single-segment`

To allow long single-segment cases:

```bash
--allow-single-segment
```

---

# Protocol 4 — Rotor drift tracking

The `rotor-track` protocol follows PS detections over time and identifies
moving rotor tracks using global one-to-one assignment.

Tracks are subsequently processed through:

```text
temporal-gap splitting
        ↓
rolling-drift splitting
        ↓
minimum-duration filtering
        ↓
accepted rotor segments
```

Both Euclidean and geodesic distances are available.

---

## 4A. Rotor tracking — Euclidean

```bash
mkdir -p examples/results/PS_results/rotor_track_euclidean

psfilter rotor-track \
    --mesh-points examples/data/example_mesh.pts \
    --points-time examples/results/PS_results/example_mesh.pts_t \
    --tcl-summary examples/results/tcl_results/global_tcl_AP.txt \
    --output-dir examples/results/PS_results/rotor_track_euclidean \
    --distance-mode euclidean \
    --track-radius 2000 \
    --max-gap-steps 5 \
    --min-segment-time 120 \
    --gap-factor 3 \
    --drift-window-time 120 \
    --drift-radius-factor 2 \
    --hard-cap-fraction 0.90
```

Main outputs include:

```text
rotor_track_counts_cycle_<N>.dat
rotor_track_counts_ptst.dat
rotor_track_counts_normalized_ptst.dat
rotor_track_hit_counts_ptst.dat
rotor_track_duration_ptst.dat
rotor_tracks_counts_report.csv
rotor_summary.csv
rotor_track_segments.csv
rotor_track_drift_splits.csv
rotor_track_points.csv
rotor_track_mapping_report.csv
rotor_tracking_run_summary.csv
```

### Figure placeholder — rotor tracking, Euclidean

> **ADD FIGURE HERE**
>
> Suggested file: `figures/07_rotor_tracking_euclidean.png`
>
> Suggested content: rotor trajectories, rotor count map, or drift paths using
> Euclidean distance.

<!-- Replace with:
![Rotor tracking using Euclidean distance](figures/07_rotor_tracking_euclidean.png)
-->

---

## 4B. Rotor tracking — geodesic

```bash
mkdir -p examples/results/PS_results/rotor_track_geodesic

psfilter rotor-track \
    --mesh-points examples/data/example_mesh.pts \
    --mesh-elements examples/data/example_mesh.elem \
    --points-time examples/results/PS_results/example_mesh.pts_t \
    --tcl-summary examples/results/tcl_results/global_tcl_AP.txt \
    --output-dir examples/results/PS_results/rotor_track_geodesic \
    --distance-mode geodesic \
    --track-radius 2000 \
    --max-gap-steps 5 \
    --min-segment-time 120 \
    --gap-factor 3 \
    --drift-window-time 120 \
    --drift-radius-factor 2 \
    --hard-cap-fraction 0.90
```

### Figure placeholder — rotor tracking, geodesic

> **ADD FIGURE HERE**
>
> Suggested file: `figures/08_rotor_tracking_geodesic.png`
>
> Suggested content: tracked rotor trajectories or drift map using geodesic
> distance.

<!-- Replace with:
![Rotor tracking using geodesic distance](figures/08_rotor_tracking_geodesic.png)
-->

---

## Main rotor-tracking parameters

### `--track-radius`

Maximum distance allowed when matching a current PS detection to an active
rotor track.

```text
default = 2000
```

### `--max-gap-steps`

Controls track identity.

It determines how long a track remains available for matching through missing
PS detections.

```text
default = 5
```

### `--gap-factor`

Controls temporal continuity after the track has already been constructed.

```text
gap tolerance = gap factor × estimated PS timestep
```

A sufficiently large gap splits one track into separate temporal segments.

Therefore:

```text
max-gap-steps
    -> Can this detection still belong to the same track?

gap-factor
    -> Is this part of the track still one continuous segment?
```

### `--drift-window-time`

Controls the rolling time interval over which rotor displacement is evaluated.

```text
default = 120
```

The standalone command also accepts:

```text
--drift-min-time
```

as an alias.

### `--drift-radius-factor`

Maximum allowed rolling displacement:

```text
drift radius factor × track radius
```

With the defaults:

```text
2 × 2000 = 4000
```

If the displacement exceeds this value over the rolling drift window, the
segment is split.

### `--min-segment-time`

Applied after temporal-gap and drift splitting.

A final segment shorter than this duration is rejected.

```text
default = 120
```

Therefore:

```text
drift-window-time
    -> controls rolling-drift evaluation and splitting

min-segment-time
    -> controls final segment acceptance
```

---

# Protocol 5 — Complete workflow with `all`

After testing the stages separately, the complete workflow can be run with
`psfilter all`.

The pipeline combines:

```text
TCL calculation
      ↓
activation-sequence analysis
      ↓
PS detection
      ↓
PS coordinate mapping
      ↓
stable PS filtering
      ↓
rotor tracking
```

A basic example is:

```bash
psfilter all \
    --work-dir examples \
    --mesh data/example_mesh \
    --ascii-vm-input data/transmembrane_v.dat \
    --vm-igb data/vm.igb \
    --tcl-output-dir results/tcl_results \
    --ps-output-dir results/PS_results \
    --dt 1.0 \
    --reference-index 7
```

Paths passed to `all` are resolved relative to `--work-dir`.

The example above therefore resolves the main inputs to:

```text
examples/data/example_mesh.pts
examples/data/example_mesh.elem
examples/data/transmembrane_v.dat
examples/data/vm.igb
```

---

## Complete pipeline — Euclidean

For the faster complete workflow:

```bash
psfilter all \
    --work-dir examples \
    --mesh data/example_mesh \
    --ascii-vm-input data/transmembrane_v.dat \
    --vm-igb data/vm.igb \
    --tcl-output-dir results/tcl_results \
    --ps-output-dir results/PS_results \
    --dt 1.0 \
    --reference-index 7 \
    --ps-distance-mode euclidean \
    --track-distance-mode euclidean
```

No geodesic graph is built.

---

## Complete pipeline — geodesic

```bash
psfilter all \
    --work-dir examples \
    --mesh data/example_mesh \
    --ascii-vm-input data/transmembrane_v.dat \
    --vm-igb data/vm.igb \
    --tcl-output-dir results/tcl_results \
    --ps-output-dir results/PS_results \
    --dt 1.0 \
    --reference-index 7 \
    --ps-distance-mode geodesic \
    --track-distance-mode geodesic
```

A single geodesic graph is reused when both PS-analysis stages require it.

---

## Complete pipeline — mixed modes

For a large mesh, a useful compromise is:

```bash
psfilter all \
    --work-dir examples \
    --mesh data/example_mesh \
    --ascii-vm-input data/transmembrane_v.dat \
    --vm-igb data/vm.igb \
    --tcl-output-dir results/tcl_results \
    --ps-output-dir results/PS_results \
    --dt 1.0 \
    --reference-index 7 \
    --ps-distance-mode euclidean \
    --track-distance-mode geodesic
```

This gives:

```text
stable PS filtering -> Euclidean
rotor tracking      -> geodesic
```

Only the rotor-tracking stage requires the graph.

---

## Reuse an existing `.pts_t`

If PS detection has already been completed:

```bash
psfilter all \
    --work-dir examples \
    --mesh data/example_mesh \
    --ascii-vm-input data/transmembrane_v.dat \
    --vm-igb data/vm.igb \
    --points-time results/PS_results/example_mesh.pts_t \
    --tcl-output-dir results/tcl_results \
    --ps-output-dir results/PS_results \
    --dt 1.0 \
    --reference-index 7 \
    --skip-ps-detection \
    --ps-distance-mode euclidean \
    --track-distance-mode euclidean
```

When `--skip-ps-detection` is used, the `.pts_t` file must already exist and
be nonempty.

---

## Recalculate existing results

PSFilter reuses existing nonempty numerical outputs by default.

Use:

```bash
--overwrite
```

to intentionally recalculate them.

Example:

```bash
psfilter all \
    --work-dir examples \
    --mesh data/example_mesh \
    --ascii-vm-input data/transmembrane_v.dat \
    --vm-igb data/vm.igb \
    --tcl-output-dir results/tcl_results \
    --ps-output-dir results/PS_results \
    --dt 1.0 \
    --reference-index 7 \
    --ps-distance-mode euclidean \
    --track-distance-mode euclidean \
    --overwrite
```

Use `--overwrite` carefully for large geodesic calculations.

### Figure placeholder — complete pipeline

> **ADD FIGURE HERE**
>
> Suggested file: `figures/09_complete_pipeline.png`
>
> Suggested content: overview of the original simulation, signals, raw PSs,
> stable PS map, and rotor tracks.

<!-- Replace with:
![Complete PSFilter workflow](figures/09_complete_pipeline.png)
-->

---

# Figure placeholder — input mesh

> **ADD FIGURE HERE**
>
> Suggested file: `figures/01_input_mesh.png`
>
> Suggested content: the example mesh, anatomical regions, and/or signal
> locations.

<!-- Replace with:
![Example input mesh](figures/01_input_mesh.png)
-->

---

# Recommended order for a new user

Run the protocols separately first:

```text
1. psfilter tcl
2. psfilter ps-detect
3. psfilter ps-coords --distance-mode euclidean
4. psfilter rotor-track --distance-mode euclidean
5. psfilter ps-coords --distance-mode geodesic
6. psfilter rotor-track --distance-mode geodesic
7. psfilter all
```

This makes it easier to identify whether a problem comes from:

```text
TCL / signal processing
PS detection
mesh mapping
stable-PS filtering
rotor tracking
geodesic graph construction
```

For large meshes, start with Euclidean distance before running the geodesic
versions.

---

# Expected final example structure

After running the individual protocols:

```text
examples/
├── README.md
├── data/
│   ├── example_mesh.pts
│   ├── example_mesh.elem
│   ├── transmembrane_v.dat
│   ├── signal_labels.txt
│   └── vm.igb
│
├── results/
│   ├── tcl_results/
│   │   ├── global_tcl_AP.txt
│   │   ├── tcl_peaks_AP.csv
│   │   ├── tcl_per_signal_AP.csv
│   │   ├── tcl_summary.csv
│   │   ├── cycle_activations.csv
│   │   ├── cycle_summary.csv
│   │   ├── global_activation_order.csv
│   │   ├── activation_sequences.csv
│   │   └── sequence_counts.csv
│   │
│   └── PS_results/
│       ├── example_mesh.pts_t
│       ├── ps_coords_euclidean/
│       │   └── ...
│       ├── ps_coords_geodesic/
│       │   └── ...
│       ├── rotor_track_euclidean/
│       │   └── ...
│       └── rotor_track_geodesic/
│           └── ...
│
└── figures/
    ├── 01_input_mesh.png
    ├── 02_tcl_results.png
    ├── 03_signal_sequences.png
    ├── 04_ps_detection.png
    ├── 05_stable_ps_euclidean.png
    ├── 06_stable_ps_geodesic.png
    ├── 07_rotor_tracking_euclidean.png
    ├── 08_rotor_tracking_geodesic.png
    └── 09_complete_pipeline.png
```

---

# Reproducibility

For each example run, record:

```text
PSFilter version / Git commit
Python version
igbhead version
igbfilament version
mesh
reference TCL
PS timestep
distance mode
stable-PS radius
track radius
history steps
max gap steps
gap factor
minimum segment time
drift window time
drift radius factor
hard-cap fraction
```

These parameters can materially affect the detected stable PSs and rotor
tracks.

---

# Notes

- PSFilter is research software.
- The example parameter values are not universal physiological thresholds.
- Geodesic calculations can be substantially slower than Euclidean
  calculations on large meshes.
- The full `vm.igb` file should be downloaded separately rather than committed
  directly to the Git repository.
- The smaller derived files allow most of the package to be tested without
  downloading the complete simulation output.
