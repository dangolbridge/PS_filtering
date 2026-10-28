# PSFilter

**PSFilter** is a Python package and command-line workflow for post-processing cardiac electrophysiology simulations. It provides tools for tachycardia cycle length (TCL) analysis, activation-order analysis, phase-singularity (PS) detection, PS-to-mesh mapping, local PS-continuity filtering, and moving rotor tracking.

> **Current release:** v1.1.0  
> **Status:** research software

## What changed in v1.1.0

Version 1.1 introduces a revised PS-filtering and rotor-tracking workflow with:

- Euclidean and mesh-geodesic distance backends;
- global one-to-one assignment for rotor tracking;
- explicit short-gap handling;
- rolling-drift splitting;
- segment-level duration filtering;
- per-cycle, aggregated-cycle, and full-range maps;
- expanded mapping, track, segment, occupancy, and drift diagnostics;
- a cleaned package and CLI structure.

## Main concepts

### `ps-coords`: local continuous PS activity

`ps-coords` identifies regions where PS activity is spatially and temporally continuous. At each timestep, a current PS is counted when it remains within a selected spatial radius of PS activity observed in recent timesteps.

This analysis answers:

> **Where is there stable or recurrent local PS activity?**

It does not assign a persistent identity to an individual rotor.

### `rotor-track`: individual rotor trajectories

`rotor-track` explicitly connects PS detections through time into individual tracks using global one-to-one assignment, a spatial matching radius, short-gap handling, temporal-gap splitting, rolling-drift splitting, and minimum-duration filtering.

This analysis answers:

> **Can an individual rotor core be followed continuously through time, and how does it move?**

For trajectory, lifetime, drift, and individual-rotor analyses, `rotor-track` is the primary method.

## Requirements

### Python

Python **3.10 or newer** is recommended.

### Python dependencies

Installed through `pyproject.toml`:

- NumPy
- pandas
- SciPy

### External command-line tools

Some stages require external tools that are not installed through `pip`:

- `meshtool`
- `igbextract`
- `igbhead`
- `igbfilament`

They must either be available in `PATH` or supplied explicitly to the relevant command.

## Installation

```bash
git clone https://github.com/dangolbridge/PS_filtering.git
cd PS_filtering

python -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -e .
```

Check the installation:

```bash
psfilter --version
psfilter --help
```

Equivalent module form:

```bash
python -m psfilter --help
```

Compatibility launcher:

```bash
python run_psfilter.py --help
```

## CLI commands

```text
psfilter tcl
psfilter ps-detect
psfilter ps-coords
psfilter rotor-track
psfilter all
```

Inspect any command with:

```bash
psfilter <command> --help
```

## Repository layout

```text
PS_filtering/
├── pyproject.toml
├── README.md
├── run_psfilter.py
└── src/
    └── psfilter/
        ├── __init__.py
        ├── __main__.py
        ├── activation.py
        ├── cli.py
        ├── commands.py
        ├── external_tools.py
        ├── io.py
        ├── mapping.py
        ├── mesh_geodesic.py
        ├── phase_singularity.py
        ├── preprocessing.py
        ├── ps_coordinates.py
        ├── ps_coordinates_geodesic.py
        ├── rotor_tracking.py
        ├── rotor_tracking_geodesic.py
        └── tcl.py
```

The Euclidean and geodesic implementations are both retained. The CLI selects the requested backend through distance-mode arguments.

## Typical workflow

```text
transmembrane voltage
        │
        ▼
TCL / activation analysis
        │
        ▼
phase-singularity detection
        │
        ▼
PS-to-mesh mapping
        │
        ├──────────────► ps-coords
        │                 local stable/recurrent PS activity
        │
        └──────────────► rotor-track
                          individual moving rotor trajectories
```

The stages can be run independently or through `psfilter all`. For first use and debugging, running the stages independently is recommended.

## 1. TCL and activation analysis

```bash
psfilter tcl \
    --input /path/to/transmembrane_v.dat \
    --output-dir /path/to/tcl_results \
    --dt 1.0 \
    --reference-index 7
```

Main outputs:

```text
tcl_peaks_AP.csv
tcl_per_signal_AP.csv
global_tcl_AP.txt
tcl_summary.csv
cycle_activations.csv
cycle_summary.csv
global_activation_order.csv
activation_sequences.csv
sequence_counts.csv
```

Useful options include:

```text
--derivative-threshold
--min-peak-distance
--labels-file
--derivative-per-time
--fraction-window
--fraction-cluster
--max-groups
--no-stop-on-max-groups
--skip-activation
```

If `transmembrane_v.dat` does not already exist, the TCL stage can use `meshtool` and `igbextract` to extract signals from an IGB file. Inspect the full options with:

```bash
psfilter tcl --help
```

## 2. Phase-singularity detection

```bash
psfilter ps-detect \
    --vm-igb /path/to/vm.igb \
    --mesh /path/to/mesh_basename \
    --output /path/to/PS_results/phase_singularities.pts_t \
    --threshold -50 \
    --filament-dt 8
```

If the required programs are not in `PATH`:

```bash
psfilter ps-detect \
    --vm-igb /path/to/vm.igb \
    --mesh /path/to/mesh_basename \
    --output /path/to/PS_results/phase_singularities.pts_t \
    --igbhead /full/path/to/igbhead \
    --igbfilament /full/path/to/igbfilament
```

`--filament-dt` controls the temporal interval used by the filament-detection procedure. It should not automatically be interpreted as the sampling interval of the resulting `.pts_t` file. PSFilter estimates the `.pts_t` timestep from timestamps stored in that file.

## 3. Coordinate-based PS filtering

```bash
psfilter ps-coords \
    --mesh-points /path/to/mesh.pts \
    --mesh-elements /path/to/mesh.elem \
    --points-time /path/to/PS_results/phase_singularities.pts_t \
    --tcl-summary /path/to/tcl_results/global_tcl_AP.txt \
    --output-dir /path/to/PS_results/ps_coords \
    --distance-mode geodesic
```

Use Euclidean distance with:

```bash
--distance-mode euclidean
```

### Current v1.1 PS-coordinate defaults

| Parameter | Default | Meaning |
|---|---:|---|
| spatial radius | 2000 | Maximum spatial continuation distance |
| `history_steps` | 10 | Number of recent PS timesteps searched for local continuity |
| `min_segment_time` | 120 ms | Minimum accepted continuous duration |
| `gap_factor` | 10 | Maximum temporal separation as a multiple of estimated PS timestep |
| `hard_cap_fraction` | 0.90 | Maximum segment duration relative to nominal simulation duration |
| single-segment handling | reject by default | Use `--allow-single-segment` to retain one sufficiently long continuous segment |

The spatial unit is inherited from the mesh. For meshes stored in micrometres, `2000` corresponds to approximately `2 mm`.

### Main PS-coordinate outputs

```text
rotor_counts_cycle_<N>.dat
rotor_counts_cycles_aggregated.dat
rotor_counts_ptst.dat
rotor_counts_normalized_ptst.dat
rotor_segments_filtered.csv
rotor_lifetimes.csv
ps_mapping_report.csv
ps_coords_summary.csv
```

Interpretation:

```text
rotor_counts_cycle_<N>.dat
    independent continuity analysis within one TCL cycle

rotor_counts_cycles_aggregated.dat
    raw sum of independently calculated cycle maps

rotor_counts_ptst.dat
    full-range continuity analysis across the complete simulation

rotor_counts_normalized_ptst.dat
    normalized full-range filtered map
```

Because cycle analyses restart temporal history at each TCL boundary, the aggregated-cycle map does not necessarily equal the full-range map.

## 4. Rotor tracking

```bash
psfilter rotor-track \
    --mesh-points /path/to/mesh.pts \
    --mesh-elements /path/to/mesh.elem \
    --points-time /path/to/PS_results/phase_singularities.pts_t \
    --tcl-summary /path/to/tcl_results/global_tcl_AP.txt \
    --output-dir /path/to/PS_results/rotor_tracking \
    --distance-mode geodesic
```

Euclidean tracking:

```bash
--distance-mode euclidean
```

### Current v1.1 rotor-tracking defaults

| Parameter | Default | Meaning |
|---|---:|---|
| track matching radius | 2000 | Maximum distance for linking a current PS to an active track |
| `max_gap_steps` | 5 | Maximum timestep-index gap allowed during track construction |
| `min_segment_time` | 120 ms | Minimum accepted final segment duration |
| `gap_factor` | 5 | Temporal-gap threshold as a multiple of estimated PS timestep |
| `drift_window_time` | 120 ms | Rolling window used to assess rotor-core displacement |
| `drift_radius_factor` | 2 | Maximum rolling drift as a multiple of matching radius |
| `hard_cap_fraction` | 1.0 | Maximum allowed segment duration relative to nominal simulation duration |

With a radius of `2000`, the default rolling-drift threshold is:

```text
2 × 2000 = 4000 mesh units
```

For meshes in micrometres:

```text
4000 µm = 4 mm
```

over a default rolling window of:

```text
120 ms
```

The default hard-cap fraction of `1.0` allows a valid rotor segment to persist for the complete nominal simulation duration.

### Tracking procedure

For every timestep, the tracker:

1. keeps tracks that are still eligible for continuation;
2. computes candidate distances between active tracks and current PS detections;
3. performs global one-to-one assignment;
4. starts new tracks from unmatched PS detections;
5. later separates tracks at excessive temporal gaps;
6. separates tracks when rolling displacement exceeds the drift threshold;
7. applies minimum-duration and optional maximum-duration criteria to the resulting segments.

Global assignment means that one current PS cannot be assigned to multiple tracks and one track cannot receive multiple PS points at the same timestep.

### Geodesic tracking

With `--distance-mode geodesic`, track matching and rolling-drift checks use shortest paths along the mesh-edge graph.

This is an edge-based approximation of surface geodesic distance rather than an exact continuous geodesic calculation.

The graph follows the topology encoded in the `.elem` file. Therefore:

- disconnected surfaces remain disconnected;
- explicit mesh connections can be used by the shortest path;
- a volumetric mesh can permit paths through the tissue volume.

Use the `.elem` file corresponding to the topology on which the desired distance should be measured.

### Main rotor-tracking outputs

```text
rotor_track_counts_cycle_<N>.dat
rotor_track_counts_cycles_aggregated.dat
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

Interpretation:

```text
rotor_track_counts_cycle_<N>.dat
    rotor-track counts obtained independently within one TCL cycle

rotor_track_counts_cycles_aggregated.dat
    sum of independently calculated cycle maps

rotor_track_counts_ptst.dat
    accepted rotor-track visitation map calculated across the full time range

rotor_track_hit_counts_ptst.dat
    number of valid PS detections at each vertex

rotor_track_duration_ptst.dat
    estimated valid occupancy duration at each vertex

rotor_track_counts_normalized_ptst.dat
    normalized full-range track-count map
```

By default, an accepted track contributes once per visited vertex. Use:

```bash
--count-every-hit
```

when every valid visit should contribute to the compatibility count map.

## 5. Complete pipeline

```bash
psfilter all \
    --work-dir /path/to/simulation \
    --mesh mesh_basename \
    --dt 1.0 \
    --reference-index 7
```

Important stage controls:

```text
--skip-activation
--skip-ps-detection
--skip-ps-coords
--skip-rotor-track
--overwrite
```

Distance modes can be selected independently:

```bash
psfilter all \
    --work-dir /path/to/simulation \
    --mesh mesh_basename \
    --dt 1.0 \
    --reference-index 7 \
    --ps-distance-mode euclidean \
    --track-distance-mode geodesic
```

The pipeline writes:

```text
psfilter_pipeline_summary.csv
```

## Mapping PS coordinates to the mesh

PS coordinates are mapped to mesh vertices before PS analysis.

Mapping diagnostics include:

- number of input PS detections;
- number of mapped detections;
- number of dropped detections;
- maximum mapping distance.

An optional maximum mapping distance can reject PS coordinates that are too far from the mesh.

## Euclidean vs geodesic distance

### Euclidean

Uses straight-line distance between mesh coordinates.

Advantages:

- faster;
- simple;
- useful for sensitivity analysis and debugging.

Limitation:

- points that are close in 3-D can be anatomically distant along a surface.

### Geodesic

Uses shortest-path distance over mesh edges.

Advantages:

- respects mesh connectivity;
- avoids direct straight-line jumps between nearby but disconnected surfaces.

Limitations:

- computationally more expensive;
- depends on supplied `.elem` topology and mesh resolution;
- is an edge-graph approximation rather than an exact continuous geodesic.

## Performance notes

Geodesic `ps-coords` can be substantially slower than geodesic `rotor-track`, especially for large meshes or dense `.pts_t` files.

For large-mesh analyses, it can be useful to run:

```bash
psfilter rotor-track ...
```

directly, or:

```bash
psfilter all ... --skip-ps-coords
```

Useful geodesic-performance controls include:

```text
--geodesic-cache-size
--edge-chunk-size
```

## Output reuse and overwrite

Existing nonempty principal outputs are reused by default.

To recalculate:

```bash
--overwrite
```

For sensitivity studies, separate output directories are strongly recommended:

```text
PS_results/
├── rotor_baseline/
├── rotor_duration_90ms/
├── rotor_duration_120ms/
├── rotor_drift_3mm/
├── rotor_drift_4mm/
└── rotor_drift_5mm/
```

## Parameter sensitivity

The defaults are working defaults for the current analysis workflow, not universal physiological constants.

For a new dataset, sensitivity analysis is recommended for:

```text
minimum segment duration
track matching radius
maximum timestep gap
rolling-drift threshold
```

A robust default should ideally lie in a region where modest parameter changes do not qualitatively change the scientific conclusions.

## Reproducibility recommendations

For each analysis, record:

- PSFilter version or Git commit;
- complete CLI command;
- mesh `.pts` and `.elem` files;
- coordinate units;
- `.pts_t` file;
- reference TCL;
- distance backend;
- spatial radius;
- temporal-gap parameters;
- minimum segment duration;
- drift-window and drift-distance parameters;
- external-tool versions used for PS detection or signal extraction.

For multi-model studies, retain `rotor_tracking_run_summary.csv` and the corresponding command for every run.

## Development checks

```bash
python -m compileall -q src run_psfilter.py
psfilter --help
psfilter ps-coords --help
psfilter rotor-track --help
```

If development dependencies are configured:

```bash
python -m pip install -e ".[dev]"
pytest
```

## Troubleshooting

### `psfilter: command not found`

```bash
python -m pip install -e .
which python
which psfilter
python -c "import psfilter; print(psfilter.__file__)"
```

### Geodesic analysis is unexpectedly slow

Check mesh size, topology, PS density per timestep, selected radius, `history_steps`, and whether the intended surface `.elem` file is being used.

For comparison:

```bash
--distance-mode euclidean
```

### Results cross anatomically separate surfaces

Inspect the `.elem` connectivity. Geodesic distance follows the graph encoded by the mesh.

### Existing results are unexpectedly reused

Use:

```bash
--overwrite
```

or a new output directory.

## Research-use notice

PSFilter is research software for computational cardiac electrophysiology.

Results can depend on phase-singularity detection settings, mesh topology and resolution, coordinate units, signal preprocessing, temporal and spatial thresholds, distance backend, and external-tool versions.

Default parameters should therefore be documented and, where appropriate, evaluated through sensitivity analysis before being used for scientific conclusions.

## Version history

### v1.1.0

Revised PS filtering and rotor tracking with Euclidean/geodesic backends, global assignment, rolling-drift splitting, revised temporal continuity handling, aggregated cycle maps, and expanded diagnostics.

### v1.0.0

Previous main version, preserved as the `v1.0.0` Git tag and the `release/1.0` branch.

## Citation

A formal software citation can be added here when available.

For scientific publications using PSFilter, report the software version or Git commit together with the relevant analysis parameters.
