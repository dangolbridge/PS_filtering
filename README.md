# PSFilter

PSFilter is a Python package and command-line workflow for tachycardia cycle
length (TCL) analysis, activation-pattern analysis, phase-singularity (PS)
detection, coordinate-based PS filtering, and moving rotor tracking on cardiac
meshes.

The package is designed around openCARP-style files and external IGB utilities.
It supports individual analysis stages as well as a complete `all` pipeline.

## Main capabilities

- Calculate TCL from extracted transmembrane-voltage signals.
- Analyse cycle-by-cycle activation order and recurring activation sequences.
- Generate phase-singularity `.pts_t` files using `igbhead` and `igbfilament`.
- Map PS coordinates to mesh vertices.
- Filter stable PS locations using either:
  - Euclidean distance; or
  - bounded mesh-geodesic distance.
- Build moving rotor tracks using:
  - global one-to-one assignment;
  - Euclidean or mesh-geodesic matching;
  - temporal-gap splitting;
  - rolling-drift splitting;
  - segment-duration filtering;
  - occupancy, movement, and gap diagnostics.
- Run all stages through one CLI.
- Reuse mesh coordinates, mapped PS data, and one geodesic graph when possible.

## Requirements

### Python

- Python 3.10 or newer.

### Python dependencies

Installed through `pyproject.toml`:

- NumPy
- pandas
- SciPy

### External command-line programs

Some stages require programs that are not installed by `pip`:

- `meshtool`
- `igbextract`
- `igbhead`
- `igbfilament`

These executables must either be available in `PATH` or be supplied explicitly
through the corresponding CLI arguments.

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

The CLI exposes both Euclidean and geodesic implementations. The distance
backend is selected at run time.

## Installation

Clone the repository:

```bash
git clone https://github.com/dangolbridge/PS_filtering.git
cd PS_filtering
```

Create an environment with a modern Python version. For example, with Conda:

```bash
conda create -n psfilter-env python=3.11 pip
conda activate psfilter-env
```

Install the package in editable mode:

```bash
python -m pip install --upgrade pip
python -m pip install -e .
```

Check the installation:

```bash
psfilter --help
```

The equivalent module command is:

```bash
python -m psfilter --help
```

The root-level launcher can also be used:

```bash
python run_psfilter.py --help
```

## CLI commands

PSFilter exposes five subcommands:

```text
psfilter tcl
psfilter ps-detect
psfilter ps-coords
psfilter rotor-track
psfilter all
```

Display detailed help for a specific command with:

```bash
psfilter tcl --help
psfilter ps-detect --help
psfilter ps-coords --help
psfilter rotor-track --help
psfilter all --help
```

## Distance backends

Both `ps-coords` and `rotor-track` support:

```text
geodesic
euclidean
```

For standalone commands, select the backend with:

```bash
--distance-mode {geodesic,euclidean}
```

The default is `geodesic`.

For the complete pipeline, the two PS stages are configured independently:

```bash
--ps-distance-mode {geodesic,euclidean}
--track-distance-mode {geodesic,euclidean}
```

Both default to `geodesic`.

### Geodesic mode

Geodesic mode builds a graph from the `.elem` connectivity. Mesh-edge weights
are the Euclidean lengths of connected mesh vertices, and bounded graph
searches are used to estimate surface-aware distances.

Use the `.elem` file corresponding to the surface on which the PS coordinates
were generated. A volumetric mesh can allow paths through the atrial wall
rather than along the intended surface.

### Euclidean mode

Euclidean mode uses direct straight-line distance between mapped mesh
coordinates and does not construct a geodesic graph.

### Performance note

Geodesic `ps-coords` can be substantially slower than geodesic rotor tracking.
`ps-coords` checks PS continuity against several previous timesteps, which can
produce many local graph-distance searches.

A practical mixed workflow for large meshes is:

```bash
psfilter all \
    --work-dir /path/to/simulation \
    --mesh Reentry_surface_iac \
    --dt 1.0 \
    --reference-index 7 \
    --skip-ps-detection \
    --ps-distance-mode euclidean \
    --track-distance-mode geodesic
```

If surface-aware distance is not needed in either stage:

```bash
psfilter all \
    --work-dir /path/to/simulation \
    --mesh Reentry_surface_iac \
    --dt 1.0 \
    --reference-index 7 \
    --skip-ps-detection \
    --ps-distance-mode euclidean \
    --track-distance-mode euclidean
```

When both PS stages are Euclidean, no geodesic graph is built.

---

# 1. TCL and activation analysis

The `tcl` command calculates TCL and, unless disabled, activation-pattern
outputs.

## Basic usage

```bash
psfilter tcl \
    --ascii-vm-input /path/to/transmembrane_v.dat \
    --tcl-output-dir /path/to/tcl_results \
    --dt 1.0 \
    --reference-index 7
```

### Required arguments

| Argument | Meaning |
|---|---|
| `--ascii-vm-input` | Path to the ASCII transmembrane-voltage signal file. |
| `--dt` | Time interval between consecutive signal samples. |
| `--reference-index` | Signal index used as the reference for the mean TCL report. |

`--ascii-vm-input` is required by the standalone `tcl` parser.

### TCL and activation options

| Argument | Default | Meaning |
|---|---:|---|
| `--tcl-output-dir` | `./tcl_results` | Output directory for TCL and activation results. |
| `--derivative-threshold` | `5.0` | Derivative threshold used for peak detection. |
| `--min-peak-distance` | `10.0` | Minimum accepted separation between detected peaks. |
| `--labels-file` | none | Optional file containing signal labels. |
| `--derivative-per-time` | off | Calculate the derivative per unit time rather than per sample. |
| `--fraction-window` | `0.9` | TCL fraction used to define the activation-sequence window. |
| `--fraction-cluster` | `0.1` | TCL fraction used when grouping similar activation sequences. |
| `--max-groups` | `20` | Maximum number of activation-sequence groups used by the stopping logic. |
| `--no-stop-on-max-groups` | off | Do not stop when `max-groups` is reached. |
| `--skip-activation` | off | Calculate TCL only and skip activation-pattern analysis. |
| `--overwrite` | off | Recalculate numerical outputs instead of reusing existing nonempty results. |

## Main TCL outputs

```text
tcl_peaks_AP.csv
tcl_per_signal_AP.csv
global_tcl_AP.txt
tcl_summary.csv
```

When activation analysis is enabled:

```text
cycle_activations.csv
cycle_summary.csv
global_activation_order.csv
activation_sequences.csv
sequence_counts.csv
```

## Extracting the ASCII signal file when it is missing

If the requested `--ascii-vm-input` file does not already exist, PSFilter can
call `meshtool` and `igbextract` to generate it.

The current implementation checks for all of these extraction inputs:

```text
--vm-igb
--transformed-points
--mesh
--query-points
--node-indices
--meshtool
--igbextract
```

The extraction search radius is controlled by:

```text
--search-radius
```

with default `50.0`.

Example:

```bash
psfilter tcl \
    --ascii-vm-input /path/to/transmembrane_v.dat \
    --tcl-output-dir /path/to/tcl_results \
    --dt 1.0 \
    --reference-index 7 \
    --vm-igb /path/to/vm.igb \
    --transformed-points /path/to/query_locations.pts \
    --mesh /path/to/Reentry_surface_iac \
    --query-points /path/to/tcl_query_points.pts \
    --node-indices /path/to/node_indices.txt \
    --meshtool /path/to/meshtool \
    --igbextract /path/to/igbextract \
    --search-radius 50
```

An existing nonempty ASCII signal file is reused even when `--overwrite` is
used for the numerical outputs.

---

# 2. Phase-singularity detection

The `ps-detect` command creates a phase-singularity `.pts_t` file from
`vm.igb`.

## Basic usage

```bash
psfilter ps-detect \
    --vm-igb /path/to/vm.igb \
    --mesh /path/to/Reentry_surface_iac \
    --output /path/to/PS_results/Reentry_surface_iac.pts_t
```

### Required arguments

```text
--vm-igb
--mesh
--output
```

### Detection options

| Argument | Default | Meaning |
|---|---:|---|
| `--igbhead` | `igbhead` | `igbhead` executable or path. |
| `--igbfilament` | `igbfilament` | `igbfilament` executable or path. |
| `--cleaned-igb` | none | Optional path for the cleaned temporary IGB file. |
| `--threshold` | `-50.0` | Voltage threshold used by the PS-detection workflow. |
| `--filament-dt` | `8.0` | Timestep used for filament/PS extraction. |
| `--overwrite` | off | Recalculate an existing PS output. |
| `--keep-cleaned-igb` | off | Keep the cleaned IGB file after PS detection. |
| `--dry-run` | off | Show the expected operation without performing the calculation. |

Example with explicit executable paths:

```bash
psfilter ps-detect \
    --vm-igb /path/to/vm.igb \
    --mesh /path/to/Reentry_surface_iac \
    --output /path/to/PS_results/Reentry_surface_iac.pts_t \
    --igbhead /path/to/igbhead \
    --igbfilament /path/to/igbfilament \
    --threshold -50 \
    --filament-dt 8
```

---

# 3. Coordinate-based stable PS filtering

The `ps-coords` command identifies spatially and temporally persistent PS
locations. It supports both Euclidean and geodesic distance.

## Geodesic example

```bash
psfilter ps-coords \
    --mesh-points /path/to/Reentry_surface_iac.pts \
    --mesh-elements /path/to/Reentry_surface_iac.elem \
    --points-time /path/to/PS_results/Reentry_surface_iac.pts_t \
    --tcl-summary /path/to/tcl_results/global_tcl_AP.txt \
    --output-dir /path/to/PS_results/ps_coords_geodesic \
    --distance-mode geodesic \
    --stable-ps-radius 2000 \
    --history-steps 5 \
    --min-segment-time 120 \
    --gap-factor 3 \
    --hard-cap-fraction 0.90
```

## Euclidean example

```bash
psfilter ps-coords \
    --mesh-points /path/to/Reentry_surface_iac.pts \
    --points-time /path/to/PS_results/Reentry_surface_iac.pts_t \
    --tcl-summary /path/to/tcl_results/global_tcl_AP.txt \
    --output-dir /path/to/PS_results/ps_coords_euclidean \
    --distance-mode euclidean \
    --stable-ps-radius 2000 \
    --history-steps 5 \
    --min-segment-time 120 \
    --gap-factor 3 \
    --hard-cap-fraction 0.90
```

### Required arguments

```text
--mesh-points
--points-time
--tcl-summary
--output-dir
```

`--mesh-elements` is optional at the parser level. In geodesic mode, if it is
not supplied, PSFilter looks for an `.elem` file with the same basename as
`--mesh-points`.

## `ps-coords` parameters

| Argument | Default | Meaning |
|---|---:|---|
| `--stable-ps-radius` | `2000.0` | Maximum spatial distance used when testing PS continuity. |
| `--history-steps` | `5` | Number of recent PS timesteps considered for continuity. |
| `--min-segment-time` | `120.0` | Minimum duration required for an accepted continuous PS segment. |
| `--gap-factor` | `3.0` | Temporal-gap tolerance expressed as a multiple of the estimated PS timestep. |
| `--hard-cap-fraction` | `0.90` | Maximum allowed segment duration relative to the nominal simulation duration. |
| `--allow-single-segment` | off | Allow a vertex represented by a single continuous segment. |
| `--max-mapping-distance` | none | Optional maximum distance allowed when mapping PS coordinates to mesh vertices. |
| `--distance-mode` | `geodesic` | Select `geodesic` or `euclidean` distance. |
| `--geodesic-cache-size` | `200000` | Maximum number of finite geodesic pair distances retained in the cache. |
| `--edge-chunk-size` | `500000` | Temporary edge-buffer size used during geodesic graph construction. |
| `--overwrite` | off | Recalculate existing nonempty outputs. |

`--geodesic-cache-size` and `--edge-chunk-size` affect only geodesic mode.

## Main `ps-coords` outputs

```text
rotor_counts_cycle_<N>.dat
rotor_counts_ptst.dat
rotor_counts_normalized_ptst.dat
rotor_segments_filtered.csv
rotor_lifetimes.csv
ps_mapping_report.csv
ps_coords_summary.csv
```

The summary records the selected distance method.

### Performance suggestions

For large meshes:

1. use `--distance-mode euclidean` when surface-aware distance is not required;
2. reduce `--history-steps` when scientifically appropriate;
3. avoid unnecessarily large `--stable-ps-radius` values;
4. use the correct surface `.elem` file for geodesic calculations;
5. use `--skip-ps-coords` in the complete pipeline when this stage is not required.

---

# 4. Moving rotor tracking

The `rotor-track` command builds moving PS tracks using global one-to-one
assignment and filters them using temporal-gap, rolling-drift, and
segment-duration criteria.

## Geodesic example

```bash
psfilter rotor-track \
    --mesh-points /path/to/Reentry_surface_iac.pts \
    --mesh-elements /path/to/Reentry_surface_iac.elem \
    --points-time /path/to/PS_results/Reentry_surface_iac.pts_t \
    --tcl-summary /path/to/tcl_results/global_tcl_AP.txt \
    --output-dir /path/to/PS_results/rotor_track_geodesic \
    --distance-mode geodesic \
    --track-radius 2000 \
    --max-gap-steps 5 \
    --min-segment-time 120 \
    --gap-factor 3 \
    --drift-window-time 120 \
    --drift-radius-factor 2 \
    --hard-cap-fraction 0.90
```

## Euclidean example

```bash
psfilter rotor-track \
    --mesh-points /path/to/Reentry_surface_iac.pts \
    --points-time /path/to/PS_results/Reentry_surface_iac.pts_t \
    --tcl-summary /path/to/tcl_results/global_tcl_AP.txt \
    --output-dir /path/to/PS_results/rotor_track_euclidean \
    --distance-mode euclidean \
    --track-radius 2000 \
    --max-gap-steps 5 \
    --min-segment-time 120 \
    --gap-factor 3 \
    --drift-window-time 120 \
    --drift-radius-factor 2 \
    --hard-cap-fraction 0.90
```

### Required arguments

```text
--mesh-points
--points-time
--tcl-summary
--output-dir
```

As with `ps-coords`, `--mesh-elements` is needed by the analysis only for
geodesic distance. If omitted, the matching `.elem` path is inferred from
`--mesh-points`.

## `rotor-track` parameters

| Argument | Default | Meaning |
|---|---:|---|
| `--track-radius` | `2000.0` | Maximum distance allowed when matching a current PS detection to an active rotor track. |
| `--max-gap-steps` | `5` | Maximum timestep-index separation for retaining a track as matchable. |
| `--min-segment-time` | `120.0` | Minimum duration required for a final rotor segment to be valid. |
| `--gap-factor` | `3.0` | Temporal-gap threshold used to split an already-built track into continuous segments. |
| `--drift-window-time` | `120.0` | Rolling time window over which spatial drift is evaluated. |
| `--drift-min-time` | alias | Alias for `--drift-window-time` in the standalone command. |
| `--drift-radius-factor` | `2.0` | Maximum rolling displacement expressed as a multiple of `track-radius`. |
| `--hard-cap-fraction` | `0.90` | Maximum segment duration relative to nominal simulation duration. |
| `--count-every-hit` | off | Count every valid rotor hit instead of unique vertex membership per track in the compatibility count map. |
| `--max-mapping-distance` | none | Optional maximum PS-to-mesh mapping distance. |
| `--distance-mode` | `geodesic` | Select `geodesic` or `euclidean` matching. |
| `--geodesic-cache-size` | `200000` | Maximum number of finite geodesic pair distances retained in the cache. |
| `--edge-chunk-size` | `500000` | Temporary edge-buffer size used during geodesic graph construction. |
| `--overwrite` | off | Recalculate existing nonempty outputs. |

A negative `--hard-cap-fraction` disables the maximum segment-duration cap for
the standalone rotor command.

## `max-gap-steps` versus `gap-factor`

These parameters both concern missing detections, but they are used at
different stages.

`--max-gap-steps` is used during **track construction**. It determines whether
a new PS detection can still be matched to a previously active track. It is
expressed in timestep-index units.

`--gap-factor` is used later during **track segmentation**. The temporal gap
threshold is:

```text
gap_tolerance_time = gap_factor × estimated_PS_timestep
```

If two consecutive detections inside an already-built track are separated by
more than this time, the track is split into separate continuous segments.

In short:

```text
max-gap-steps  -> track identity tolerance
gap-factor     -> continuous-segment tolerance
```

## `drift-window-time` versus `min-segment-time`

`--drift-window-time` controls how far backward in time the tracker looks when
evaluating rotor movement.

If the rolling displacement exceeds:

```text
drift-radius-factor × track-radius
```

the segment is split.

`--min-segment-time` is applied after temporal-gap and rolling-drift splitting.
A final segment shorter than this value is rejected as too short.

In short:

```text
drift-window-time -> controls rolling-drift evaluation and splitting
min-segment-time  -> controls final segment acceptance
```

## Main rotor-tracking outputs

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

---

# 5. Complete pipeline

The `all` command runs the complete analysis workflow.

A basic command is:

```bash
psfilter all \
    --work-dir /path/to/simulation \
    --mesh Reentry_surface_iac \
    --dt 1.0 \
    --reference-index 7
```

The individual commands are useful for checking and tuning each stage before
using the complete pipeline.

## Default paths

For the example above, the path resolver uses:

```text
/path/to/simulation/transmembrane_v.dat
/path/to/simulation/vm.igb
/path/to/simulation/Reentry_surface_iac.pts
/path/to/simulation/Reentry_surface_iac.elem
/path/to/simulation/tcl_results/
/path/to/simulation/PS_results/
/path/to/simulation/PS_results/Reentry_surface_iac.pts_t
```

Relative paths supplied to `all` are resolved relative to `--work-dir`.

## Main `all` path arguments

| Argument | Default |
|---|---|
| `--work-dir` | `.` |
| `--mesh` | required |
| `--ascii-vm-input` | `<work-dir>/transmembrane_v.dat` when omitted |
| `--vm-igb` | `<work-dir>/vm.igb` when omitted |
| `--mesh-points` | `<mesh>.pts` when omitted |
| `--mesh-elements` | `<mesh>.elem` when omitted |
| `--points-time` | `<ps-output-dir>/<mesh-name>.pts_t` when omitted |
| `--tcl-output-dir` | `<work-dir>/tcl_results` |
| `--ps-output-dir` | `<work-dir>/PS_results` |
| `--ps-coords-output-dir` | mode-specific subdirectory under `PS_results` |
| `--rotor-output-dir` | mode-specific subdirectory under `PS_results` |

Default PS-analysis output directories are:

```text
PS_results/ps_coords_geodesic/
PS_results/rotor_track_geodesic/
```

When Euclidean modes are selected:

```text
PS_results/ps_coords_euclidean/
PS_results/rotor_track_euclidean/
```

## TCL options in `all`

```text
--dt
--reference-index
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

`--dt` and `--reference-index` are required.

Extraction-related options are:

```text
--transformed-points
--query-points
--node-indices
--meshtool
--igbextract
--search-radius
```

## PS-detection options in `all`

```text
--igbhead
--igbfilament
--cleaned-igb
--threshold
--filament-dt
--keep-cleaned-igb
--skip-ps-detection
```

When `--skip-ps-detection` is used, an existing nonempty `.pts_t` file must
already be present at the resolved `--points-time` path.

## Stable-PS options in `all`

| Argument | Default |
|---|---:|
| `--stable-ps-radius` | `2000.0` |
| `--history-steps` | `5` |
| `--ps-min-segment-time` | `120.0` |
| `--ps-gap-factor` | `3.0` |
| `--ps-hard-cap-fraction` | `0.90` |
| `--allow-single-segment` | off |
| `--skip-ps-coords` | off |
| `--ps-distance-mode` | `geodesic` |
| `--ps-coords-output-dir` | mode-specific default |

## Rotor-tracking options in `all`

| Argument | Default |
|---|---:|
| `--track-radius` | `2000.0` |
| `--max-gap-steps` | `5` |
| `--track-min-segment-time` | `120.0` |
| `--track-gap-factor` | `3.0` |
| `--drift-window-time` | `120.0` |
| `--drift-radius-factor` | `2.0` |
| `--track-hard-cap-fraction` | `0.90` |
| `--count-every-hit` | off |
| `--skip-rotor-track` | off |
| `--track-distance-mode` | `geodesic` |
| `--rotor-output-dir` | mode-specific default |

Shared PS-analysis options are:

```text
--max-mapping-distance
--geodesic-cache-size
--edge-chunk-size
```

The whole pipeline can be forced to recalculate numerical outputs with:

```bash
--overwrite
```

## Output reuse

If TCL outputs already exist and are nonempty, they are reused unless
`--overwrite` is specified.

If stable-PS or rotor-tracking principal outputs already exist and are
nonempty, those stages are reused unless `--overwrite` is specified.

An existing nonempty extracted `transmembrane_v.dat` is reused even when
`--overwrite` is supplied.

## Geodesic graph reuse in `all`

The pipeline builds a geodesic graph only when at least one PS stage actually
needs geodesic distance.

For example:

```text
ps-coords   = euclidean
rotor-track = geodesic
```

builds one graph for rotor tracking.

If both stages are Euclidean, no graph is created.

## Pipeline summary

The `all` command writes:

```text
psfilter_pipeline_summary.csv
```

Stage names include the selected backend, for example:

```text
tcl
ps-detect
ps-coords-euclidean
rotor-track-geodesic
```

---

# Core parameter summary

| CLI parameter | Default | Main role |
|---|---:|---|
| `--stable-ps-radius` | `2000` | Spatial PS-continuity radius in `ps-coords`. |
| `--history-steps` | `5` | Previous timesteps used for stable-PS continuity. |
| `--track-radius` | `2000` | Maximum PS-to-track matching distance. |
| `--max-gap-steps` | `5` | Track-identity tolerance in timestep-index units. |
| `--gap-factor` | `3` | Continuous-segment gap tolerance in standalone commands. |
| `--min-segment-time` | `120` | Minimum valid-segment duration in standalone commands. |
| `--drift-window-time` | `120` | Rolling time window used to evaluate rotor drift. |
| `--drift-radius-factor` | `2` | Maximum rolling drift as a multiple of `track-radius`. |
| `--hard-cap-fraction` | `0.90` | Maximum single-segment duration relative to nominal simulation duration. |
| `--max-mapping-distance` | none | Optional PS-to-mesh mapping-distance limit. |
| `--geodesic-cache-size` | `200000` | Maximum geodesic distance-cache size. |
| `--edge-chunk-size` | `500000` | Edge chunk size during graph construction. |

In `all`, stage-specific versions are used for several parameters:

```text
--ps-min-segment-time
--ps-gap-factor
--ps-hard-cap-fraction

--track-min-segment-time
--track-gap-factor
--track-hard-cap-fraction
```

---

# Troubleshooting

## `psfilter: command not found`

Activate the environment in which PSFilter was installed:

```bash
conda activate psfilter-env
```

Then verify:

```bash
which python
python --version
python -m pip show psfilter
which psfilter
```

For an editable installation:

```bash
cd /path/to/PS_filtering
python -m pip install -e .
hash -r
```

The module form can also be tested:

```bash
python -m psfilter --help
```

## Help output does not show the latest arguments

Check the imported CLI module:

```bash
python -c "import psfilter.cli; print(psfilter.cli.__file__)"
```

To run directly from the repository source:

```bash
cd /path/to/PS_filtering
PYTHONPATH="$PWD/src" python -m psfilter rotor-track --help
```

The standalone PS commands should show:

```text
--distance-mode {geodesic,euclidean}
```

The complete pipeline should show:

```text
--ps-distance-mode {geodesic,euclidean}
--track-distance-mode {geodesic,euclidean}
```

## Geodesic `ps-coords` is very slow

Check:

- mesh size;
- PS density;
- `--history-steps`;
- `--stable-ps-radius`;
- whether the correct surface `.elem` file is being used.

A faster alternative is:

```bash
psfilter ps-coords ... --distance-mode euclidean
```

or:

```bash
psfilter all ... \
    --ps-distance-mode euclidean \
    --track-distance-mode geodesic
```

## Geodesic results cross an unexpected anatomical connection

The graph follows the topology encoded in the `.elem` file. Every explicit
edge or connection in that graph is available to the geodesic search.

Use the surface mesh whose topology matches the intended analysis.

## Existing outputs are unexpectedly reused

Use:

```bash
--overwrite
```

to force recalculation.

## `--skip-ps-detection` fails

When this option is used, the resolved `.pts_t` file must already exist and be
nonempty.

Check `--points-time`, or the default:

```text
<ps-output-dir>/<mesh-name>.pts_t
```

---

# Development

Install development dependencies:

```bash
python -m pip install -e ".[dev]"
```

Run a syntax check:

```bash
python -m compileall -q src run_psfilter.py
```

Run tests:

```bash
pytest
```

After changing CLI argument names, check every command help page:

```bash
python -m psfilter tcl --help
python -m psfilter ps-detect --help
python -m psfilter ps-coords --help
python -m psfilter rotor-track --help
python -m psfilter all --help
```

---

# Research-use notice

PSFilter is research software. Analysis parameters, PS sampling interval, mesh
topology, external-tool versions, and signal preprocessing choices can
materially affect the results.

For reproducibility, record:

- the complete CLI command;
- PSFilter version or Git commit;
- mesh files;
- external-tool versions;
- key distance, gap, drift, and duration parameters.

Before a public software release, add an appropriate license and project
citation.
