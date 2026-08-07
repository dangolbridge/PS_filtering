# PSFilter

PSFilter is a Python package and command-line workflow for tachycardia cycle
length analysis, activation-order analysis, phase-singularity detection,
coordinate-based PS filtering, and moving rotor tracking on cardiac meshes.

The package is designed around openCARP-style files and external IGB utilities.
It supports both individual analysis stages and a complete pipeline.

## Main capabilities

- Calculate tachycardia cycle length (TCL) from extracted transmembrane-voltage signals.
- Analyse cycle-by-cycle activation order and recurring activation sequences.
- Run `igbhead` and `igbfilament` to generate phase-singularity `.pts_t` files.
- Map phase-singularity coordinates to mesh vertices.
- Filter coordinate-based PS detections using either:
  - Euclidean distance; or
  - bounded mesh-geodesic distance.
- Build moving rotor tracks with:
  - global one-to-one assignment;
  - Euclidean or mesh-geodesic matching;
  - rolling drift splitting;
  - segment-level filtering;
  - occupancy and time-gap diagnostics.
- Run all stages through one CLI while reusing loaded mesh coordinates, parsed
  `.pts_t` data, and PS-to-vertex mapping.
- Build and reuse a `MeshGeodesicGraph` only when at least one requested stage
  uses geodesic distance.

## Distance backends and performance

Both `ps-coords` and `rotor-track` support two distance backends:

```text
geodesic
euclidean
```

The default is `geodesic`.

For the standalone commands, select the backend with:

```bash
--distance-mode geodesic
```

or:

```bash
--distance-mode euclidean
```

For the complete `all` pipeline, the two stages can be configured independently:

```bash
--ps-distance-mode {geodesic,euclidean}
--track-distance-mode {geodesic,euclidean}
```

The geodesic `PS_coords` analysis can be substantially slower than geodesic
rotor tracking. For every current PS vertex, `PS_coords` checks continuity
against PS vertices from several recent timesteps. Depending on PS density,
this can produce many bounded graph-distance queries, and the analysis is
performed for individual TCL cycles as well as for the complete time range.

Rotor tracking normally performs fewer geodesic queries because it compares
only active tracks with PS points in the current timestep.

For large meshes or dense `.pts_t` files, a practical mixed-backend workflow is:

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

This keeps the faster Euclidean backend for coordinate-history filtering while
retaining surface-aware geodesic matching for rotor tracking.

When surface topology is not required for either analysis, both stages can use
Euclidean distance:

```bash
psfilter all ... \
    --ps-distance-mode euclidean \
    --track-distance-mode euclidean
```

In that case, no geodesic graph is built.

## Requirements

### Python

- Python 3.10 or newer

### Python dependencies

Installed automatically through `pyproject.toml`:

- NumPy
- pandas
- SciPy

### External command-line programs

Some stages require tools that are not installed through `pip`:

- `meshtool`
- `igbextract`
- `igbhead`
- `igbfilament`

These programs must be installed separately and either available in `PATH` or
provided explicitly through CLI arguments.

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

The CLI exposes both Euclidean and geodesic implementations. The backend is
selected at run time rather than by changing modules manually.

## Installation

Clone the repository:

```bash
git clone https://github.com/dangolbridge/PS_filtering.git
cd PS_filtering
```

Create and activate a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
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

A root-level compatibility launcher is also available:

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

Inspect a command before running it:

```bash
psfilter rotor-track --help
```

## Quick start

### Run a complete analysis with existing input files

This example assumes that `transmembrane_v.dat` and the PS `.pts_t` file
already exist:

```bash
psfilter all \
    --work-dir /path/to/simulation \
    --mesh Reentry_surface_iac \
    --dt 1.0 \
    --reference-index 7 \
    --skip-ps-detection
```

Because the default distance backend is geodesic for both PS stages, this also
uses:

```text
/path/to/simulation/Reentry_surface_iac.pts
/path/to/simulation/Reentry_surface_iac.elem
```

To use Euclidean distance for both PS stages:

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

To use Euclidean `PS_coords` and geodesic rotor tracking:

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

### Default paths used by `all`

For:

```bash
psfilter all \
    --work-dir /path/to/simulation \
    --mesh Reentry_surface_iac \
    --dt 1.0 \
    --reference-index 7
```

the default input paths are:

```text
/path/to/simulation/vm.igb
/path/to/simulation/transmembrane_v.dat
/path/to/simulation/Reentry_surface_iac.pts
/path/to/simulation/Reentry_surface_iac.elem
/path/to/simulation/tcl_results/
/path/to/simulation/PS_results/Reentry_surface_iac.pts_t
```

The default PS-analysis output directories are mode-specific:

```text
/path/to/simulation/PS_results/ps_coords_geodesic/
/path/to/simulation/PS_results/rotor_track_geodesic/
```

If Euclidean modes are selected, the corresponding defaults become:

```text
/path/to/simulation/PS_results/ps_coords_euclidean/
/path/to/simulation/PS_results/rotor_track_euclidean/
```

Custom directories can be selected with:

```text
--ps-coords-output-dir
--rotor-output-dir
```

Relative paths passed to `all` are resolved relative to `--work-dir`.

## TCL and activation analysis

Run TCL and activation-pattern analysis from an existing signal file:

```bash
psfilter tcl \
    --input /path/to/transmembrane_v.dat \
    --output-dir /path/to/tcl_results \
    --dt 1.0 \
    --reference-index 7
```

Useful options:

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

### Extract `transmembrane_v.dat` when it is missing

```bash
psfilter tcl \
    --input /path/to/transmembrane_v.dat \
    --output-dir /path/to/tcl_results \
    --dt 1.0 \
    --reference-index 7 \
    --vm-igb /path/to/vm.igb \
    --transformed-points /path/to/query_locations.pts \
    --mesh /path/to/Reentry_surface_iac \
    --query-points /path/to/tcl_query_points.pts \
    --node-indices /path/to/node_indices.txt \
    --meshtool /path/to/meshtool \
    --igbextract /path/to/igbextract
```

## Phase-singularity detection

Create a `.pts_t` file from `vm.igb`:

```bash
psfilter ps-detect \
    --vm-igb /path/to/vm.igb \
    --mesh /path/to/Reentry_surface_iac \
    --output /path/to/PS_results/Reentry_surface_iac.pts_t \
    --threshold -50 \
    --filament-dt 8
```

Specify executable paths when they are not available in `PATH`:

```bash
psfilter ps-detect \
    --vm-igb /path/to/vm.igb \
    --mesh /path/to/Reentry_surface_iac \
    --output /path/to/PS_results/Reentry_surface_iac.pts_t \
    --igbhead /path/to/igbhead \
    --igbfilament /path/to/igbfilament
```

Additional options:

```text
--cleaned-igb
--overwrite
--keep-cleaned-igb
--dry-run
```

## Coordinate-based PS filtering

`ps-coords` supports both Euclidean and mesh-geodesic continuity.

### Geodesic mode

```bash
psfilter ps-coords \
    --mesh-points /path/to/Reentry_surface_iac.pts \
    --mesh-elements /path/to/Reentry_surface_iac.elem \
    --points-time /path/to/PS_results/Reentry_surface_iac.pts_t \
    --tcl-summary /path/to/tcl_results/global_tcl_AP.txt \
    --output-dir /path/to/PS_results/ps_coords_geodesic \
    --distance-mode geodesic \
    --radius 2000 \
    --history-steps 5 \
    --min-segment-time 120 \
    --gap-factor 3 \
    --hard-cap-fraction 0.90
```

In geodesic mode, the current PS vertex is counted when it is within the
configured mesh-geodesic radius of a PS vertex from one of the previous
`history_steps` timesteps.

`--mesh-elements` is required for the graph calculation. If it is omitted, the
CLI looks for an `.elem` file with the same basename as `--mesh-points`.

### Euclidean mode

```bash
psfilter ps-coords \
    --mesh-points /path/to/Reentry_surface_iac.pts \
    --points-time /path/to/PS_results/Reentry_surface_iac.pts_t \
    --tcl-summary /path/to/tcl_results/global_tcl_AP.txt \
    --output-dir /path/to/PS_results/ps_coords_euclidean \
    --distance-mode euclidean \
    --radius 2000 \
    --history-steps 5 \
    --min-segment-time 120 \
    --gap-factor 3 \
    --hard-cap-fraction 0.90
```

Euclidean mode uses straight-line coordinate distance and does not build a
geodesic graph.

Main outputs:

```text
rotor_counts_cycle_<N>.dat
rotor_counts_ptst.dat
rotor_counts_normalized_ptst.dat
rotor_segments_filtered.csv
rotor_lifetimes.csv
ps_mapping_report.csv
ps_coords_summary.csv
```

The summary file records the selected backend in `distance_method`.

By default, the historical single-segment rejection behavior is retained.
Allow a vertex represented by one continuous segment with:

```bash
--allow-single-segment
```

### Runtime controls

The most influential parameters are:

```text
--distance-mode
--radius
--history-steps
--geodesic-cache-size
--edge-chunk-size
```

`--geodesic-cache-size` and `--edge-chunk-size` affect only geodesic mode.

To reduce runtime:

1. use `--distance-mode euclidean` when surface-aware distance is not required;
2. reduce `--history-steps`;
3. verify that `--radius` is not unnecessarily large;
4. use the correct surface mesh rather than a volumetric mesh in geodesic mode;
5. run rotor tracking instead when coordinate-history maps are not required;
6. skip this stage in the complete pipeline with `--skip-ps-coords`.

## Rotor tracking

`rotor-track` supports both Euclidean and mesh-geodesic matching while retaining
the same global assignment, rolling-drift, gap, and segment-filtering logic.

### Geodesic mode

```bash
psfilter rotor-track \
    --mesh-points /path/to/Reentry_surface_iac.pts \
    --mesh-elements /path/to/Reentry_surface_iac.elem \
    --points-time /path/to/PS_results/Reentry_surface_iac.pts_t \
    --tcl-summary /path/to/tcl_results/global_tcl_AP.txt \
    --output-dir /path/to/PS_results/rotor_track_geodesic \
    --distance-mode geodesic \
    --radius 2000 \
    --max-gap-steps 5 \
    --min-segment-time 120 \
    --gap-factor 3 \
    --drift-window-time 120 \
    --drift-radius-factor 2 \
    --hard-cap-fraction 0.90
```

### Euclidean mode

```bash
psfilter rotor-track \
    --mesh-points /path/to/Reentry_surface_iac.pts \
    --points-time /path/to/PS_results/Reentry_surface_iac.pts_t \
    --tcl-summary /path/to/tcl_results/global_tcl_AP.txt \
    --output-dir /path/to/PS_results/rotor_track_euclidean \
    --distance-mode euclidean \
    --radius 2000 \
    --max-gap-steps 5 \
    --min-segment-time 120 \
    --gap-factor 3 \
    --drift-window-time 120 \
    --drift-radius-factor 2 \
    --hard-cap-fraction 0.90
```

Euclidean mode does not require an `.elem` file and does not construct a
geodesic graph.

The tracker:

1. maps PS coordinates to mesh vertices;
2. globally assigns current PS points to active tracks;
3. evaluates candidate matching distance with the selected backend;
4. splits tracks at temporal gaps;
5. splits segments when rolling drift exceeds the configured limit;
6. filters and counts valid segments only.

Main outputs:

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

The run summary records the selected backend in `distance_method`.

Disable the maximum segment-duration filter:

```bash
--hard-cap-fraction -1
```

Count every valid track visit in the compatibility count map:

```bash
--count-every-hit
```

Separate hit-count and duration maps are written regardless of this option.

## Complete pipeline

```bash
psfilter all \
    --work-dir /path/to/simulation \
    --mesh Reentry_surface_iac \
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

Distance-backend controls:

```text
--ps-distance-mode {geodesic,euclidean}
--track-distance-mode {geodesic,euclidean}
```

Both default to `geodesic`, and they can be selected independently.

Mode-specific output-directory controls:

```text
--ps-coords-output-dir
--rotor-output-dir
```

The `all` command reuses:

- mesh coordinates;
- parsed `.pts_t` data;
- PS-to-vertex mapping.

A single `MeshGeodesicGraph` is built and reused only when at least one PS
analysis stage actually requires geodesic distance. If both stages use
Euclidean distance, no graph is built.

A stage summary is written to:

```text
psfilter_pipeline_summary.csv
```

The stage names include the selected backend, for example:

```text
ps-coords-euclidean
rotor-track-geodesic
```

## Geodesic-distance implementation

This section applies when `geodesic` is selected for `ps-coords` or
`rotor-track`.

The mesh graph is built from the edges in the openCARP `.elem` file. Edge
weights are the Euclidean lengths of connected mesh vertices.

Distances are evaluated with bounded, target-aware Dijkstra searches. The code
does not build an all-pairs distance matrix.

The resulting distance is an edge-based approximation of surface geodesic
distance. Use the `.elem` file corresponding to the surface on which the PS
coordinates were generated.

Do not unintentionally use a volumetric tetrahedral mesh when the desired
distance is along the atrial surface. A volumetric graph can permit paths
through the wall rather than along the surface.

The graph follows the topology encoded in the `.elem` file. Explicit
inter-layer or interatrial connections are therefore available to the
geodesic search.

For `euclidean` mode, distances are calculated directly from mesh coordinates
and no graph is constructed.

## Core analysis parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `distance_mode` | `geodesic` | Distance backend for standalone `ps-coords` or `rotor-track` |
| `ps_distance_mode` | `geodesic` | Distance backend for `ps-coords` inside `all` |
| `track_distance_mode` | `geodesic` | Distance backend for rotor tracking inside `all` |
| `radius` | 2000 | Continuation or track-matching radius in the selected distance backend |
| `history_steps` | 5 | Previous timesteps used by `PS_coords` |
| `max_gap_steps` | 5 | Timestep-index gap allowed for track continuation |
| `min_segment_time` | 120 | Minimum accepted continuous-segment duration |
| `gap_factor` | 3 | Time-gap threshold as a multiple of PS timestep |
| `drift_window_time` | 120 | Rolling time window for drift evaluation |
| `drift_radius_factor` | 2 | Maximum rolling drift as a multiple of radius |
| `hard_cap_fraction` | 0.90 | Maximum segment duration relative to nominal simulation duration |
| `geodesic_cache_size` | 200000 | Maximum cached finite pair distances in geodesic mode |
| `edge_chunk_size` | 500000 | Temporary edge-buffer size during geodesic graph construction |

The nominal simulation duration used by the segment-duration cap is:

```text
number of TCL cycles × reference mean TCL
```

## Output reuse and overwrite

Existing nonempty principal output files are reused by default.

Use:

```bash
--overwrite
```

to recalculate numerical outputs.

For safety and efficiency, an existing nonempty extracted
`transmembrane_v.dat` is reused even when numerical outputs are overwritten.

## Troubleshooting

### `psfilter: command not found`

Activate the virtual environment and reinstall:

```bash
source .venv/bin/activate
python -m pip install -e .
```

The module form should still work:

```bash
python -m psfilter --help
```

### A required external executable is missing

Provide its full path:

```bash
--igbhead /full/path/to/igbhead
--igbfilament /full/path/to/igbfilament
--meshtool /full/path/to/meshtool
--igbextract /full/path/to/igbextract
```

### Geodesic `PS_coords` appears to hang

It may simply be performing a large number of bounded graph searches. Check:

- mesh size and topology;
- number of PS points per timestep;
- `--history-steps`;
- geodesic radius;
- whether the correct surface `.elem` file is being used.

For a faster coordinate-history analysis, use:

```bash
psfilter ps-coords ... --distance-mode euclidean
```

For the complete pipeline, a practical mixed-backend configuration is:

```bash
psfilter all ... \
    --ps-distance-mode euclidean \
    --track-distance-mode geodesic
```

You can also skip coordinate-based filtering entirely:

```bash
psfilter all ... --skip-ps-coords
```

### Results are joined across anatomically separate surfaces

Confirm that the `.elem` file contains the intended topology. Geodesic
tracking cannot cross disconnected mesh components, but it can cross every
explicit edge or connection present in the graph.

### Euclidean mode still asks for or builds an `.elem` graph

Confirm that the intended mode is actually selected:

```bash
psfilter ps-coords --help
psfilter rotor-track --help
```

The standalone commands should show:

```text
--distance-mode {geodesic,euclidean}
```

For the complete pipeline, use:

```text
--ps-distance-mode {geodesic,euclidean}
--track-distance-mode {geodesic,euclidean}
```

Also confirm that the active Python environment is loading the expected package:

```bash
python -c "import psfilter.cli; print(psfilter.cli.__file__)"
```

### Existing outputs are skipped unexpectedly

Use:

```bash
--overwrite
```

Also confirm that the package version of `required_outputs_exist()` treats
zero-byte files as incomplete.

## Development installation

Install development tools:

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

## Research-use notice

This package is research software. Analysis parameters, mesh topology,
external-tool versions, and signal preprocessing choices can materially affect
the results. Record the complete CLI command and software versions used for
each analysis.

Before public release, add an appropriate license file and project citation.
