# Atlas Explorer MVP

Minimal static web scaffold for inspecting Atlas-0 baseline morphology and Atlas-1 transition structure.

## Folder structure

```text
atlas_explorer/
  index.html
  app.js
  three-grid.js
  style.css
  README.md
  vendor/
    three.module.js
```

## What it does

- Loads Atlas-0 baseline CSV + JSON with no backend.
- Renders a sortable and filterable run table.
- Draws a deterministic regime map from `anisotropy_score` and `coherence_score`.
- Adds an interactive 3D parameter space with axes `central_k`, `bandwidth`, and `amplitude`.
- Supports a wave / particle proxy color mode derived deterministically from recorded spread and compactness metrics.
- Reconstructs a 3D lattice replay with node intensity and edge-flow proxies from stored `metrics.centers`, `widths`, `energies`, and `coherences`.
- Renders the lattice replay through a local vendored `three.js` WebGL scene, with canvas fallback if WebGL is unavailable.
- Provides play/pause playback for the lattice replay timeline.
- Shows run-level metrics plus artifact images for the selected baseline run.
- Loads Atlas-1 transition CSV + JSON when present and renders a transition matrix plus persistence bar chart.

## Default data paths

The default config lives at the top of [`app.js`](./app.js):

```js
const defaultAtlasConfig = {
  csvPath: "../data/20260314_143113_stage10_atlas0_baseline_morphology.csv",
  jsonPath: "../data/20260314_143113_stage10_atlas0_baseline_morphology.json",
  transitionCsvPath: "../data/20260314_153736_stage10_atlas1_perturbation_resilience.csv",
  transitionJsonPath: "../data/20260314_153736_stage10_atlas1_perturbation_resilience.json",
  plotBasePath: "../plots/",
  latticeSide: 16,
  trajectoryRenderer: "three",
};
```

If you generate a newer Atlas run, update those paths to the new stamped files.

## Run locally

Serve the `HAOS-IIP` directory with a simple static server, then open `atlas_explorer/`.

```bash
cd /Volumes/Samsung\ T5/2026/HAOS/HAOS\ DOCS/HAOS-IIP
python3 -m http.server 8000
```

Open:

- `http://localhost:8000/atlas_explorer/`

## Notes

- Open the app through HTTP, not `file://`, because the browser needs `fetch()` access to CSV and JSON files.
- The detail panel resolves plot filenames deterministically from the configured Stage 10 timestamps and `run_id`.
- The grid-flow view is a deterministic proxy built from the recorded packet center, width, energy, and coherence traces, not a raw per-edge current dump.
- `three.module.js` is vendored locally under `atlas_explorer/vendor/` so the app stays a static local site and does not depend on a live CDN.
- Set `trajectoryRenderer: "canvas"` in `app.js` if you want to force the previous lightweight fallback.
- The wave / particle readout is a proxy view built from existing morphology metrics, not a literal physical classifier.
- Atlas-2 and later layers can plug in by extending the normalization functions in `app.js` rather than rewriting the UI shell.
