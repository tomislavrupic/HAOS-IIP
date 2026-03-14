const defaultAtlasConfig = {
  csvPath: "../data/20260314_143113_stage10_atlas0_baseline_morphology.csv",
  jsonPath: "../data/20260314_143113_stage10_atlas0_baseline_morphology.json",
  transitionCsvPath: "../data/20260314_153736_stage10_atlas1_perturbation_resilience.csv",
  transitionJsonPath: "../data/20260314_153736_stage10_atlas1_perturbation_resilience.json",
  plotBasePath: "../plots/",
  latticeSide: 16,
  trajectoryRenderer: "three",
};

const atlasConfig = {
  ...defaultAtlasConfig,
  ...(window.atlasConfig || {}),
};

const TABLE_COLUMNS = [
  { key: "run_id", label: "run_id", type: "string" },
  { key: "atlas_phase", label: "atlas_phase", type: "string" },
  { key: "regime_label", label: "regime_label", type: "string" },
  { key: "sector_label", label: "sector_label", type: "string" },
  { key: "proto_spacetime_score", label: "proto_spacetime_score", type: "number" },
  { key: "coherence_score", label: "coherence_score", type: "number" },
  { key: "anisotropy_score", label: "anisotropy_score", type: "number" },
  { key: "boundary_type", label: "boundary_type", type: "string" },
  { key: "operator_sector", label: "operator_sector", type: "string" },
];

const PLOT_DEFINITIONS = [
  { key: "field_snapshot", label: "Field snapshot" },
  { key: "spectrum_snapshot", label: "Spectrum snapshot" },
  { key: "width_trace", label: "Width trace" },
  { key: "constraint_trace", label: "Constraint trace" },
];

const PREFERRED_REGIME_ORDER = [
  "Ballistic coherent",
  "Ballistic dispersive",
  "Diffusive",
  "Localized",
  "Oscillatory trapped",
  "Chaotic or irregular",
  "Fragmenting",
  "Metastable structured",
  "Constraint-failing",
];

const REGIME_COLORS = {
  "Ballistic coherent": "#7ea2a8",
  "Ballistic dispersive": "#b49463",
  Diffusive: "#708da4",
  Localized: "#8b7a68",
  "Oscillatory trapped": "#9b8d64",
  "Chaotic or irregular": "#b36e60",
  Fragmenting: "#899267",
  "Metastable structured": "#6f9682",
  "Constraint-failing": "#8e6b6b",
};

const SECTOR_SHAPES = {
  "Scalar-dominant": "circle",
  "Transverse-dominant": "diamond",
  "Constraint-failing": "square",
  "Mixed-sector": "triangle",
};

const AXIS_INDEX = {
  x: 0,
  y: 1,
  z: 2,
};

const state = {
  atlas0: {
    rows: [],
    rowById: new Map(),
    executedByRunId: new Map(),
    metricExtents: null,
    stats: null,
    timestamp: null,
  },
  atlas1: {
    rows: [],
    executedByRunId: new Map(),
    transitionsByBaseline: new Map(),
    stats: null,
    timestamp: null,
    loaded: false,
  },
  filters: {
    search: "",
    regime: "all",
    sector: "all",
    boundary: "all",
  },
  sort: {
    key: "coherence_score",
    direction: "desc",
  },
  filteredRows: [],
  selectedRunId: null,
  detailToken: 0,
  plotResolutionCache: new Map(),
  assetProbeCache: new Map(),
  latticeCache: new Map(),
  view3d: {
    parameter: {
      rotationX: -0.52,
      rotationY: 0.82,
      zoom: 1.16,
      colorMode: "regime",
      hoverRunId: null,
      projected: [],
      drag: null,
    },
    trajectory: {
      rotationX: -0.34,
      rotationY: 0.72,
      zoom: 1.2,
      frameIndex: 0,
      projected: [],
      visibleNodeCount: 0,
      activeEdgeCount: 0,
      drag: null,
      lastRunId: null,
      playIntervalId: null,
    },
  },
};

const elements = {};

document.addEventListener("DOMContentLoaded", () => {
  cacheElements();
  bindEvents();
  void init();
});

async function init() {
  renderDataSources();
  try {
    await loadAtlas0();
    await loadAtlas1();
    populateFilterOptions();
    applyFilters();
    renderStaticPanels();
    renderTransitionPanel();
    renderSummary();
    renderSelection();
    renderParameterSpace();
    renderTrajectoryPlayback();
    await renderRunDetail();
  } catch (error) {
    handleFatalError(error);
  }
}

function cacheElements() {
  elements.atlas0Status = document.getElementById("atlas0-status");
  elements.atlas1Status = document.getElementById("atlas1-status");
  elements.dataSourceList = document.getElementById("data-source-list");
  elements.summaryChips = document.getElementById("summary-chips");
  elements.selectionBrief = document.getElementById("selection-brief");
  elements.searchInput = document.getElementById("search-input");
  elements.regimeFilter = document.getElementById("regime-filter");
  elements.sectorFilter = document.getElementById("sector-filter");
  elements.boundaryFilter = document.getElementById("boundary-filter");
  elements.clearFilters = document.getElementById("clear-filters");
  elements.resultCount = document.getElementById("result-count");
  elements.sortIndicator = document.getElementById("sort-indicator");
  elements.runTableBody = document.getElementById("run-table-body");
  elements.regimeMap = document.getElementById("regime-map");
  elements.mapLegend = document.getElementById("map-legend");
  elements.parameterSpaceCanvas = document.getElementById("parameter-space-canvas");
  elements.parameterSpaceInfo = document.getElementById("parameter-space-info");
  elements.parameterColorMode = document.getElementById("parameter-color-mode");
  elements.reset3dView = document.getElementById("reset-3d-view");
  elements.detailEmpty = document.getElementById("detail-empty");
  elements.detailContent = document.getElementById("detail-content");
  elements.detailMetrics = document.getElementById("detail-metrics");
  elements.detailMetadata = document.getElementById("detail-metadata");
  elements.detailTransitions = document.getElementById("detail-transitions");
  elements.detailImages = document.getElementById("detail-images");
  elements.trajectoryCanvas = document.getElementById("trajectory-canvas");
  elements.trajectoryPlayToggle = document.getElementById("trajectory-play-toggle");
  elements.trajectoryScrub = document.getElementById("trajectory-scrub");
  elements.trajectoryFrameInfo = document.getElementById("trajectory-frame-info");
  elements.resetTrajectoryView = document.getElementById("reset-trajectory-view");
  elements.transitionEmpty = document.getElementById("transition-empty");
  elements.transitionContent = document.getElementById("transition-content");
  elements.transitionMatrix = document.getElementById("transition-matrix");
  elements.persistenceBar = document.getElementById("persistence-bar");
  elements.transitionSummary = document.getElementById("transition-summary");
}

function bindEvents() {
  elements.searchInput.addEventListener("input", (event) => {
    state.filters.search = event.target.value.trim().toLowerCase();
    applyFilters();
  });

  elements.regimeFilter.addEventListener("change", (event) => {
    state.filters.regime = event.target.value;
    applyFilters();
  });

  elements.sectorFilter.addEventListener("change", (event) => {
    state.filters.sector = event.target.value;
    applyFilters();
  });

  elements.boundaryFilter.addEventListener("change", (event) => {
    state.filters.boundary = event.target.value;
    applyFilters();
  });

  elements.clearFilters.addEventListener("click", () => {
    state.filters = {
      search: "",
      regime: "all",
      sector: "all",
      boundary: "all",
    };
    elements.searchInput.value = "";
    elements.regimeFilter.value = "all";
    elements.sectorFilter.value = "all";
    elements.boundaryFilter.value = "all";
    applyFilters();
  });

  elements.runTableBody.addEventListener("click", (event) => {
    const row = event.target.closest("tr[data-run-id]");
    if (!row) {
      return;
    }
    selectRun(row.dataset.runId);
  });

  elements.regimeMap.addEventListener("click", (event) => {
    const point = event.target.closest("[data-run-id]");
    if (!point) {
      return;
    }
    selectRun(point.dataset.runId);
  });

  elements.parameterColorMode.addEventListener("change", (event) => {
    state.view3d.parameter.colorMode = event.target.value;
    renderParameterSpace();
  });

  elements.reset3dView.addEventListener("click", () => {
    state.view3d.parameter.rotationX = -0.52;
    state.view3d.parameter.rotationY = 0.82;
    renderParameterSpace();
  });

  elements.trajectoryPlayToggle.addEventListener("click", () => {
    toggleTrajectoryPlayback();
  });

  elements.resetTrajectoryView.addEventListener("click", () => {
    state.view3d.trajectory.rotationX = -0.34;
    state.view3d.trajectory.rotationY = 0.72;
    renderTrajectoryPlayback();
  });

  elements.trajectoryScrub.addEventListener("input", (event) => {
    stopTrajectoryPlayback();
    state.view3d.trajectory.frameIndex = toNumber(event.target.value);
    renderTrajectoryPlayback();
  });

  document.querySelectorAll(".sort-button").forEach((button) => {
    button.addEventListener("click", () => {
      const { sortKey } = button.dataset;
      if (state.sort.key === sortKey) {
        state.sort.direction = state.sort.direction === "asc" ? "desc" : "asc";
      } else {
        state.sort.key = sortKey;
        state.sort.direction = defaultDirectionForKey(sortKey);
      }
      applyFilters();
    });
  });

  bind3dCanvasInteractions();
  window.addEventListener("resize", () => {
    renderParameterSpace();
    renderTrajectoryPlayback();
  });
  window.addEventListener("atlas-three-ready", () => {
    renderTrajectoryPlayback();
  });
}

async function loadAtlas0() {
  elements.atlas0Status.textContent = "Loading";
  const [csvText, jsonPayload] = await Promise.all([
    fetchText(atlasConfig.csvPath),
    fetchJson(atlasConfig.jsonPath),
  ]);

  const rows = parseCsv(csvText);
  state.atlas0.executedByRunId = new Map(
    (jsonPayload.executed_runs || []).map((item) => [item.run.run_id, item]),
  );
  state.atlas0.timestamp = extractTimestamp(atlasConfig.csvPath) || extractTimestamp(atlasConfig.jsonPath);
  state.atlas0.rows = rows.map((row) => createRunRecord(row, state.atlas0.executedByRunId.get(row.run_id) || null));
  state.atlas0.rowById = new Map(state.atlas0.rows.map((row) => [row.run_id, row]));
  state.atlas0.metricExtents = createMetricExtents(state.atlas0.rows);
  state.atlas0.stats = createRegimeStats(state.atlas0.rows);
  state.selectedRunId = null;
  elements.atlas0Status.textContent = `${state.atlas0.rows.length} runs`;
}

async function loadAtlas1() {
  elements.atlas1Status.textContent = "Loading";
  try {
    const [csvText, jsonPayload] = await Promise.all([
      fetchText(atlasConfig.transitionCsvPath),
      fetchJson(atlasConfig.transitionJsonPath).catch(() => null),
    ]);
    const rows = parseCsv(csvText);
    const executedByRunId = new Map(
      ((jsonPayload && jsonPayload.executed_runs) || []).map((item) => [item.run.run_id, item]),
    );
    state.atlas1.executedByRunId = executedByRunId;
    state.atlas1.timestamp =
      extractTimestamp(atlasConfig.transitionCsvPath) || extractTimestamp(atlasConfig.transitionJsonPath);
    state.atlas1.rows = rows.map((row) =>
      createTransitionRecord(row, executedByRunId.get(row.run_id) || null),
    );
    state.atlas1.transitionsByBaseline = groupBy(state.atlas1.rows, "baseline_run_id");
    state.atlas1.stats = createTransitionStats(state.atlas1.rows);
    state.atlas1.loaded = state.atlas1.rows.length > 0;
    elements.atlas1Status.textContent = state.atlas1.loaded
      ? `${state.atlas1.rows.length} transitions`
      : "Not loaded";
  } catch (error) {
    state.atlas1.loaded = false;
    state.atlas1.rows = [];
    state.atlas1.transitionsByBaseline = new Map();
    state.atlas1.stats = null;
    elements.atlas1Status.textContent = "Optional";
    console.warn("Atlas-1 transition artifacts not loaded.", error);
  }
}

function createRunRecord(row, executedRun) {
  const metrics = executedRun?.metrics || {};
  const labels = executedRun?.labels || {};
  const integration = executedRun?.integration || {};
  return {
    run_id: row.run_id,
    atlas_phase: row.atlas_phase,
    graph_type: row.graph_type,
    kernel_type: row.kernel_type,
    operator_sector: row.operator_sector,
    boundary_type: row.boundary_type,
    initial_seed_type: row.initial_seed_type,
    central_k: toNumber(row.central_k),
    bandwidth: toNumber(row.bandwidth),
    amplitude: toNumber(row.amplitude),
    phase_pattern: row.phase_pattern,
    packet_count: toNumber(row.packet_count),
    random_seed: toNumber(row.random_seed),
    constraint_max: toNumber(row.constraint_max),
    sector_leakage: toNumber(row.sector_leakage),
    norm_drift: toNumber(row.norm_drift),
    anisotropy_score: toNumber(row.anisotropy_score),
    coherence_score: toNumber(row.coherence_score),
    proto_spacetime_score: toNumber(row.proto_spacetime_score),
    regime_label: row.regime_label,
    sector_label: row.sector_label,
    notes: row.notes || "",
    center_shift: toNumber(row.center_shift),
    width_ratio: toNumber(row.width_ratio),
    spectral_centroid_final: toNumber(row.spectral_centroid_final),
    spectral_spread_final: toNumber(row.spectral_spread_final),
    metrics,
    labels,
    integration,
    plotCandidates: buildAtlas0PlotCandidates(row.run_id),
  };
}

function createTransitionRecord(row, executedRun) {
  return {
    run_id: row.run_id,
    atlas_phase: row.atlas_phase,
    baseline_run_id: row.baseline_run_id,
    perturbation_axis: row.perturbation_axis,
    perturbation_strength: toNumber(row.perturbation_strength),
    operator_sector: row.operator_sector,
    boundary_type: row.boundary_type,
    random_seed: toNumber(row.random_seed),
    baseline_regime: row.baseline_regime,
    perturbed_regime: row.perturbed_regime,
    perturbed_regime_perturbed_constraint: row.perturbed_regime_perturbed_constraint,
    baseline_sector_label: row.baseline_sector_label,
    sector_label: row.sector_label,
    sector_label_perturbed_constraint: row.sector_label_perturbed_constraint,
    persistence_score: toNumber(row.persistence_score),
    persistence_score_perturbed_constraint: toNumber(row.persistence_score_perturbed_constraint),
    constraint_max: toNumber(row.constraint_max),
    constraint_max_perturbed: toNumber(row.constraint_max_perturbed),
    baseline_constraint_max: toNumber(row.baseline_constraint_max),
    sector_leakage: toNumber(row.sector_leakage),
    baseline_sector_leakage: toNumber(row.baseline_sector_leakage),
    norm_drift: toNumber(row.norm_drift),
    anisotropy_score: toNumber(row.anisotropy_score),
    coherence_score: toNumber(row.coherence_score),
    baseline_coherence_score: toNumber(row.baseline_coherence_score),
    coherence_delta: toNumber(row.coherence_delta),
    center_shift: toNumber(row.center_shift),
    width_ratio: toNumber(row.width_ratio),
    notes: row.notes || "",
    labels: executedRun?.labels || null,
    labels_perturbed_constraint: executedRun?.labels_perturbed_constraint || null,
    perturbation_metadata: executedRun?.perturbation_metadata || null,
  };
}

function createRegimeStats(rows) {
  const regimeCounts = countBy(rows, "regime_label");
  const sectorCounts = countBy(rows, "sector_label");
  const boundaryCounts = countBy(rows, "boundary_type");
  const meanCoherence =
    rows.reduce((sum, row) => sum + row.coherence_score, 0) / Math.max(rows.length, 1);

  return {
    totalRuns: rows.length,
    regimeCounts,
    sectorCounts,
    boundaryCounts,
    meanCoherence,
    topRegime: firstEntry(regimeCounts),
    topSector: firstEntry(sectorCounts),
  };
}

function createMetricExtents(rows) {
  return {
    coherence_score: dataExtent(rows.map((row) => row.coherence_score)),
    proto_spacetime_score: dataExtent(rows.map((row) => row.proto_spacetime_score)),
    width_ratio: dataExtent(rows.map((row) => row.width_ratio)),
    spectral_spread_final: dataExtent(rows.map((row) => row.spectral_spread_final)),
  };
}

function getLatticeGeometry(boundaryType) {
  const key = `${boundaryType}|${atlasConfig.latticeSide}`;
  if (state.latticeCache.has(key)) {
    return state.latticeCache.get(key);
  }

  const geometry = buildLatticeGeometry(Number(atlasConfig.latticeSide) || 16, boundaryType);
  state.latticeCache.set(key, geometry);
  return geometry;
}

function buildLatticeGeometry(nSide, boundaryType) {
  const periodic = boundaryType === "periodic";
  const scale = periodic ? nSide : Math.max(nSide - 1, 1);
  const points = [];
  const edges = [];
  const nodeIndex = (i, j, k) => i * nSide * nSide + j * nSide + k;

  for (let i = 0; i < nSide; i += 1) {
    for (let j = 0; j < nSide; j += 1) {
      for (let k = 0; k < nSide; k += 1) {
        points.push([i / scale, j / scale, k / scale]);
      }
    }
  }

  for (let i = 0; i < nSide; i += 1) {
    for (let j = 0; j < nSide; j += 1) {
      for (let k = 0; k < nSide; k += 1) {
        addEdge("x", i, j, k, i + 1, j, k);
        addEdge("y", i, j, k, i, j + 1, k);
        addEdge("z", i, j, k, i, j, k + 1);
      }
    }
  }

  function addEdge(axis, i, j, k, ni, nj, nk) {
    if (!periodic && (ni >= nSide || nj >= nSide || nk >= nSide)) {
      return;
    }

    const wrap = periodic && (ni >= nSide || nj >= nSide || nk >= nSide);
    const wrappedI = periodic ? ni % nSide : ni;
    const wrappedJ = periodic ? nj % nSide : nj;
    const wrappedK = periodic ? nk % nSide : nk;
    edges.push({
      axis,
      start: nodeIndex(i, j, k),
      end: nodeIndex(wrappedI, wrappedJ, wrappedK),
      wrap,
    });
  }

  return { nSide, boundaryType, points, edges };
}

function createTransitionStats(rows) {
  const baselineRegimes = orderedCategories(rows.map((row) => row.baseline_regime));
  const perturbedRegimes = orderedCategories(rows.map((row) => row.perturbed_regime));
  const matrix = new Map();

  rows.forEach((row) => {
    const key = `${row.baseline_regime}|||${row.perturbed_regime}`;
    matrix.set(key, (matrix.get(key) || 0) + 1);
  });

  const persistenceByBaseline = baselineRegimes.map((regime) => {
    const matches = rows.filter((row) => row.baseline_regime === regime);
    const mean =
      matches.reduce((sum, row) => sum + row.persistence_score, 0) / Math.max(matches.length, 1);
    return {
      regime,
      meanPersistence: mean,
      count: matches.length,
    };
  });

  const stabilityByRun = Array.from(groupBy(rows, "baseline_run_id").entries())
    .map(([baselineRunId, items]) => {
      const meanPersistence =
        items.reduce((sum, item) => sum + item.persistence_score, 0) / Math.max(items.length, 1);
      const baselineRow = state.atlas0.rowById.get(baselineRunId);
      return {
        baselineRunId,
        meanPersistence,
        count: items.length,
        regimeLabel: baselineRow?.regime_label || items[0]?.baseline_regime || "Unknown",
        sectorLabel: baselineRow?.sector_label || items[0]?.baseline_sector_label || "Unknown",
      };
    })
    .sort((left, right) => {
      if (right.meanPersistence !== left.meanPersistence) {
        return right.meanPersistence - left.meanPersistence;
      }
      return left.baselineRunId.localeCompare(right.baselineRunId);
    });

  return {
    totalTransitions: rows.length,
    persistentTransitions: rows.filter((row) => row.persistence_score >= 1).length,
    baselineRegimes,
    perturbedRegimes,
    matrix,
    persistenceByBaseline,
    stabilityByRun,
  };
}

function populateFilterOptions() {
  populateSelect(elements.regimeFilter, uniqueSorted(state.atlas0.rows.map((row) => row.regime_label)), "All regimes");
  populateSelect(elements.sectorFilter, uniqueSorted(state.atlas0.rows.map((row) => row.sector_label)), "All sectors");
  populateSelect(
    elements.boundaryFilter,
    uniqueSorted(state.atlas0.rows.map((row) => row.boundary_type)),
    "All boundaries",
  );
}

function populateSelect(element, values, label) {
  element.innerHTML = "";
  const defaultOption = document.createElement("option");
  defaultOption.value = "all";
  defaultOption.textContent = label;
  element.appendChild(defaultOption);

  values.forEach((value) => {
    const option = document.createElement("option");
    option.value = value;
    option.textContent = value;
    element.appendChild(option);
  });
}

function applyFilters() {
  const filtered = state.atlas0.rows.filter((row) => {
    const matchesSearch =
      !state.filters.search || row.run_id.toLowerCase().includes(state.filters.search);
    const matchesRegime =
      state.filters.regime === "all" || row.regime_label === state.filters.regime;
    const matchesSector =
      state.filters.sector === "all" || row.sector_label === state.filters.sector;
    const matchesBoundary =
      state.filters.boundary === "all" || row.boundary_type === state.filters.boundary;
    return matchesSearch && matchesRegime && matchesSector && matchesBoundary;
  });

  filtered.sort((left, right) => compareRows(left, right, state.sort.key, state.sort.direction));
  state.filteredRows = filtered;

  if (!filtered.some((row) => row.run_id === state.selectedRunId)) {
    state.selectedRunId = filtered[0]?.run_id || null;
  }

  renderRunTable();
  renderRegimeMap();
  renderSummary();
  renderSelection();
  renderParameterSpace();
  renderTrajectoryPlayback();
  void renderRunDetail();
}

function renderStaticPanels() {
  renderDataSources();
  renderRunTable();
  renderRegimeMap();
  renderParameterSpace();
  renderTrajectoryPlayback();
}

function renderDataSources() {
  elements.dataSourceList.innerHTML = [
    ["Atlas-0 CSV", atlasConfig.csvPath],
    ["Atlas-0 JSON", atlasConfig.jsonPath],
    ["Atlas-1 CSV", atlasConfig.transitionCsvPath],
    ["Atlas-1 JSON", atlasConfig.transitionJsonPath],
    ["Plots", atlasConfig.plotBasePath],
    ["Atlas-0 rows", state.atlas0.rows.length ? `${state.atlas0.rows.length} loaded` : "Pending"],
    ["Atlas-1 rows", state.atlas1.loaded ? `${state.atlas1.rows.length} loaded` : "Optional"],
  ]
    .map(
      ([label, value]) =>
        `<dt>${escapeHtml(label)}</dt><dd>${escapeHtml(String(value))}</dd>`,
    )
    .join("");
}

function renderSummary() {
  const stats = state.atlas0.stats;
  if (!stats) {
    elements.summaryChips.innerHTML = "";
    return;
  }

  const chips = [
    {
      label: "Baseline runs",
      value: String(stats.totalRuns),
    },
    {
      label: "Top regime",
      value: stats.topRegime ? `${stats.topRegime[0]} (${stats.topRegime[1]})` : "N/A",
    },
    {
      label: "Top sector",
      value: stats.topSector ? `${stats.topSector[0]} (${stats.topSector[1]})` : "N/A",
    },
    {
      label: "Mean coherence",
      value: formatNumber(stats.meanCoherence),
    },
  ];

  if (state.atlas1.stats?.stabilityByRun?.length) {
    const topStable = state.atlas1.stats.stabilityByRun[0];
    chips.push({
      label: "Most stable seed",
      value: `${topStable.baselineRunId} (${formatPercent(topStable.meanPersistence)})`,
    });
  }

  elements.summaryChips.innerHTML = chips
    .map(
      (chip) =>
        `<div class="summary-chip"><strong>${escapeHtml(chip.label)}</strong><span>${escapeHtml(chip.value)}</span></div>`,
    )
    .join("");
}

function renderSelection() {
  const run = getSelectedRun();
  if (!run) {
    elements.selectionBrief.textContent = "No run selected.";
    return;
  }

  const transitions = state.atlas1.transitionsByBaseline.get(run.run_id) || [];
  const meanPersistence = transitions.length
    ? transitions.reduce((sum, item) => sum + item.persistence_score, 0) / transitions.length
    : null;
  const persistenceText =
    meanPersistence == null ? "No Atlas-1 perturbations loaded." : `Mean persistence ${formatPercent(meanPersistence)} across ${transitions.length} perturbations.`;

  elements.selectionBrief.innerHTML = `
    <div class="pill"><span class="mono">${escapeHtml(run.run_id)}</span></div>
    <p><strong>${escapeHtml(run.regime_label)}</strong> in ${escapeHtml(run.sector_label)}.</p>
    <p>coherence_score ${formatNumber(run.coherence_score)} | anisotropy_score ${formatNumber(run.anisotropy_score)} | proto_spacetime_score ${formatInteger(run.proto_spacetime_score)}</p>
    <p>${escapeHtml(persistenceText)}</p>
  `;
}

function renderRunTable() {
  const rowCountText = `${state.filteredRows.length} of ${state.atlas0.rows.length} runs shown`;
  elements.resultCount.textContent = rowCountText;
  elements.sortIndicator.textContent = `Sorted by ${state.sort.key} ${state.sort.direction === "asc" ? "ascending" : "descending"}`;

  document.querySelectorAll(".sort-button").forEach((button) => {
    const isActive = button.dataset.sortKey === state.sort.key;
    button.classList.toggle("active", isActive);
    button.textContent = `${button.dataset.sortKey}${isActive ? state.sort.direction === "asc" ? " ↑" : " ↓" : ""}`;
  });

  if (!state.filteredRows.length) {
    elements.runTableBody.innerHTML = `
      <tr>
        <td colspan="${TABLE_COLUMNS.length}" class="faint">No runs match the current filters.</td>
      </tr>
    `;
    return;
  }

  elements.runTableBody.innerHTML = state.filteredRows
    .map((row) => {
      const selectedClass = row.run_id === state.selectedRunId ? "selected" : "";
      return `
        <tr class="${selectedClass}" data-run-id="${escapeHtml(row.run_id)}">
          <td class="mono">${escapeHtml(row.run_id)}</td>
          <td>${escapeHtml(row.atlas_phase)}</td>
          <td>${escapeHtml(row.regime_label)}</td>
          <td>${escapeHtml(row.sector_label)}</td>
          <td>${formatInteger(row.proto_spacetime_score)}</td>
          <td>${formatNumber(row.coherence_score)}</td>
          <td>${formatNumber(row.anisotropy_score)}</td>
          <td>${escapeHtml(row.boundary_type)}</td>
          <td>${escapeHtml(row.operator_sector)}</td>
        </tr>
      `;
    })
    .join("");
}

function renderRegimeMap() {
  const svg = elements.regimeMap;
  clearSvg(svg);

  if (!state.atlas0.rows.length) {
    return;
  }

  const width = 640;
  const height = 420;
  const margin = { top: 26, right: 22, bottom: 58, left: 64 };
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;
  const domainRows = state.atlas0.rows;
  const drawRows = state.filteredRows;
  const xValues = domainRows.map((row) => row.anisotropy_score);
  const yValues = domainRows.map((row) => row.coherence_score);
  const xDomain = paddedExtent(xValues);
  const yDomain = paddedExtent(yValues);

  drawLinearGrid(svg, {
    width,
    height,
    margin,
    xDomain,
    yDomain,
    xLabel: "anisotropy_score",
    yLabel: "coherence_score",
  });

  drawRows.forEach((row) => {
    const x = scaleLinear(row.anisotropy_score, xDomain[0], xDomain[1], margin.left, width - margin.right);
    const y = scaleLinear(row.coherence_score, yDomain[0], yDomain[1], height - margin.bottom, margin.top);
    const shape = createSectorShape(row.sector_label, x, y, 7);
    shape.setAttribute("fill", regimeColor(row.regime_label));
    shape.setAttribute("opacity", row.run_id === state.selectedRunId ? "1" : "0.86");
    shape.setAttribute("data-run-id", row.run_id);
    shape.classList.add("map-point");
    if (row.run_id === state.selectedRunId) {
      shape.classList.add("is-selected");
    }

    const tooltip = document.createElementNS("http://www.w3.org/2000/svg", "title");
    tooltip.textContent = `${row.run_id}\n${row.regime_label}\n${row.sector_label}\ncoherence ${formatNumber(row.coherence_score)} | anisotropy ${formatNumber(row.anisotropy_score)}`;
    shape.appendChild(tooltip);
    svg.appendChild(shape);
  });

  renderMapLegend();
}

function renderMapLegend() {
  const regimeItems = orderedCategories(state.atlas0.rows.map((row) => row.regime_label))
    .map(
      (regime) =>
        `<div class="legend-item"><span class="legend-swatch" style="background:${regimeColor(regime)}"></span>${escapeHtml(regime)}</div>`,
    )
    .join("");

  const sectorItems = uniqueSorted(state.atlas0.rows.map((row) => row.sector_label))
    .map(
      (sector) =>
        `<div class="legend-item"><span class="legend-shape">${renderShapeGlyph(sector)}</span>${escapeHtml(sector)}</div>`,
    )
    .join("");

  elements.mapLegend.innerHTML = `
    <div>
      <h3>Color by regime_label</h3>
      <div class="legend-section">${regimeItems}</div>
    </div>
    <div>
      <h3>Shape by sector_label</h3>
      <div class="legend-section">${sectorItems}</div>
    </div>
  `;
}

function renderParameterSpace() {
  const canvas = elements.parameterSpaceCanvas;
  if (!canvas) {
    return;
  }

  const { ctx, width, height } = prepareCanvas(canvas);
  if (!ctx || !state.atlas0.rows.length) {
    elements.parameterSpaceInfo.textContent = "No Atlas-0 parameter data available.";
    return;
  }

  const visibleRows = state.filteredRows.length ? state.filteredRows : state.atlas0.rows;
  const extents = {
    x: paddedExtent(state.atlas0.rows.map((row) => row.central_k)),
    y: paddedExtent(state.atlas0.rows.map((row) => row.bandwidth)),
    z: paddedExtent(state.atlas0.rows.map((row) => row.amplitude)),
  };

  draw3dScaffold(ctx, width, height, state.view3d.parameter, {
    x: "central_k",
    y: "bandwidth",
    z: "amplitude",
  });

  const projected = visibleRows
    .map((row) => {
      const proxy = waveParticleProxy(row);
      const point = {
        x: normalizeToCube(row.central_k, extents.x),
        y: normalizeToCube(row.bandwidth, extents.y),
        z: normalizeToCube(row.amplitude, extents.z),
      };
      const projection = projectScenePoint(point, state.view3d.parameter, width, height, 0.9);
      return { row, proxy, point, ...projection };
    })
    .sort((left, right) => left.depth - right.depth);

  projected.forEach((item) => {
    const isSelected = item.row.run_id === state.selectedRunId;
    const isHovered = item.row.run_id === state.view3d.parameter.hoverRunId;
    const fill =
      state.view3d.parameter.colorMode === "wave-particle"
        ? waveParticleColor(item.proxy)
        : regimeColor(item.row.regime_label);
    drawCanvasShape(ctx, item.row.sector_label, item.x, item.y, pointRadius(item, 7), fill, isSelected, isHovered);
  });

  state.view3d.parameter.projected = projected;

  const focusRun =
    (state.view3d.parameter.hoverRunId && state.atlas0.rowById.get(state.view3d.parameter.hoverRunId)) ||
    getSelectedRun();

  if (!focusRun) {
    elements.parameterSpaceInfo.textContent = "Drag to rotate. Click a point to inspect a run.";
    return;
  }

  const focusProxy = waveParticleProxy(focusRun);
  const colorModeLabel =
    state.view3d.parameter.colorMode === "wave-particle" ? "wave / particle proxy" : "regime";

  elements.parameterSpaceInfo.innerHTML = `
    <strong class="mono">${escapeHtml(focusRun.run_id)}</strong><br>
    ${escapeHtml(focusRun.regime_label)} | ${escapeHtml(focusRun.sector_label)} | ${escapeHtml(propagationClassForRun(focusRun))}<br>
    central_k ${formatNumber(focusRun.central_k)} | bandwidth ${formatNumber(focusRun.bandwidth)} | amplitude ${formatNumber(focusRun.amplitude)}<br>
    coherence_score ${formatNumber(focusRun.coherence_score)} | anisotropy_score ${formatNumber(focusRun.anisotropy_score)}<br>
    color mode ${escapeHtml(colorModeLabel)} | wave proxy ${formatNumber(focusProxy.wave)} | particle proxy ${formatNumber(focusProxy.particle)} | ${escapeHtml(focusProxy.label)}
  `;
}

function getThreeTrajectoryRenderer() {
  if (atlasConfig.trajectoryRenderer === "canvas") {
    return null;
  }

  const renderer = window.AtlasThreeGridRenderer;
  if (!renderer || renderer.available === false || typeof renderer.render !== "function") {
    return null;
  }

  return renderer;
}

function clearTrajectorySurface() {
  const canvas = elements.trajectoryCanvas;
  if (!canvas) {
    return;
  }

  const renderer = getThreeTrajectoryRenderer();
  if (renderer && typeof renderer.clear === "function") {
    renderer.clear(canvas);
    return;
  }

  const { ctx, width, height } = prepareCanvas(canvas);
  if (ctx) {
    ctx.clearRect(0, 0, width, height);
  }
}

function centerToScenePoint(center) {
  return {
    x: toNumber(center?.[0]) - 0.5,
    y: toNumber(center?.[1]) - 0.5,
    z: toNumber(center?.[2]) - 0.5,
  };
}

function buildTrajectoryFrame(run) {
  const centers = Array.isArray(run.metrics?.centers) ? run.metrics.centers : [];
  const times = Array.isArray(run.metrics?.times) ? run.metrics.times : [];
  const widths = Array.isArray(run.metrics?.widths) ? run.metrics.widths : [];
  const energies = Array.isArray(run.metrics?.energies) ? run.metrics.energies : [];
  const coherences = Array.isArray(run.metrics?.coherences) ? run.metrics.coherences : [];

  if (!centers.length) {
    return null;
  }

  if (state.view3d.trajectory.lastRunId !== run.run_id) {
    state.view3d.trajectory.lastRunId = run.run_id;
    state.view3d.trajectory.frameIndex = Math.max(centers.length - 1, 0);
  }

  const frameIndex = clamp(state.view3d.trajectory.frameIndex, 0, centers.length - 1);
  state.view3d.trajectory.frameIndex = frameIndex;
  elements.trajectoryScrub.max = String(Math.max(centers.length - 1, 0));
  elements.trajectoryScrub.value = String(frameIndex);

  const lattice = getLatticeGeometry(run.boundary_type);
  const currentCenter = centers[frameIndex];
  const previousCenter = centers[Math.max(frameIndex - 1, 0)] || currentCenter;
  const currentWidth = widths[frameIndex] ?? widths[widths.length - 1] ?? 0;
  const currentTime = times[frameIndex] ?? frameIndex;
  const currentCoherence = coherences[frameIndex] ?? run.coherence_score;
  const currentEnergy = energies[frameIndex] ?? energies[energies.length - 1] ?? 0;
  const maxEnergy = Math.max(currentEnergy, ...energies, Number.EPSILON);
  const sigma = flowSigmaFromWidth(currentWidth);
  const velocity = wrappedVelocityVector(previousCenter, currentCenter, run.boundary_type);
  const speed = Math.hypot(velocity[0], velocity[1], velocity[2]);
  const velocityUnit = speed > 0 ? velocity.map((value) => value / speed) : [1, 0, 0];
  const energyNorm = clamp(currentEnergy / maxEnergy, 0, 1);
  const pulsePhaseSeed = ((Date.now() / 420) + frameIndex * 0.19) % 1;
  const regimeFill = regimeColor(run.regime_label);

  const centerTrail = centers.map((center, index) => ({
    index,
    point: centerToScenePoint(center),
    widthValue: widths[index] ?? widths[widths.length - 1] ?? 0.08,
    coherenceValue: coherences[index] ?? run.coherence_score,
  }));

  const nodeIntensities = new Array(lattice.points.length);
  const nodes = lattice.points.map((rawPoint, index) => {
    const delta = wrappedDisplacement(rawPoint, currentCenter, run.boundary_type);
    const dist2 = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
    const intensity =
      Math.exp(-dist2 / Math.max(2 * sigma * sigma, Number.EPSILON)) *
      (0.34 + energyNorm * 0.46) *
      (0.42 + clamp(run.amplitude, 0, 1.2) * 0.58);

    nodeIntensities[index] = intensity;

    return {
      index,
      point: {
        x: rawPoint[0] - 0.5,
        y: rawPoint[1] - 0.5,
        z: rawPoint[2] - 0.5,
      },
      intensity,
      sampleVisible: intensity > 0.03 || index % 7 === 0,
    };
  });

  const activeEdges = [];
  lattice.edges.forEach((edge, edgeIndex) => {
    if (edge.wrap) {
      return;
    }

    const meanIntensity = (nodeIntensities[edge.start] + nodeIntensities[edge.end]) * 0.5;
    const signedAlignment = velocityUnit[AXIS_INDEX[edge.axis]] || 0;
    const flowStrength = meanIntensity * (0.22 + 0.78 * Math.abs(signedAlignment));
    if (flowStrength < 0.025) {
      return;
    }

    activeEdges.push({
      edgeIndex,
      axis: edge.axis,
      startIndex: edge.start,
      endIndex: edge.end,
      flowStrength,
      signedAlignment,
    });
  });

  return {
    run,
    lattice,
    frameIndex,
    centers,
    centerTrail,
    times,
    widths,
    coherences,
    energies,
    currentCenter,
    currentCenterPoint: centerToScenePoint(currentCenter),
    currentWidth,
    currentTime,
    currentCoherence,
    currentEnergy,
    energyNorm,
    sigma,
    velocity,
    velocityUnit,
    pulsePhaseSeed,
    nodeIntensities,
    nodes,
    activeEdges,
    visibleNodeCount: nodes.filter((item) => item.intensity > 0.03).length,
    regimeFill,
    modeProxy: waveParticleProxy(run),
  };
}

function renderTrajectoryPlaybackCanvas(frame, ctx, width, height) {
  draw3dScaffold(ctx, width, height, state.view3d.trajectory, {
    x: "x",
    y: "y",
    z: "z",
  });

  const projected = frame.centerTrail.map((item) => ({
    ...item,
    ...projectScenePoint(item.point, state.view3d.trajectory, width, height, 0.82),
  }));

  const projectedNodes = frame.nodes.map((item) => ({
    ...item,
    ...projectScenePoint(item.point, state.view3d.trajectory, width, height, 0.84),
  }));

  projectedNodes
    .filter((item) => item.sampleVisible)
    .sort((left, right) => left.depth - right.depth)
    .forEach((item) => {
      const intensity = clamp(item.intensity, 0, 1);
      const radius = intensity > 0.06 ? 0.9 + intensity * 8.6 : 0.7;
      const alpha = intensity > 0.06 ? 0.1 + intensity * 0.7 : 0.05;

      ctx.save();
      ctx.beginPath();
      ctx.fillStyle =
        intensity > 0.08 ? colorWithAlpha(frame.regimeFill, alpha) : `rgba(157, 177, 196, ${alpha})`;
      ctx.arc(item.x, item.y, radius, 0, Math.PI * 2);
      ctx.fill();
      ctx.restore();
    });

  const projectedNodeByIndex = new Map(projectedNodes.map((item) => [item.index, item]));
  const projectedActiveEdges = [];

  frame.activeEdges.forEach((edge) => {
    const from = projectedNodeByIndex.get(edge.startIndex);
    const to = projectedNodeByIndex.get(edge.endIndex);
    if (!from || !to) {
      return;
    }

    projectedActiveEdges.push({
      ...edge,
      from,
      to,
      depth: (from.depth + to.depth) * 0.5,
    });
  });

  projectedActiveEdges
    .sort((left, right) => left.depth - right.depth)
    .forEach((edge) => {
      ctx.save();
      ctx.lineWidth = 0.45 + edge.flowStrength * 5.2;
      ctx.strokeStyle =
        edge.signedAlignment >= 0
          ? colorWithAlpha("#7ea2a8", 0.08 + edge.flowStrength * 0.58)
          : colorWithAlpha("#b49463", 0.08 + edge.flowStrength * 0.58);
      ctx.beginPath();
      ctx.moveTo(edge.from.x, edge.from.y);
      ctx.lineTo(edge.to.x, edge.to.y);
      ctx.stroke();
      ctx.restore();

      const directedPhase = edge.signedAlignment >= 0 ? frame.pulsePhaseSeed : 1 - frame.pulsePhaseSeed;
      const pulseX = edge.from.x + (edge.to.x - edge.from.x) * directedPhase;
      const pulseY = edge.from.y + (edge.to.y - edge.from.y) * directedPhase;
      const pulseRadius = 0.8 + edge.flowStrength * 3.2;

      ctx.save();
      ctx.beginPath();
      ctx.fillStyle =
        edge.signedAlignment >= 0
          ? colorWithAlpha("#f0f7fb", 0.16 + edge.flowStrength * 0.84)
          : colorWithAlpha("#f3d9b0", 0.16 + edge.flowStrength * 0.84);
      ctx.arc(pulseX, pulseY, pulseRadius, 0, Math.PI * 2);
      ctx.fill();
      ctx.restore();
    });

  ctx.save();
  ctx.lineWidth = 1.6;
  ctx.strokeStyle = "rgba(126, 162, 168, 0.52)";
  ctx.beginPath();
  projected.forEach((item, index) => {
    if (index === 0) {
      ctx.moveTo(item.x, item.y);
    } else {
      ctx.lineTo(item.x, item.y);
    }
  });
  ctx.stroke();
  ctx.restore();

  projected.forEach((item) => {
    const isCurrent = item.index === frame.frameIndex;
    const radius = isCurrent ? 4.8 + item.widthValue * 22 : 2;
    const fill = isCurrent ? frame.regimeFill : coherenceColor(item.coherenceValue);
    drawCanvasShape(ctx, frame.run.sector_label, item.x, item.y, radius, fill, isCurrent, false);
  });

  state.view3d.trajectory.projected = projected;
}

function renderTrajectoryPlayback() {
  const canvas = elements.trajectoryCanvas;
  if (!canvas) {
    return;
  }

  const run = getSelectedRun();
  if (!run) {
    stopTrajectoryPlayback();
    clearTrajectorySurface();
    elements.trajectoryFrameInfo.textContent = "Select a run to inspect its lattice replay.";
    return;
  }

  const frame = buildTrajectoryFrame(run);
  if (!frame) {
    stopTrajectoryPlayback();
    clearTrajectorySurface();
    elements.trajectoryFrameInfo.textContent = "This run does not expose recorded center traces in the current JSON payload.";
    return;
  }

  const renderer = getThreeTrajectoryRenderer();
  const usedThree =
    renderer &&
    renderer.render({
      canvas,
      frame,
      viewState: state.view3d.trajectory,
    });

  if (!usedThree) {
    const { ctx, width, height } = prepareCanvas(canvas);
    if (!ctx) {
      return;
    }
    renderTrajectoryPlaybackCanvas(frame, ctx, width, height);
  } else {
    state.view3d.trajectory.projected = [];
  }

  state.view3d.trajectory.visibleNodeCount = frame.visibleNodeCount;
  state.view3d.trajectory.activeEdgeCount = frame.activeEdges.length;
  const boundaryNote =
    run.boundary_type === "periodic"
      ? "Periodic wrap links are omitted from the cube view."
      : "Open-boundary links are shown directly.";
  const latticeLabel = `${atlasConfig.latticeSide}x${atlasConfig.latticeSide}x${atlasConfig.latticeSide}`;
  const rendererLabel = usedThree ? "three.js webgl" : "canvas fallback";

  elements.trajectoryFrameInfo.innerHTML = `
    <strong class="mono">${escapeHtml(run.run_id)}</strong><br>
    lattice ${escapeHtml(latticeLabel)} | ${escapeHtml(run.boundary_type)} boundary | active nodes ${state.view3d.trajectory.visibleNodeCount} | active links ${state.view3d.trajectory.activeEdgeCount}<br>
    frame ${frame.frameIndex + 1} / ${frame.centers.length} | time ${formatNumber(frame.currentTime)} | width ${formatNumber(frame.currentWidth)} | energy ${formatNumber(frame.currentEnergy)} | coherence ${formatNumber(frame.currentCoherence)}<br>
    center (${formatNumber(frame.currentCenter[0])}, ${formatNumber(frame.currentCenter[1])}, ${formatNumber(frame.currentCenter[2])})<br>
    wave proxy ${formatNumber(frame.modeProxy.wave)} | particle proxy ${formatNumber(frame.modeProxy.particle)} | ${escapeHtml(frame.modeProxy.label)}<br>
    renderer ${escapeHtml(rendererLabel)} | Edge flow is a deterministic proxy from the recorded center motion and packet spread. ${escapeHtml(boundaryNote)}
  `;
}

function toggleTrajectoryPlayback() {
  if (state.view3d.trajectory.playIntervalId) {
    stopTrajectoryPlayback();
    return;
  }

  const run = getSelectedRun();
  const centers = Array.isArray(run?.metrics?.centers) ? run.metrics.centers : [];
  if (!run || !centers.length) {
    renderTrajectoryPlayback();
    return;
  }

  startTrajectoryPlayback(centers.length);
}

function startTrajectoryPlayback(frameCount) {
  stopTrajectoryPlayback(false);
  elements.trajectoryPlayToggle.textContent = "Pause";

  state.view3d.trajectory.playIntervalId = window.setInterval(() => {
    if (frameCount <= 1) {
      stopTrajectoryPlayback();
      return;
    }
    state.view3d.trajectory.frameIndex =
      state.view3d.trajectory.frameIndex >= frameCount - 1 ? 0 : state.view3d.trajectory.frameIndex + 1;
    renderTrajectoryPlayback();
  }, 220);
}

function stopTrajectoryPlayback(updateLabel = true) {
  if (state.view3d.trajectory.playIntervalId) {
    window.clearInterval(state.view3d.trajectory.playIntervalId);
    state.view3d.trajectory.playIntervalId = null;
  }

  if (updateLabel && elements.trajectoryPlayToggle) {
    elements.trajectoryPlayToggle.textContent = "Play";
  }
}

function bind3dCanvasInteractions() {
  bindRotationCanvas(elements.parameterSpaceCanvas, state.view3d.parameter, {
    render: renderParameterSpace,
    onMove(x, y) {
      const nearest = findNearestProjectedPoint(state.view3d.parameter.projected, x, y, 18);
      const nextRunId = nearest ? nearest.row.run_id : null;
      if (nextRunId !== state.view3d.parameter.hoverRunId) {
        state.view3d.parameter.hoverRunId = nextRunId;
        renderParameterSpace();
      }
    },
    onLeave() {
      if (state.view3d.parameter.hoverRunId) {
        state.view3d.parameter.hoverRunId = null;
        renderParameterSpace();
      }
    },
    onClick(x, y) {
      const nearest = findNearestProjectedPoint(state.view3d.parameter.projected, x, y, 18);
      if (nearest) {
        selectRun(nearest.row.run_id);
      }
    },
  });

  bindRotationCanvas(elements.trajectoryCanvas, state.view3d.trajectory, {
    render: renderTrajectoryPlayback,
  });
}

function bindRotationCanvas(canvas, targetState, handlers) {
  if (!canvas) {
    return;
  }

  canvas.addEventListener("pointerdown", (event) => {
    const position = canvasPointerPosition(canvas, event);
    targetState.drag = {
      pointerId: event.pointerId,
      startX: position.x,
      startY: position.y,
      rotationX: targetState.rotationX,
      rotationY: targetState.rotationY,
      moved: false,
    };
    canvas.setPointerCapture(event.pointerId);
  });

  canvas.addEventListener("pointermove", (event) => {
    const position = canvasPointerPosition(canvas, event);
    if (targetState.drag) {
      const dx = position.x - targetState.drag.startX;
      const dy = position.y - targetState.drag.startY;
      if (Math.abs(dx) > 1 || Math.abs(dy) > 1) {
        targetState.drag.moved = true;
      }
      targetState.rotationY = targetState.drag.rotationY + dx * 0.012;
      targetState.rotationX = clamp(targetState.drag.rotationX + dy * 0.012, -1.35, 1.35);
      handlers.render();
      return;
    }

    if (handlers.onMove) {
      handlers.onMove(position.x, position.y);
    }
  });

  canvas.addEventListener("pointerup", (event) => {
    const position = canvasPointerPosition(canvas, event);
    const drag = targetState.drag;
    targetState.drag = null;
    if (!drag) {
      return;
    }
    if (!drag.moved && handlers.onClick) {
      handlers.onClick(position.x, position.y);
    }
  });

  canvas.addEventListener("pointerleave", () => {
    if (!targetState.drag && handlers.onLeave) {
      handlers.onLeave();
    }
  });
}

async function renderRunDetail() {
  const run = getSelectedRun();
  if (!run) {
    elements.detailEmpty.classList.remove("hidden");
    elements.detailContent.classList.add("hidden");
    return;
  }

  elements.detailEmpty.classList.add("hidden");
  elements.detailContent.classList.remove("hidden");

  elements.detailMetrics.innerHTML = [
    ["regime_label", run.regime_label],
    ["sector_label", run.sector_label],
    ["proto_spacetime_score", formatInteger(run.proto_spacetime_score)],
    ["coherence_score", formatNumber(run.coherence_score)],
    ["anisotropy_score", formatNumber(run.anisotropy_score)],
    ["boundary_type", run.boundary_type],
  ]
    .map(
      ([label, value]) =>
        `<div class="metric-card"><strong>${escapeHtml(label)}</strong><span>${escapeHtml(String(value))}</span></div>`,
    )
    .join("");

  elements.detailMetadata.innerHTML = [
    ["run_id", run.run_id],
    ["atlas_phase", run.atlas_phase],
    ["operator_sector", run.operator_sector],
    ["initial_seed_type", run.initial_seed_type],
    ["central_k", formatNumber(run.central_k)],
    ["bandwidth", formatNumber(run.bandwidth)],
    ["amplitude", formatNumber(run.amplitude)],
    ["phase_pattern", run.phase_pattern],
    ["random_seed", formatInteger(run.random_seed)],
    ["center_shift", formatNumber(run.center_shift)],
    ["width_ratio", formatNumber(run.width_ratio)],
    ["notes", run.notes || "N/A"],
  ]
    .map(
      ([label, value]) =>
        `<dt>${escapeHtml(label)}</dt><dd class="${label === "run_id" ? "mono" : ""}">${escapeHtml(String(value))}</dd>`,
    )
    .join("");

  renderTransitionDetail(run);

  const token = ++state.detailToken;
  const cached = state.plotResolutionCache.get(run.run_id);
  if (cached) {
    renderImageCards(cached);
    return;
  }

  renderImageCards(
    PLOT_DEFINITIONS.map((plot) => ({
      ...plot,
      url: null,
      loading: true,
    })),
  );

  const resolved = await Promise.all(
    PLOT_DEFINITIONS.map(async (plot) => ({
      ...plot,
      url: await pickExistingUrl(run.plotCandidates[plot.key]),
      loading: false,
    })),
  );

  if (token !== state.detailToken) {
    return;
  }

  state.plotResolutionCache.set(run.run_id, resolved);
  renderImageCards(resolved);
}

function renderTransitionDetail(run) {
  const transitions = state.atlas1.transitionsByBaseline.get(run.run_id) || [];
  if (!transitions.length) {
    elements.detailTransitions.innerHTML = `<p class="faint">No Atlas-1 perturbation rows loaded for this baseline seed.</p>`;
    return;
  }

  const meanPersistence =
    transitions.reduce((sum, item) => sum + item.persistence_score, 0) / transitions.length;
  const stableCount = transitions.filter((item) => item.persistence_score >= 1).length;
  const rows = [...transitions].sort((left, right) => left.perturbation_axis.localeCompare(right.perturbation_axis));

  elements.detailTransitions.innerHTML = `
    <div class="summary-chips">
      <div class="summary-chip"><strong>Mean persistence</strong><span>${formatPercent(meanPersistence)}</span></div>
      <div class="summary-chip"><strong>Stable transitions</strong><span>${stableCount} / ${transitions.length}</span></div>
    </div>
    <table class="mini-table">
      <thead>
        <tr>
          <th>perturbation_axis</th>
          <th>perturbed_regime</th>
          <th>persistence_score</th>
          <th>coherence_delta</th>
        </tr>
      </thead>
      <tbody>
        ${rows
          .map(
            (item) => `
              <tr>
                <td class="mono">${escapeHtml(item.perturbation_axis)}</td>
                <td>${escapeHtml(item.perturbed_regime)}</td>
                <td>${formatInteger(item.persistence_score)}</td>
                <td>${formatSignedNumber(item.coherence_delta)}</td>
              </tr>
            `,
          )
          .join("")}
      </tbody>
    </table>
  `;
}

function renderImageCards(items) {
  elements.detailImages.innerHTML = items
    .map((item) => {
      const body = item.loading
        ? `<div class="image-frame">Resolving plot path</div>`
        : item.url
          ? `<div class="image-frame"><img src="${escapeHtml(item.url)}" alt="${escapeHtml(item.label)}"></div>`
          : `<div class="image-frame">Artifact not found in plotBasePath.</div>`;
      return `
        <article class="image-card">
          <header><h3>${escapeHtml(item.label)}</h3></header>
          ${body}
        </article>
      `;
    })
    .join("");
}

function renderTransitionPanel() {
  if (!state.atlas1.loaded || !state.atlas1.stats) {
    elements.transitionEmpty.classList.remove("hidden");
    elements.transitionContent.classList.add("hidden");
    return;
  }

  elements.transitionEmpty.classList.add("hidden");
  elements.transitionContent.classList.remove("hidden");
  drawTransitionMatrix();
  drawPersistenceBar();
  renderTransitionSummary();
}

function drawTransitionMatrix() {
  const svg = elements.transitionMatrix;
  clearSvg(svg);

  const width = 680;
  const height = 420;
  const margin = { top: 72, right: 24, bottom: 110, left: 176 };
  const baselineRegimes = state.atlas1.stats.baselineRegimes;
  const perturbedRegimes = state.atlas1.stats.perturbedRegimes;
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;
  const cellWidth = innerWidth / Math.max(perturbedRegimes.length, 1);
  const cellHeight = innerHeight / Math.max(baselineRegimes.length, 1);
  const maxValue = Math.max(
    ...baselineRegimes.flatMap((baseline) =>
      perturbedRegimes.map((perturbed) => state.atlas1.stats.matrix.get(`${baseline}|||${perturbed}`) || 0),
    ),
    1,
  );

  baselineRegimes.forEach((baseline, rowIndex) => {
    const y = margin.top + rowIndex * cellHeight;
    const label = svgText(12, y + cellHeight * 0.62, baseline, "tick-text");
    label.setAttribute("font-size", "12");
    svg.appendChild(label);

    perturbedRegimes.forEach((perturbed, colIndex) => {
      const x = margin.left + colIndex * cellWidth;
      const value = state.atlas1.stats.matrix.get(`${baseline}|||${perturbed}`) || 0;
      const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
      rect.setAttribute("x", x);
      rect.setAttribute("y", y);
      rect.setAttribute("width", cellWidth - 4);
      rect.setAttribute("height", cellHeight - 4);
      rect.setAttribute("rx", "8");
      rect.setAttribute("fill", heatColor(value, maxValue));
      rect.setAttribute("stroke", "rgba(230, 237, 244, 0.12)");
      svg.appendChild(rect);

      const countLabel = svgText(x + cellWidth / 2 - 2, y + cellHeight / 2 + 4, String(value), "viz-label");
      countLabel.setAttribute("font-size", "12");
      countLabel.setAttribute("text-anchor", "middle");
      countLabel.setAttribute("fill", value > maxValue * 0.55 ? "#081016" : "#dce5ed");
      svg.appendChild(countLabel);
    });
  });

  perturbedRegimes.forEach((perturbed, index) => {
    const x = margin.left + index * cellWidth + cellWidth * 0.5;
    const label = svgText(x, height - 28, perturbed, "tick-text");
    label.setAttribute("font-size", "12");
    label.setAttribute("text-anchor", "end");
    label.setAttribute("transform", `rotate(-35 ${x} ${height - 28})`);
    svg.appendChild(label);
  });

  svg.appendChild(svgText(width / 2 - 60, 24, "perturbed_regime", "axis-text"));
  const baselineLabel = svgText(26, height / 2 + 20, "baseline_regime", "axis-text");
  baselineLabel.setAttribute("transform", `rotate(-90 26 ${height / 2 + 20})`);
  svg.appendChild(baselineLabel);
}

function drawPersistenceBar() {
  const svg = elements.persistenceBar;
  clearSvg(svg);

  const width = 560;
  const height = 420;
  const margin = { top: 24, right: 18, bottom: 110, left: 56 };
  const items = state.atlas1.stats.persistenceByBaseline;
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;
  const barWidth = innerWidth / Math.max(items.length, 1);

  drawBarGrid(svg, width, height, margin);

  items.forEach((item, index) => {
    const x = margin.left + index * barWidth + 10;
    const y = scaleLinear(item.meanPersistence, 0, 1, height - margin.bottom, margin.top);
    const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
    rect.setAttribute("x", x);
    rect.setAttribute("y", y);
    rect.setAttribute("width", Math.max(barWidth - 18, 12));
    rect.setAttribute("height", height - margin.bottom - y);
    rect.setAttribute("rx", "10");
    rect.setAttribute("fill", regimeColor(item.regime));
    rect.setAttribute("opacity", "0.92");
    svg.appendChild(rect);

    const valueLabel = svgText(x + 8, Math.max(y - 8, 18), formatPercent(item.meanPersistence), "viz-label");
    valueLabel.setAttribute("font-size", "12");
    svg.appendChild(valueLabel);

    const regimeLabel = svgText(x + Math.max(barWidth - 18, 12) / 2, height - 18, item.regime, "tick-text");
    regimeLabel.setAttribute("font-size", "12");
    regimeLabel.setAttribute("text-anchor", "end");
    regimeLabel.setAttribute("transform", `rotate(-35 ${x + Math.max(barWidth - 18, 12) / 2} ${height - 18})`);
    svg.appendChild(regimeLabel);
  });
}

function renderTransitionSummary() {
  const stats = state.atlas1.stats;
  const topRuns = stats.stabilityByRun.slice(0, 3);
  const chips = [
    {
      label: "Transition rows",
      value: String(stats.totalTransitions),
    },
    {
      label: "Persistent fraction",
      value: formatPercent(stats.persistentTransitions / Math.max(stats.totalTransitions, 1)),
    },
  ];

  topRuns.forEach((item, index) => {
    chips.push({
      label: `Stable seed ${index + 1}`,
      value: `${item.baselineRunId} (${formatPercent(item.meanPersistence)})`,
    });
  });

  elements.transitionSummary.innerHTML = chips
    .map(
      (chip) =>
        `<div class="summary-chip"><strong>${escapeHtml(chip.label)}</strong><span>${escapeHtml(chip.value)}</span></div>`,
    )
    .join("");
}

function selectRun(runId) {
  stopTrajectoryPlayback();
  state.selectedRunId = runId;
  renderRunTable();
  renderRegimeMap();
  renderSelection();
  renderParameterSpace();
  renderTrajectoryPlayback();
  void renderRunDetail();
}

function getSelectedRun() {
  return state.selectedRunId ? state.atlas0.rowById.get(state.selectedRunId) || null : null;
}

function buildAtlas0PlotCandidates(runId) {
  const prefix = state.atlas0.timestamp ? `${state.atlas0.timestamp}_` : "";
  return Object.fromEntries(
    PLOT_DEFINITIONS.map((plot) => [
      plot.key,
      [
        toAbsoluteAssetUrl(`${prefix}stage10_${runId}_${plot.key}.png`),
        toAbsoluteAssetUrl(`stage10_${runId}_${plot.key}.png`),
      ],
    ]),
  );
}

async function pickExistingUrl(candidates) {
  for (const candidate of candidates) {
    if (await probeImage(candidate)) {
      return candidate;
    }
  }
  return null;
}

function probeImage(url) {
  if (state.assetProbeCache.has(url)) {
    return state.assetProbeCache.get(url);
  }

  const promise = new Promise((resolve) => {
    const image = new Image();
    image.onload = () => resolve(true);
    image.onerror = () => resolve(false);
    image.src = url;
  });
  state.assetProbeCache.set(url, promise);
  return promise;
}

function toAbsoluteAssetUrl(filename) {
  const base = atlasConfig.plotBasePath.endsWith("/") ? atlasConfig.plotBasePath : `${atlasConfig.plotBasePath}/`;
  return new URL(`${base}${filename}`, document.baseURI).toString();
}

async function fetchText(path) {
  const response = await fetch(path, { cache: "no-store" });
  if (!response.ok) {
    throw new Error(`Failed to load ${path}: ${response.status}`);
  }
  return response.text();
}

async function fetchJson(path) {
  const response = await fetch(path, { cache: "no-store" });
  if (!response.ok) {
    throw new Error(`Failed to load ${path}: ${response.status}`);
  }
  return response.json();
}

function parseCsv(text) {
  text = text.replace(/^\uFEFF/, "");
  const rows = [];
  let current = "";
  let row = [];
  let inQuotes = false;

  for (let index = 0; index < text.length; index += 1) {
    const char = text[index];
    const next = text[index + 1];

    if (char === '"') {
      if (inQuotes && next === '"') {
        current += '"';
        index += 1;
      } else {
        inQuotes = !inQuotes;
      }
      continue;
    }

    if (char === "," && !inQuotes) {
      row.push(current);
      current = "";
      continue;
    }

    if ((char === "\n" || char === "\r") && !inQuotes) {
      if (char === "\r" && next === "\n") {
        index += 1;
      }
      row.push(current);
      current = "";
      if (row.some((cell) => cell !== "")) {
        rows.push(row);
      }
      row = [];
      continue;
    }

    current += char;
  }

  if (current.length || row.length) {
    row.push(current);
    rows.push(row);
  }

  const [header = [], ...body] = rows;
  return body.map((cells) => {
    const entry = {};
    header.forEach((column, index) => {
      entry[column] = cells[index] ?? "";
    });
    return entry;
  });
}

function countBy(rows, key) {
  return rows.reduce((counts, row) => {
    const value = row[key];
    counts[value] = (counts[value] || 0) + 1;
    return counts;
  }, {});
}

function groupBy(rows, key) {
  const groups = new Map();
  rows.forEach((row) => {
    const value = row[key];
    if (!groups.has(value)) {
      groups.set(value, []);
    }
    groups.get(value).push(row);
  });
  return groups;
}

function compareRows(left, right, key, direction) {
  const column = TABLE_COLUMNS.find((item) => item.key === key);
  const multiplier = direction === "asc" ? 1 : -1;
  if (!column || column.type === "string") {
    return multiplier * String(left[key]).localeCompare(String(right[key]));
  }
  return multiplier * ((left[key] || 0) - (right[key] || 0));
}

function defaultDirectionForKey(key) {
  const column = TABLE_COLUMNS.find((item) => item.key === key);
  return column?.type === "number" ? "desc" : "asc";
}

function uniqueSorted(values) {
  return Array.from(new Set(values)).sort((left, right) => left.localeCompare(right));
}

function orderedCategories(values) {
  const seen = new Set(values);
  const ordered = PREFERRED_REGIME_ORDER.filter((value) => seen.has(value));
  const rest = Array.from(seen).filter((value) => !ordered.includes(value)).sort((left, right) => left.localeCompare(right));
  return [...ordered, ...rest];
}

function dataExtent(values) {
  const min = Math.min(...values);
  const max = Math.max(...values);
  return [min, max];
}

function firstEntry(objectMap) {
  return Object.entries(objectMap).sort((left, right) => right[1] - left[1])[0] || null;
}

function toNumber(value) {
  const numeric = Number(value);
  return Number.isFinite(numeric) ? numeric : 0;
}

function formatNumber(value) {
  return Number.isFinite(value) ? value.toFixed(3) : "N/A";
}

function formatSignedNumber(value) {
  if (!Number.isFinite(value)) {
    return "N/A";
  }
  return value >= 0 ? `+${value.toFixed(3)}` : value.toFixed(3);
}

function formatPercent(value) {
  if (!Number.isFinite(value)) {
    return "N/A";
  }
  return `${(value * 100).toFixed(0)}%`;
}

function formatInteger(value) {
  return Number.isFinite(value) ? String(Math.round(value)) : "N/A";
}

function extractTimestamp(path) {
  const match = String(path).match(/(\d{8}_\d{6})_stage10_/);
  return match ? match[1] : null;
}

function paddedExtent(values) {
  const min = Math.min(...values);
  const max = Math.max(...values);
  if (min === max) {
    return [min - 1, max + 1];
  }
  const pad = (max - min) * 0.08;
  return [min - pad, max + pad];
}

function scaleLinear(value, domainMin, domainMax, rangeMin, rangeMax) {
  const ratio = (value - domainMin) / Math.max(domainMax - domainMin, Number.EPSILON);
  return rangeMin + ratio * (rangeMax - rangeMin);
}

function drawLinearGrid(svg, config) {
  const { width, height, margin, xDomain, yDomain, xLabel, yLabel } = config;
  const yTicks = 5;
  const xTicks = 5;

  for (let index = 0; index <= yTicks; index += 1) {
    const value = yDomain[0] + ((yDomain[1] - yDomain[0]) * index) / yTicks;
    const y = scaleLinear(value, yDomain[0], yDomain[1], height - margin.bottom, margin.top);
    const line = svgLine(margin.left, y, width - margin.right, y, "grid-line");
    svg.appendChild(line);
    const text = svgText(16, y + 4, formatNumber(value), "tick-text");
    text.setAttribute("font-size", "12");
    svg.appendChild(text);
  }

  for (let index = 0; index <= xTicks; index += 1) {
    const value = xDomain[0] + ((xDomain[1] - xDomain[0]) * index) / xTicks;
    const x = scaleLinear(value, xDomain[0], xDomain[1], margin.left, width - margin.right);
    const line = svgLine(x, margin.top, x, height - margin.bottom, "grid-line");
    svg.appendChild(line);
    const text = svgText(x - 12, height - margin.bottom + 24, formatNumber(value), "tick-text");
    text.setAttribute("font-size", "12");
    svg.appendChild(text);
  }

  svg.appendChild(svgLine(margin.left, height - margin.bottom, width - margin.right, height - margin.bottom, "axis-line"));
  svg.appendChild(svgLine(margin.left, margin.top, margin.left, height - margin.bottom, "axis-line"));
  svg.appendChild(svgText(width / 2 - 46, height - 10, xLabel, "axis-text"));

  const yText = svgText(18, height / 2 + 24, yLabel, "axis-text");
  yText.setAttribute("transform", `rotate(-90 18 ${height / 2 + 24})`);
  svg.appendChild(yText);
}

function drawBarGrid(svg, width, height, margin) {
  for (let index = 0; index <= 4; index += 1) {
    const value = index / 4;
    const y = scaleLinear(value, 0, 1, height - margin.bottom, margin.top);
    svg.appendChild(svgLine(margin.left, y, width - margin.right, y, "grid-line"));
    const text = svgText(12, y + 4, formatPercent(value), "tick-text");
    text.setAttribute("font-size", "12");
    svg.appendChild(text);
  }
  svg.appendChild(svgLine(margin.left, margin.top, margin.left, height - margin.bottom, "axis-line"));
  svg.appendChild(svgLine(margin.left, height - margin.bottom, width - margin.right, height - margin.bottom, "axis-line"));
}

function createSectorShape(sectorLabel, x, y, size) {
  const type = SECTOR_SHAPES[sectorLabel] || "triangle";
  if (type === "circle") {
    const circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
    circle.setAttribute("cx", x);
    circle.setAttribute("cy", y);
    circle.setAttribute("r", size);
    return circle;
  }
  if (type === "diamond") {
    const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
    path.setAttribute("d", `M ${x} ${y - size} L ${x + size} ${y} L ${x} ${y + size} L ${x - size} ${y} Z`);
    return path;
  }
  if (type === "square") {
    const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
    rect.setAttribute("x", x - size);
    rect.setAttribute("y", y - size);
    rect.setAttribute("width", size * 2);
    rect.setAttribute("height", size * 2);
    rect.setAttribute("rx", "2");
    return rect;
  }
  const triangle = document.createElementNS("http://www.w3.org/2000/svg", "path");
  triangle.setAttribute(
    "d",
    `M ${x} ${y - size} L ${x + size} ${y + size} L ${x - size} ${y + size} Z`,
  );
  return triangle;
}

function renderShapeGlyph(sectorLabel) {
  const type = SECTOR_SHAPES[sectorLabel] || "triangle";
  if (type === "circle") {
    return "●";
  }
  if (type === "diamond") {
    return "◆";
  }
  if (type === "square") {
    return "■";
  }
  return "▲";
}

function regimeColor(regime) {
  return REGIME_COLORS[regime] || "#8291a2";
}

function heatColor(value, maxValue) {
  const ratio = value / Math.max(maxValue, 1);
  const low = { r: 35, g: 43, b: 54 };
  const high = { r: 132, g: 170, b: 177 };
  const r = Math.round(low.r + ratio * (high.r - low.r));
  const g = Math.round(low.g + ratio * (high.g - low.g));
  const b = Math.round(low.b + ratio * (high.b - low.b));
  return `rgb(${r}, ${g}, ${b})`;
}

function coherenceColor(value) {
  const ratio = clamp(value, 0, 1);
  const low = { r: 179, g: 110, b: 96 };
  const high = { r: 126, g: 162, b: 168 };
  const r = Math.round(low.r + ratio * (high.r - low.r));
  const g = Math.round(low.g + ratio * (high.g - low.g));
  const b = Math.round(low.b + ratio * (high.b - low.b));
  return `rgb(${r}, ${g}, ${b})`;
}

function waveParticleColor(proxy) {
  const wave = { r: 190, g: 142, b: 88 };
  const particle = { r: 118, g: 168, b: 174 };
  const total = Math.max(proxy.wave + proxy.particle, Number.EPSILON);
  const particleMix = proxy.particle / total;
  const r = Math.round(wave.r + particleMix * (particle.r - wave.r));
  const g = Math.round(wave.g + particleMix * (particle.g - wave.g));
  const b = Math.round(wave.b + particleMix * (particle.b - wave.b));
  return `rgb(${r}, ${g}, ${b})`;
}

function waveParticleProxy(run) {
  const extents = state.atlas0.metricExtents;
  if (!extents) {
    return { wave: 0.5, particle: 0.5, label: "mixed proxy" };
  }

  const widthScore = normalizedMetric(run.width_ratio, extents.width_ratio);
  const spreadScore = normalizedMetric(run.spectral_spread_final, extents.spectral_spread_final);
  const coherenceScore = normalizedMetric(run.coherence_score, extents.coherence_score);
  const protoScore = normalizedMetric(run.proto_spacetime_score, extents.proto_spacetime_score);
  const wave = clamp(widthScore * 0.45 + spreadScore * 0.35 + (1 - protoScore) * 0.2, 0, 1);
  const particle = clamp(
    coherenceScore * 0.4 + protoScore * 0.35 + (1 - widthScore) * 0.15 + (1 - spreadScore) * 0.1,
    0,
    1,
  );

  let label = "mixed proxy";
  if (wave > particle + 0.08) {
    label = "wave-leaning proxy";
  } else if (particle > wave + 0.08) {
    label = "particle-leaning proxy";
  }

  return { wave, particle, label };
}

function propagationClassForRun(run) {
  const mapping = {
    "Ballistic coherent": "compact propagation",
    "Ballistic dispersive": "structured spread",
    Diffusive: "diffuse spread",
    Fragmenting: "fragmented spread",
    "Chaotic or irregular": "irregular spread",
    Localized: "compact trapping",
    "Oscillatory trapped": "bounded oscillation",
    "Metastable structured": "metastable propagation",
  };
  return mapping[run.regime_label] || "mixed propagation";
}

function prepareCanvas(canvas) {
  const dpr = window.devicePixelRatio || 1;
  const width = Math.max(canvas.clientWidth || canvas.width || 640, 320);
  const height = Math.max(canvas.clientHeight || canvas.height || 360, 240);
  const nextWidth = Math.round(width * dpr);
  const nextHeight = Math.round(height * dpr);

  if (canvas.width !== nextWidth || canvas.height !== nextHeight) {
    canvas.width = nextWidth;
    canvas.height = nextHeight;
  }

  const ctx = canvas.getContext("2d");
  if (!ctx) {
    return { ctx: null, width, height };
  }

  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, width, height);
  return { ctx, width, height };
}

function wrappedDisplacement(point, center, boundaryType) {
  const delta = [
    (point[0] || 0) - (center[0] || 0),
    (point[1] || 0) - (center[1] || 0),
    (point[2] || 0) - (center[2] || 0),
  ];

  if (boundaryType !== "periodic") {
    return delta;
  }

  return delta.map((value) => ((((value + 0.5) % 1) + 1) % 1) - 0.5);
}

function wrappedVelocityVector(previousCenter, currentCenter, boundaryType) {
  return wrappedDisplacement(currentCenter, previousCenter, boundaryType);
}

function flowSigmaFromWidth(width) {
  return clamp(width * 0.46, 0.02, 0.18);
}

function normalizedMetric(value, extent) {
  const [min, max] = extent || [0, 1];
  if (min === max) {
    return 0.5;
  }
  return clamp((value - min) / (max - min), 0, 1);
}

function normalizeToCube(value, extent) {
  const [min, max] = extent;
  if (min === max) {
    return 0;
  }
  return ((value - min) / (max - min) - 0.5) * 1.2;
}

function rotate3d(point, rotationX, rotationY) {
  const cosY = Math.cos(rotationY);
  const sinY = Math.sin(rotationY);
  const x1 = point.x * cosY + point.z * sinY;
  const z1 = -point.x * sinY + point.z * cosY;

  const cosX = Math.cos(rotationX);
  const sinX = Math.sin(rotationX);
  const y2 = point.y * cosX - z1 * sinX;
  const z2 = point.y * sinX + z1 * cosX;

  return { x: x1, y: y2, z: z2 };
}

function projectScenePoint(point, viewState, width, height, scaleFactor = 1) {
  const rotated = rotate3d(point, viewState.rotationX, viewState.rotationY);
  const cameraDistance = 2.65;
  const perspective = cameraDistance / Math.max(cameraDistance - rotated.z, 0.6);
  const scale = Math.min(width, height) * 0.38 * viewState.zoom * scaleFactor;
  return {
    x: width * 0.5 + rotated.x * scale * perspective,
    y: height * 0.5 - rotated.y * scale * perspective,
    depth: rotated.z,
    perspective,
  };
}

function colorWithAlpha(color, alpha) {
  const hex = String(color).replace("#", "");
  if (hex.length !== 6) {
    return color;
  }

  const red = Number.parseInt(hex.slice(0, 2), 16);
  const green = Number.parseInt(hex.slice(2, 4), 16);
  const blue = Number.parseInt(hex.slice(4, 6), 16);
  return `rgba(${red}, ${green}, ${blue}, ${clamp(alpha, 0, 1)})`;
}

function draw3dScaffold(ctx, width, height, viewState, labels) {
  const corners = [
    { x: -0.65, y: -0.65, z: -0.65 },
    { x: 0.65, y: -0.65, z: -0.65 },
    { x: -0.65, y: 0.65, z: -0.65 },
    { x: 0.65, y: 0.65, z: -0.65 },
    { x: -0.65, y: -0.65, z: 0.65 },
    { x: 0.65, y: -0.65, z: 0.65 },
    { x: -0.65, y: 0.65, z: 0.65 },
    { x: 0.65, y: 0.65, z: 0.65 },
  ].map((corner) => projectScenePoint(corner, viewState, width, height));

  const edges = [
    [0, 1], [0, 2], [1, 3], [2, 3],
    [4, 5], [4, 6], [5, 7], [6, 7],
    [0, 4], [1, 5], [2, 6], [3, 7],
  ];

  ctx.save();
  ctx.lineWidth = 1;
  ctx.strokeStyle = "rgba(157, 177, 196, 0.18)";
  edges.forEach(([start, end]) => {
    ctx.beginPath();
    ctx.moveTo(corners[start].x, corners[start].y);
    ctx.lineTo(corners[end].x, corners[end].y);
    ctx.stroke();
  });

  ctx.fillStyle = "#afbbc8";
  ctx.font = "12px IBM Plex Sans, sans-serif";
  const xLabel = projectScenePoint({ x: 0.82, y: -0.72, z: -0.72 }, viewState, width, height);
  const yLabel = projectScenePoint({ x: -0.76, y: 0.82, z: -0.72 }, viewState, width, height);
  const zLabel = projectScenePoint({ x: -0.74, y: -0.74, z: 0.82 }, viewState, width, height);
  ctx.fillText(labels.x, xLabel.x, xLabel.y);
  ctx.fillText(labels.y, yLabel.x, yLabel.y);
  ctx.fillText(labels.z, zLabel.x, zLabel.y);
  ctx.restore();
}

function drawCanvasShape(ctx, sectorLabel, x, y, radius, fill, isSelected, isHovered) {
  ctx.save();
  ctx.fillStyle = fill;
  ctx.globalAlpha = isHovered ? 1 : 0.92;
  ctx.beginPath();

  const type = SECTOR_SHAPES[sectorLabel] || "triangle";
  if (type === "circle") {
    ctx.arc(x, y, radius, 0, Math.PI * 2);
  } else if (type === "diamond") {
    ctx.moveTo(x, y - radius);
    ctx.lineTo(x + radius, y);
    ctx.lineTo(x, y + radius);
    ctx.lineTo(x - radius, y);
    ctx.closePath();
  } else if (type === "square") {
    ctx.rect(x - radius, y - radius, radius * 2, radius * 2);
  } else {
    ctx.moveTo(x, y - radius);
    ctx.lineTo(x + radius, y + radius);
    ctx.lineTo(x - radius, y + radius);
    ctx.closePath();
  }
  ctx.fill();

  if (isSelected || isHovered) {
    ctx.strokeStyle = isSelected ? "#f4f8fb" : "rgba(244, 248, 251, 0.72)";
    ctx.lineWidth = isSelected ? 2.2 : 1.4;
    ctx.stroke();
  }
  ctx.restore();
}

function pointRadius(item, baseRadius) {
  return baseRadius + item.row.proto_spacetime_score * 0.7 + item.perspective * 0.6;
}

function canvasPointerPosition(canvas, event) {
  const rect = canvas.getBoundingClientRect();
  return {
    x: event.clientX - rect.left,
    y: event.clientY - rect.top,
  };
}

function findNearestProjectedPoint(points, x, y, threshold) {
  let nearest = null;
  let bestDistance = threshold;

  points.forEach((item) => {
    const distance = Math.hypot(item.x - x, item.y - y);
    if (distance <= bestDistance) {
      bestDistance = distance;
      nearest = item;
    }
  });

  return nearest;
}

function clamp(value, min, max) {
  return Math.min(max, Math.max(min, value));
}

function clearSvg(svg) {
  while (svg.firstChild) {
    svg.removeChild(svg.firstChild);
  }
}

function svgLine(x1, y1, x2, y2, className) {
  const line = document.createElementNS("http://www.w3.org/2000/svg", "line");
  line.setAttribute("x1", x1);
  line.setAttribute("y1", y1);
  line.setAttribute("x2", x2);
  line.setAttribute("y2", y2);
  line.setAttribute("class", className);
  return line;
}

function svgText(x, y, textContent, className) {
  const text = document.createElementNS("http://www.w3.org/2000/svg", "text");
  text.setAttribute("x", x);
  text.setAttribute("y", y);
  text.setAttribute("class", className);
  text.textContent = textContent;
  return text;
}

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}

function handleFatalError(error) {
  console.error(error);
  elements.atlas0Status.textContent = "Load failed";
  elements.runTableBody.innerHTML = `
    <tr>
      <td colspan="${TABLE_COLUMNS.length}" class="faint">
        Atlas Explorer could not load its static inputs. Start a local HTTP server from the HAOS-IIP directory and verify the paths at the top of app.js.
      </td>
    </tr>
  `;
  elements.selectionBrief.textContent = error instanceof Error ? error.message : "Unknown load failure.";
}
