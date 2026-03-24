const THRESHOLDS = {
  maxWidthGrowth: 0.10,
  minConcentration: 0.88,
  maxParticipationGrowth: 0.20,
  minOverlap: 0.90,
  minRecoveryScore: 0.90,
};

const CLASS_COLORS = {
  stable: "#0f766e",
  marginal: "#b45309",
  unstable: "#b91c1c",
};

const scenarioOrder = ["baseline", "perturbed", "fragmented"];

const architectureNodes = [
  {
    id: "state",
    label: "State layer",
    status: "Frozen",
    summary:
      "Seed trajectories live as deterministic JSON scenarios. This is the public entry point for small stability experiments without touching frozen phase bundles.",
    points: [
      "Scenarios stay human-readable and machine-readable at the same time.",
      "The playground starts from these seeds, then applies deterministic front-end modifiers.",
      "Because the base scenarios are static assets, the site can run entirely over a simple local server.",
    ],
    refs: [
      { label: "examples/scenarios/baseline.json", href: "../examples/scenarios/baseline.json" },
      { label: "examples/scenarios/perturbed.json", href: "../examples/scenarios/perturbed.json" },
      { label: "examples/scenarios/fragmented.json", href: "../examples/scenarios/fragmented.json" },
    ],
  },
  {
    id: "metrics",
    label: "Metrics layer",
    status: "Frozen",
    summary:
      "The core readout comes from the frozen telemetry layer: persistence, ordering, causal-depth drift, and distance coherence.",
    points: [
      "These are the metrics that let the oracle ask whether coherence survives interaction.",
      "The public aliases on the site mirror the CLI naming: structural retention, temporal consistency, causal deformation, and geometric integrity.",
      "The threshold structure is fixed in the public demo, which keeps comparisons stable.",
    ],
    refs: [
      { label: "telemetry/frozen_metrics.py", href: "../telemetry/frozen_metrics.py" },
      { label: "docs note: structural stability oracle", href: "../docs/notes/applications/A_Minimal_Structural_Stability_Oracle_Based_on_Frozen_HAOS_IIP_Telemetry_v1.md" },
    ],
  },
  {
    id: "policy",
    label: "Policy layer",
    status: "Public",
    summary:
      "Policy is the small public decision rule that compresses frozen telemetry into stable, marginal, or unstable.",
    points: [
      "It stays deliberately narrow: classify structural survival, not ontology.",
      "The ladder is explicit and auditable in code, which makes the oracle inspectable instead of mystical.",
      "This layer is also where the site derives coherence score and confidence for a cold reader.",
    ],
    refs: [
      { label: "haos_iip/demo.py", href: "../haos_iip/demo.py" },
      { label: "paper 45.1", href: "../papers/pdf_releases/45.1%20A%20Minimal%20Structural-Stability%20Oracle%20Based%20on%20Frozen%20HAOS-IIP%20Telemetry.pdf" },
    ],
  },
  {
    id: "router",
    label: "Router layer",
    status: "Planned",
    summary:
      "The next public step is a thin routing surface that can dispatch scenarios, policies, and domain wrappers through the same telemetry contract.",
    points: [
      "This is the natural bridge between the current CLI tools and a broader web or API surface.",
      "The site treats routing as planned work, not as something already solved in the repo.",
      "A FastAPI wrapper or WASM port can replace the current client-side evaluator later without changing the page structure.",
    ],
    refs: [
      { label: "API_CONTRACT.md", href: "../API_CONTRACT.md" },
      { label: "run_phase.py", href: "../run_phase.py" },
      { label: "examples/quick_reproduce.py", href: "../examples/quick_reproduce.py" },
    ],
  },
  {
    id: "future",
    label: "Future domains",
    status: "Planned",
    summary:
      "Future domains are where the oracle could move from a single demo surface into agents, domain-specific probes, and practical stability checks.",
    points: [
      "The repo already contains agent-oriented scaffolding and workflow notes that can connect to this surface.",
      "A router plus domain adapters would let the same public metrics assess richer systems.",
      "The right next move is controlled expansion, not a jump to broad claims.",
    ],
    refs: [
      { label: "ai/agents", href: "../ai/agents/" },
      { label: "PROJECT_STATUS.md", href: "../PROJECT_STATUS.md" },
      { label: "RESEARCH_MAP.md", href: "../RESEARCH_MAP.md" },
    ],
  },
];

const papers = [
  {
    id: "oracle",
    code: "45.1",
    tag: "Tool surface",
    title: "A Minimal Structural-Stability Oracle Based on Frozen HAOS-IIP Telemetry",
    href: "../papers/pdf_releases/45.1%20A%20Minimal%20Structural-Stability%20Oracle%20Based%20on%20Frozen%20HAOS-IIP%20Telemetry.pdf",
    conceptual:
      "This is the shortest bridge from the research stack to a public instrument. It says: take the frozen telemetry layer and make it callable.",
    technical:
      "The paper compresses the public demo around persistence, ordering, causal-depth drift, and distance coherence, then exposes the stable/marginal/unstable ladder as a small oracle.",
  },
  {
    id: "master",
    code: "46.1",
    tag: "Synthesis",
    title: "HAOS-IIP Master Synthesis: Program-Level Compression of the Paper Spine",
    href: "../papers/pdf_releases/46.1%20HAOS-IIP%20Master%20Synthesis%20Program-Level%20Compression%20of%20the%20Paper%20Spine.pdf",
    conceptual:
      "If someone wants the whole arc without reading every release first, this is the map that compresses architecture, limits, and public tool surface into one document.",
    technical:
      "It preserves frozen claims while summarizing the authority chain from operator architecture through temporal ordering, causal closure, distance surrogate feasibility, and the public oracle.",
  },
  {
    id: "readout",
    code: "28.1",
    tag: "Readout authority",
    title: "Deterministic Coherence-Recovery Readout, Robustness, and Control Discrimination on a Frozen Phase IV Sector",
    href: "../papers/pdf_releases/28.1%20Deterministic%20Coherence-Recovery%20Readout,%20Robustness,%20and%20Control%20Discrimination%20on%20a%20Frozen%20Phase%20IV%20Sector.pdf",
    conceptual:
      "This is where the project sharpens from visual pattern hunting into a readout discipline. Recovery has to be discriminative, not just interesting.",
    technical:
      "The release establishes deterministic coherence-recovery readout, robustness checks, and control discrimination on the frozen Phase IV sector, which is the direct ancestor of the public oracle surface.",
  },
  {
    id: "continuum",
    code: "33.1",
    tag: "Bridge discipline",
    title: "Cautious Continuum-Bridge Feasibility of a Frozen Branch-Local Cochain-Laplacian Hierarchy",
    href: "../papers/pdf_releases/33.1%20Cautious%20Continuum-Bridge%20Feasibility%20of%20a%20Frozen%20Branch-Local%20Cochain-Laplacian%20Hierarchy.pdf",
    conceptual:
      "This paper matters because it makes restraint visible. The bridge is cautious by design, which is exactly why it adds credibility.",
    technical:
      "It frames scaling inspection as feasibility only and keeps the program bounded to reproducible post-processing rather than ontological escalation.",
  },
  {
    id: "ordering",
    code: "43.1",
    tag: "Late-stage feasibility",
    title: "Ordering, Causal Closure, and Proto-Geometric Distance-Surrogate Feasibility",
    href: "../papers/pdf_releases/43.1%20Ordering,%20Causal%20Closure,%20and%20Proto-Geometric%20Distance-Surrogate%20Feasibility%20on%20a%20Frozen%20Branch-Local%20Cochain-Laplacian%20Hierarchy.pdf",
    conceptual:
      "This is the point where the stack starts to feel system-like rather than phase-by-phase, because temporal order, causal closure, and distance surrogate live in one frame.",
    technical:
      "It integrates the frozen late-stage feasibility results that the public oracle later compresses into a more accessible diagnostic surface.",
  },
  {
    id: "atlas",
    code: "09.1",
    tag: "Atlas instrument",
    title: "A Pre-Geometric Atlas of Morphology, Resilience, and Transition Structure",
    href: "../papers/pdf_releases/09.1%20A%20Pre-Geometric%20Atlas%20of%20Morphology,%20Resilience,%20and%20Transition%20Structure%20on%20a%20Frozen%20Kernel-Induced%20Cochain%20Operator%20Architecture.pdf",
    conceptual:
      "This is the paper for seeing how the project learned to classify before interpreting. It gives the repository a visual survey language.",
    technical:
      "Atlas-0, Atlas-1, and Atlas-2 establish a reproducible morphology, resilience, and transition instrument that still informs the newer public explorer surface.",
  },
];

const atlasArtifacts = [
  {
    id: "overview",
    label: "Oracle overview",
    title: "Public stability ladder overview",
    src: "../stability_demo_overview.svg",
    summary:
      "The exact default CLI bundle rendered as a compact overview: baseline -> stable, perturbed -> marginal, fragmented -> unstable.",
  },
  {
    id: "scan",
    label: "Scan heatmap",
    title: "Baseline scan across noise and connectivity drop",
    src: "../examples/output/stability_baseline_scan_heatmap.svg",
    summary:
      "The public scan surface is tiny on purpose. It shows how quickly recoverability collapses once connectivity distortion enters the picture.",
  },
  {
    id: "robustness",
    label: "Robustness map",
    title: "Stage 23.8 coherence robustness heatmap",
    src: "../plots/20260315_185337_stage23_8_coherence_robustness_heatmap.png",
    summary:
      "A real repository artifact that gives the abstract language of coherence a concrete visual map.",
  },
  {
    id: "collapse",
    label: "Collapse panel",
    title: "Stage 23.13 coherence collapse panel",
    src: "../plots/20260316_110226_stage23_13_coherence_collapse_panel.png",
    summary:
      "This is where the viewer can see collapse as a diagnostic pattern rather than a dramatic metaphor.",
  },
  {
    id: "prediction",
    label: "Prediction check",
    title: "Stage 24.2 prediction versus truth panel",
    src: "../plots/20260316_115346_stage24_2_prediction_vs_truth_panel.png",
    summary:
      "A useful credibility panel: predicted structure is compared directly against observed output instead of being left as prose.",
  },
  {
    id: "mechanism",
    label: "Mechanism matrix",
    title: "Stage 23.13 mechanism branch comparison matrix",
    src: "../plots/20260316_110320_stage23_13_mechanism_branch_comparison_matrix.png",
    summary:
      "A fast way to show branching behavior, comparison structure, and why the system is starting to look like an engineering instrument.",
  },
];

const roadmapSteps = [
  {
    id: "phase-iii",
    title: "Phase III opener",
    status: "Shipped",
    goal: "Open the DK collision and braid topology line while staying in a controlled, classification-first regime.",
    result:
      "The project gained a richer texture vocabulary and real protection structure without pretending the whole system was solved.",
    artifacts: [
      { label: "paper 23.1", href: "../papers/pdf_releases/23.1%20Dirac-Kahler%20Collision%20Texture%20and%20Grade-Exchange%20Baseline%20on%20a%20Kernel-Induced%20Cochain%20Operator%20Architecture.pdf" },
      { label: "paper 23.2", href: "../papers/pdf_releases/23.2%20Dirac-Kahler%20Braid%20Topology%20and%20Protection%20Structure%20on%20a%20Kernel-Induced%20Cochain%20Operator%20Architecture.pdf" },
      { label: "robustness heatmap artifact", href: "../plots/20260315_185337_stage23_8_coherence_robustness_heatmap.png" },
    ],
    questions: [
      "Which texture classes stay robust under stronger controls?",
      "Where does protection stop being descriptive and start becoming reusable readout?",
    ],
  },
  {
    id: "phase-iv",
    title: "Phase IV sector",
    status: "Shipped",
    goal: "Freeze a discriminative sector where coherence-recovery can be read out deterministically.",
    result:
      "This phase created the practical bridge from experimental morphology toward a callable readout discipline.",
    artifacts: [
      { label: "phase IV paper 27.1", href: "../papers/pdf_releases/27.1%20phase_iv_paper_numbered.pdf" },
      { label: "Phase V diagnostic image", href: "../Images/HAOS-IIP%20Phase%20V%20Diagnostic%20Data%20Readout.png" },
    ],
    questions: [
      "How thin can the readout surface become before it loses discrimination?",
      "Which controls most cleanly separate survival from collapse?",
    ],
  },
  {
    id: "phase-v",
    title: "Phase V readout authority",
    status: "Shipped",
    goal: "Turn coherence-recovery into an auditable deterministic readout with robustness and control discrimination.",
    result:
      "The public explorer is directly downstream of this move because it depends on a stable vocabulary for evaluating recovery.",
    artifacts: [
      { label: "paper 28.1", href: "../papers/pdf_releases/28.1%20Deterministic%20Coherence-Recovery%20Readout,%20Robustness,%20and%20Control%20Discrimination%20on%20a%20Frozen%20Phase%20IV%20Sector.pdf" },
      { label: "API contract", href: "../API_CONTRACT.md" },
    ],
    questions: [
      "How do we expose readout without overselling the ontology?",
      "What is the thinnest reliable public interface?",
    ],
  },
  {
    id: "phase-vi",
    title: "Phase VI operator freeze",
    status: "Shipped",
    goal: "Freeze the branch-local cochain-Laplacian family so later tools can reference a stable operator layer.",
    result:
      "This gave the program a stronger reproducibility spine and made later public tooling easier to anchor cleanly.",
    artifacts: [
      { label: "paper 29.1", href: "../papers/pdf_releases/29.1%20Frozen%20Branch-Local%20Cochain-Laplacian%20Operator%20Family%20and%20Spectral%20Feasibility%20Across%20a%20Deterministic%20Refinement%20Hierarchy.pdf" },
      { label: "README quick start", href: "../README.md" },
    ],
    questions: [
      "Which later public tools should reference the operator layer directly versus only its telemetry products?",
      "How much operator detail belongs in the public narrative?",
    ],
  },
  {
    id: "oracle",
    title: "Stability oracle",
    status: "Shipped",
    goal: "Expose a public diagnostic layer that compresses frozen telemetry into a callable stability ladder.",
    result:
      "The repository now has an actual tool surface that can be invoked from CLI, summarized in JSON, and turned into an interactive web experience.",
    artifacts: [
      { label: "paper 45.1", href: "../papers/pdf_releases/45.1%20A%20Minimal%20Structural-Stability%20Oracle%20Based%20on%20Frozen%20HAOS-IIP%20Telemetry.pdf" },
      { label: "haos_iip/demo.py", href: "../haos_iip/demo.py" },
      { label: "stability demo overview", href: "../stability_demo_overview.svg" },
    ],
    questions: [
      "Which domain wrappers can reuse the ladder without bending the meaning of the metrics?",
      "Where should confidence and aggregate coherence be computed: client, API, or both?",
    ],
  },
  {
    id: "router",
    title: "Router surface",
    status: "Planned",
    goal: "Route scenarios, policies, and future domain adapters through one thin public interface.",
    result:
      "Not shipped yet. The explorer presents the router as the natural next implementation layer after the current static and CLI surfaces.",
    artifacts: [
      { label: "API contract", href: "../API_CONTRACT.md" },
      { label: "examples/quick_reproduce.py", href: "../examples/quick_reproduce.py" },
      { label: "project status", href: "../PROJECT_STATUS.md" },
    ],
    questions: [
      "Should the first router be HTTP-first, CLI-first, or both?",
      "Which parts of the scenario grammar should remain frozen once external users depend on them?",
    ],
  },
  {
    id: "agents",
    title: "Agent integration",
    status: "Planned",
    goal: "Let future agents call the oracle, inspect artifacts, and route stability checks without bypassing the frozen telemetry vocabulary.",
    result:
      "This remains a forward path. The repo already contains agent-oriented directories, but the site marks integration as upcoming rather than implied.",
    artifacts: [
      { label: "ai/agents", href: "../ai/agents/" },
      { label: "haos_iip/demo.py", href: "../haos_iip/demo.py" },
    ],
    questions: [
      "How should an agent report uncertainty while staying grounded in recoverable coherence?",
      "What safeguards keep the system instrument-like instead of becoming a belief engine?",
    ],
  },
];

const state = {
  scenarios: {},
  presetResults: {},
  activeScenario: "baseline",
  paperMode: "conceptual",
  activeArchitecture: "state",
  activeArtifact: "overview",
  activeRoadmap: "phase-iii",
  modifiers: {
    noise: 0,
    connectivityDrop: 0,
    perturbation: 0,
    temporalDistortion: 0,
  },
  evaluation: null,
  animationTick: 0,
  animationFrameId: null,
};

const elements = {};

document.addEventListener("DOMContentLoaded", async () => {
  cacheElements();
  bindEvents();
  await loadScenarios();
  syncModifierInputs();
  renderPresetCards();
  renderArchitecture();
  renderPapers();
  renderAtlas();
  renderRoadmap();
  renderDeveloperPayload();
  evaluateCurrentScenario();
  startSignalAnimation();
});

function cacheElements() {
  elements.presetGrid = document.getElementById("preset-grid");
  elements.noiseSlider = document.getElementById("noise-slider");
  elements.connectivitySlider = document.getElementById("connectivity-slider");
  elements.perturbationSlider = document.getElementById("perturbation-slider");
  elements.temporalSlider = document.getElementById("temporal-slider");
  elements.noiseValue = document.getElementById("noise-value");
  elements.connectivityValue = document.getElementById("connectivity-value");
  elements.perturbationValue = document.getElementById("perturbation-value");
  elements.temporalValue = document.getElementById("temporal-value");
  elements.evaluateButton = document.getElementById("evaluate-button");
  elements.resetButton = document.getElementById("reset-button");
  elements.resultTitle = document.getElementById("result-title");
  elements.resultSubtitle = document.getElementById("result-subtitle");
  elements.resultBadge = document.getElementById("result-badge");
  elements.resultStats = document.getElementById("result-stats");
  elements.metricBars = document.getElementById("metric-bars");
  elements.heatmapGrid = document.getElementById("heatmap-grid");
  elements.heatmapCaption = document.getElementById("heatmap-caption");
  elements.oracleJson = document.getElementById("oracle-json");
  elements.signalCanvas = document.getElementById("signal-canvas");
  elements.architectureLanes = document.getElementById("architecture-lanes");
  elements.architectureTitle = document.getElementById("architecture-title");
  elements.architectureStatus = document.getElementById("architecture-status");
  elements.architectureSummary = document.getElementById("architecture-summary");
  elements.architecturePoints = document.getElementById("architecture-points");
  elements.architectureRefs = document.getElementById("architecture-refs");
  elements.conceptualMode = document.getElementById("conceptual-mode");
  elements.technicalMode = document.getElementById("technical-mode");
  elements.papersGrid = document.getElementById("papers-grid");
  elements.atlasImage = document.getElementById("atlas-image");
  elements.atlasLabel = document.getElementById("atlas-label");
  elements.atlasTitle = document.getElementById("atlas-title");
  elements.atlasSummary = document.getElementById("atlas-summary");
  elements.atlasLink = document.getElementById("atlas-link");
  elements.atlasThumbs = document.getElementById("atlas-thumbs");
  elements.roadmapTimeline = document.getElementById("roadmap-timeline");
  elements.roadmapTitle = document.getElementById("roadmap-title");
  elements.roadmapStatus = document.getElementById("roadmap-status");
  elements.roadmapGoal = document.getElementById("roadmap-goal");
  elements.roadmapResult = document.getElementById("roadmap-result");
  elements.roadmapArtifacts = document.getElementById("roadmap-artifacts");
  elements.roadmapQuestions = document.getElementById("roadmap-questions");
  elements.developerJson = document.getElementById("developer-json");
}

function bindEvents() {
  const sliderBindings = [
    [elements.noiseSlider, "noise", elements.noiseValue],
    [elements.connectivitySlider, "connectivityDrop", elements.connectivityValue],
    [elements.perturbationSlider, "perturbation", elements.perturbationValue],
    [elements.temporalSlider, "temporalDistortion", elements.temporalValue],
  ];

  for (const [input, key, output] of sliderBindings) {
    input.addEventListener("input", () => {
      state.modifiers[key] = toNumber(input.value);
      output.textContent = formatFixed(state.modifiers[key], key === "perturbation" || key === "temporalDistortion" ? 2 : 2);
    });
  }

  elements.evaluateButton.addEventListener("click", () => {
    evaluateCurrentScenario();
  });

  elements.resetButton.addEventListener("click", () => {
    state.modifiers = {
      noise: 0,
      connectivityDrop: 0,
      perturbation: 0,
      temporalDistortion: 0,
    };
    syncModifierInputs();
    evaluateCurrentScenario();
  });

  elements.conceptualMode.addEventListener("click", () => {
    state.paperMode = "conceptual";
    renderPapers();
  });

  elements.technicalMode.addEventListener("click", () => {
    state.paperMode = "technical";
    renderPapers();
  });

  window.addEventListener("resize", () => {
    drawSignalCanvas();
  });
}

async function loadScenarios() {
  const loaded = await Promise.all(
    scenarioOrder.map(async (name) => {
      const response = await fetch(`../examples/scenarios/${name}.json`);
      if (!response.ok) {
        throw new Error(`Unable to load scenario ${name}`);
      }
      return [name, await response.json()];
    }),
  );

  state.scenarios = Object.fromEntries(loaded);

  for (const name of scenarioOrder) {
    state.presetResults[name] = evaluateScenario(state.scenarios[name]);
  }
}

function renderPresetCards() {
  const fragment = document.createDocumentFragment();
  for (const name of scenarioOrder) {
    const result = state.presetResults[name];
    const button = document.createElement("button");
    button.type = "button";
    button.className = `preset-card ${state.activeScenario === name ? "active" : ""}`;
    button.addEventListener("click", () => {
      state.activeScenario = name;
      state.modifiers = {
        noise: 0,
        connectivityDrop: 0,
        perturbation: 0,
        temporalDistortion: 0,
      };
      syncModifierInputs();
      renderPresetCards();
      evaluateCurrentScenario();
    });

    button.innerHTML = `
      <div class="preset-head">
        <div class="preset-title">${escapeHtml(name)}</div>
        <span class="oracle-badge ${result.classification}">${escapeHtml(result.classification)}</span>
      </div>
      <p class="preset-summary">${escapeHtml(result.description)}</p>
      <div class="preset-metrics">
        <span>retention ${formatMetric(result.structural_retention)}</span>
        <span>consistency ${formatMetric(result.temporal_consistency)}</span>
      </div>
    `;
    fragment.appendChild(button);
  }
  elements.presetGrid.replaceChildren(fragment);
}

function renderArchitecture() {
  const fragment = document.createDocumentFragment();
  for (const node of architectureNodes) {
    const button = document.createElement("button");
    button.type = "button";
    button.className = `architecture-node ${state.activeArchitecture === node.id ? "active" : ""}`;
    button.innerHTML = `<strong>${escapeHtml(node.label)}</strong><span>${escapeHtml(node.summary)}</span>`;
    button.addEventListener("click", () => {
      state.activeArchitecture = node.id;
      renderArchitecture();
      renderArchitectureDetail();
    });
    fragment.appendChild(button);
  }
  elements.architectureLanes.replaceChildren(fragment);
  renderArchitectureDetail();
}

function renderArchitectureDetail() {
  const node = architectureNodes.find((item) => item.id === state.activeArchitecture) || architectureNodes[0];
  elements.architectureTitle.textContent = node.label;
  elements.architectureStatus.textContent = node.status;
  elements.architectureSummary.textContent = node.summary;
  elements.architecturePoints.replaceChildren(...node.points.map((point) => createListItem(point)));
  elements.architectureRefs.replaceChildren(...node.refs.map(createReferenceLink));
}

function renderPapers() {
  elements.conceptualMode.classList.toggle("active", state.paperMode === "conceptual");
  elements.technicalMode.classList.toggle("active", state.paperMode === "technical");

  const fragment = document.createDocumentFragment();
  for (const paper of papers) {
    const article = document.createElement("article");
    article.className = "paper-card";
    article.innerHTML = `
      <div class="paper-meta">
        <span class="paper-tag">${escapeHtml(paper.tag)}</span>
        <span class="paper-code">${escapeHtml(paper.code)}</span>
      </div>
      <h3>${escapeHtml(paper.title)}</h3>
      <p class="paper-summary">${escapeHtml(paper[state.paperMode])}</p>
      <div class="paper-actions">
        <a class="button button-secondary" href="${paper.href}" target="_blank" rel="noreferrer">Open PDF</a>
      </div>
    `;
    fragment.appendChild(article);
  }
  elements.papersGrid.replaceChildren(fragment);
}

function renderAtlas() {
  const active = atlasArtifacts.find((item) => item.id === state.activeArtifact) || atlasArtifacts[0];
  elements.atlasImage.src = active.src;
  elements.atlasImage.alt = active.title;
  elements.atlasLabel.textContent = active.label;
  elements.atlasTitle.textContent = active.title;
  elements.atlasSummary.textContent = active.summary;
  elements.atlasLink.href = active.src;

  const fragment = document.createDocumentFragment();
  for (const artifact of atlasArtifacts) {
    const button = document.createElement("button");
    button.type = "button";
    button.className = `atlas-thumb ${state.activeArtifact === artifact.id ? "active" : ""}`;
    button.innerHTML = `<strong>${escapeHtml(artifact.label)}</strong><span>${escapeHtml(artifact.title)}</span>`;
    button.addEventListener("click", () => {
      state.activeArtifact = artifact.id;
      renderAtlas();
    });
    fragment.appendChild(button);
  }
  elements.atlasThumbs.replaceChildren(fragment);
}

function renderRoadmap() {
  const fragment = document.createDocumentFragment();
  for (const step of roadmapSteps) {
    const button = document.createElement("button");
    button.type = "button";
    button.className = `roadmap-step ${state.activeRoadmap === step.id ? "active" : ""}`;
    button.innerHTML = `
      <div class="roadmap-step-head">
        <strong>${escapeHtml(step.title)}</strong>
        <span>${escapeHtml(step.status)}</span>
      </div>
      <span>${escapeHtml(step.goal)}</span>
    `;
    button.addEventListener("click", () => {
      state.activeRoadmap = step.id;
      renderRoadmap();
      renderRoadmapDetail();
    });
    fragment.appendChild(button);
  }
  elements.roadmapTimeline.replaceChildren(fragment);
  renderRoadmapDetail();
}

function renderRoadmapDetail() {
  const step = roadmapSteps.find((item) => item.id === state.activeRoadmap) || roadmapSteps[0];
  elements.roadmapTitle.textContent = step.title;
  elements.roadmapStatus.textContent = step.status;
  elements.roadmapGoal.textContent = step.goal;
  elements.roadmapResult.textContent = step.result;
  elements.roadmapArtifacts.replaceChildren(...step.artifacts.map(createReferenceLink));
  elements.roadmapQuestions.replaceChildren(...step.questions.map((question) => createListItem(question)));
}

function renderDeveloperPayload() {
  const baseline = state.presetResults.baseline;
  const payload = {
    scenario: baseline.scenario,
    classification: baseline.classification,
    structural_retention: roundMetric(baseline.structural_retention),
    temporal_consistency: roundMetric(baseline.temporal_consistency),
    causal_deformation: roundMetric(baseline.causal_deformation),
    geometric_integrity: roundMetric(baseline.geometric_integrity),
  };
  elements.developerJson.textContent = JSON.stringify(payload, null, 2);
}

function evaluateCurrentScenario() {
  const sourceScenario = state.scenarios[state.activeScenario];
  if (!sourceScenario) {
    return;
  }

  const modifiedScenario =
    isZeroModifier(state.modifiers) ? deepClone(sourceScenario) : applyGenerator(sourceScenario, state.modifiers);

  const evaluation = evaluateScenario(modifiedScenario);
  const coherenceScore = deriveCoherenceScore(evaluation);
  const confidence = deriveConfidence(evaluation);

  state.evaluation = {
    ...evaluation,
    scenarioSeed: state.activeScenario,
    modifiers: { ...state.modifiers },
    coherence_score: coherenceScore,
    confidence,
  };

  renderEvaluation();
  drawSignalCanvas();
}

function renderEvaluation() {
  const result = state.evaluation;
  if (!result) {
    return;
  }

  elements.resultTitle.textContent = `${result.scenarioSeed} -> ${result.classification}`;
  elements.resultSubtitle.textContent = result.description;
  elements.resultBadge.textContent = result.classification;
  elements.resultBadge.className = `oracle-badge ${result.classification}`;

  const statCards = [
    { label: "Classification", value: capitalize(result.classification) },
    { label: "Confidence", value: formatMetric(result.confidence) },
    { label: "Coherence score", value: formatMetric(result.coherence_score) },
    { label: "Scenario", value: result.scenarioSeed },
  ];
  elements.resultStats.replaceChildren(...statCards.map(createStatCard));

  const metrics = [
    {
      label: "Structural retention",
      short: "persistence",
      value: result.structural_retention,
      display: result.structural_retention,
    },
    {
      label: "Temporal consistency",
      short: "ordering",
      value: result.temporal_consistency,
      display: result.temporal_consistency,
    },
    {
      label: "Causal deformation",
      short: "depth drift",
      value: 1 - clamp(result.causal_deformation),
      display: result.causal_deformation,
    },
    {
      label: "Geometric integrity",
      short: "distance coherence",
      value: result.geometric_integrity,
      display: result.geometric_integrity,
    },
    {
      label: "Aggregate coherence",
      short: "derived score",
      value: result.coherence_score,
      display: result.coherence_score,
    },
  ];
  elements.metricBars.replaceChildren(...metrics.map(createMetricBar));

  renderHeatmap();

  const payload = {
    scenario: result.scenario,
    seed: result.scenarioSeed,
    classification: result.classification,
    confidence: roundMetric(result.confidence),
    coherence_score: roundMetric(result.coherence_score),
    structural_retention: roundMetric(result.structural_retention),
    temporal_consistency: roundMetric(result.temporal_consistency),
    causal_deformation: roundMetric(result.causal_deformation),
    geometric_integrity: roundMetric(result.geometric_integrity),
    modifiers: {
      noise: roundMetric(result.modifiers.noise),
      connectivity_drop: roundMetric(result.modifiers.connectivityDrop),
      perturbation_strength: roundMetric(result.modifiers.perturbation),
      temporal_distortion: roundMetric(result.modifiers.temporalDistortion),
    },
  };
  elements.oracleJson.textContent = JSON.stringify(payload, null, 2);
}

function renderHeatmap() {
  const result = state.evaluation;
  if (!result) {
    return;
  }

  const noiseValues = createAxisSeries(result.modifiers.noise, 0.03, 5, 0, 0.15);
  const connectivityValues = createAxisSeries(result.modifiers.connectivityDrop, 0.06, 5, 0, 0.30);
  const sourceScenario = state.scenarios[state.activeScenario];
  const fragment = document.createDocumentFragment();

  for (const noise of noiseValues) {
    for (const connectivityDrop of connectivityValues) {
      const probeScenario = applyGenerator(sourceScenario, {
        ...result.modifiers,
        noise,
        connectivityDrop,
      });
      const probe = evaluateScenario(probeScenario);
      const cell = document.createElement("div");
      const isCurrent = almostEqual(noise, result.modifiers.noise) && almostEqual(connectivityDrop, result.modifiers.connectivityDrop);
      cell.className = `heatmap-cell ${isCurrent ? "current" : ""}`;
      cell.style.background = CLASS_COLORS[probe.classification];
      cell.innerHTML = `
        <span class="heatmap-class">${escapeHtml(probe.classification)}</span>
        <span class="heatmap-score">score ${formatMetric(deriveCoherenceScore(probe))}</span>
        <span class="heatmap-coords">n ${formatMetric(noise)} | c ${formatMetric(connectivityDrop)}</span>
      `;
      fragment.appendChild(cell);
    }
  }

  elements.heatmapGrid.replaceChildren(fragment);
  elements.heatmapCaption.textContent =
    "Current modifiers are outlined. The heatmap sweeps only noise and connectivity so visitors can see the local stability basin around the current state.";
}

function startSignalAnimation() {
  const frame = () => {
    state.animationTick += 0.012;
    drawSignalCanvas();
    state.animationFrameId = requestAnimationFrame(frame);
  };
  if (state.animationFrameId) {
    cancelAnimationFrame(state.animationFrameId);
  }
  state.animationFrameId = requestAnimationFrame(frame);
}

function drawSignalCanvas() {
  const canvas = elements.signalCanvas;
  const result = state.evaluation;
  if (!canvas || !result) {
    return;
  }

  const context = canvas.getContext("2d");
  if (!context) {
    return;
  }

  const rect = canvas.getBoundingClientRect();
  const dpr = window.devicePixelRatio || 1;
  const width = Math.max(320, Math.floor(rect.width * dpr));
  const height = Math.max(240, Math.floor(rect.height * dpr));
  if (canvas.width !== width || canvas.height !== height) {
    canvas.width = width;
    canvas.height = height;
  }

  context.save();
  context.scale(dpr, dpr);
  const drawWidth = width / dpr;
  const drawHeight = height / dpr;
  context.clearRect(0, 0, drawWidth, drawHeight);

  const modifiers = result.modifiers;
  const ordering = result.temporal_consistency;
  const integrity = result.geometric_integrity;
  const deformation = result.causal_deformation;
  const positions = [
    { id: "source", baseX: 0.12, depth: 0 },
    { id: "near", baseX: 0.33, depth: 1 },
    { id: "mid", baseX: 0.58, depth: 2 },
    { id: "far", baseX: 0.84, depth: 3 },
  ].map((node, index) => {
    const xShift = modifiers.connectivityDrop * node.depth * 0.07 - modifiers.perturbation * (node.depth >= 2 ? 0.05 : 0);
    const yWave =
      Math.sin(state.animationTick * 2.5 + index * 0.9) * (10 + 22 * modifiers.noise) +
      modifiers.temporalDistortion * node.depth * 16 +
      modifiers.perturbation * (node.depth === 3 ? 20 : 0);
    return {
      ...node,
      x: (node.baseX + xShift) * drawWidth,
      y: drawHeight * 0.52 + yWave,
    };
  });

  const glow = context.createLinearGradient(0, 0, drawWidth, drawHeight);
  glow.addColorStop(0, "rgba(139, 199, 196, 0.18)");
  glow.addColorStop(1, "rgba(209, 154, 86, 0.06)");
  context.fillStyle = glow;
  context.fillRect(0, 0, drawWidth, drawHeight);

  context.lineWidth = 2;
  for (let index = 0; index < positions.length - 1; index += 1) {
    const left = positions[index];
    const right = positions[index + 1];
    const edgeStrength = clamp((ordering + integrity) / 2 - modifiers.perturbation * 0.24);
    context.strokeStyle = `rgba(139, 199, 196, ${0.20 + edgeStrength * 0.45})`;
    context.beginPath();
    context.moveTo(left.x, left.y);
    context.bezierCurveTo(
      left.x + (right.x - left.x) * 0.34,
      left.y - 26 * (1 - deformation),
      left.x + (right.x - left.x) * 0.66,
      right.y + 26 * deformation,
      right.x,
      right.y,
    );
    context.stroke();
  }

  const pulseTravel = (state.animationTick % 1) * (positions.length - 1);
  for (let index = 0; index < positions.length; index += 1) {
    const node = positions[index];
    const activity = clamp(1 - Math.abs(index - pulseTravel) * 0.7 + result.structural_retention * 0.2);
    const radius = 14 + activity * 7;
    context.beginPath();
    context.fillStyle = `rgba(209, 154, 86, ${0.18 + activity * 0.25})`;
    context.arc(node.x, node.y, radius * 1.8, 0, Math.PI * 2);
    context.fill();

    context.beginPath();
    context.fillStyle = `rgba(139, 199, 196, ${0.70 + activity * 0.18})`;
    context.arc(node.x, node.y, radius, 0, Math.PI * 2);
    context.fill();

    context.fillStyle = "rgba(15, 20, 24, 0.92)";
    context.font = '600 12px "Avenir Next", "IBM Plex Sans", sans-serif';
    context.textAlign = "center";
    context.fillText(node.id, node.x, node.y + 4);
  }

  context.textAlign = "left";
  context.fillStyle = "rgba(243, 239, 232, 0.92)";
  context.font = '600 14px "Avenir Next", "IBM Plex Sans", sans-serif';
  context.fillText(`retention ${formatMetric(result.structural_retention)}`, 18, 28);
  context.fillText(`ordering ${formatMetric(result.temporal_consistency)}`, 18, 50);
  context.fillText(`deformation ${formatMetric(result.causal_deformation)}`, 18, 72);
  context.fillText(`integrity ${formatMetric(result.geometric_integrity)}`, 18, 94);

  context.restore();
}

function applyGenerator(scenario, modifiers) {
  const noise = clamp(modifiers.noise, 0, 0.15);
  const connectivityDrop = clamp(modifiers.connectivityDrop, 0, 0.30);
  const perturbation = clamp(modifiers.perturbation, 0, 1);
  const temporalDistortion = clamp(modifiers.temporalDistortion, 0, 1);

  if (noise <= 0 && connectivityDrop <= 0 && perturbation <= 0 && temporalDistortion <= 0) {
    return deepClone(scenario);
  }

  const generated = deepClone(scenario);
  const modifierBits = [];
  if (noise > 0) {
    modifierBits.push(`noise${noise.toFixed(2)}`);
  }
  if (connectivityDrop > 0) {
    modifierBits.push(`drop${connectivityDrop.toFixed(2)}`);
  }
  if (perturbation > 0) {
    modifierBits.push(`perturb${perturbation.toFixed(2)}`);
  }
  if (temporalDistortion > 0) {
    modifierBits.push(`time${temporalDistortion.toFixed(2)}`);
  }

  generated.name = `${generated.name}_${modifierBits.map((part) => part.replaceAll(".", "p")).join("_")}`;
  generated.description = `${generated.description} Deterministic public modifiers: ${modifierBits.join(", ")}.`;

  const history = generated.persistence_history;
  for (let index = 0; index < history.length; index += 1) {
    const record = history[index];
    const progress = (index + 1) / Math.max(history.length, 1);
    const fragmentationLoad = perturbation * 0.16 * progress;
    const temporalLoad = temporalDistortion * 0.10 * progress;
    record.width_growth = roundMetric(
      toNumber(record.width_growth) +
        noise * 0.20 * progress +
        connectivityDrop * 0.40 * progress +
        fragmentationLoad +
        temporalLoad * 0.35,
    );
    record.concentration_retention = roundMetric(
      clamp(
        toNumber(record.concentration_retention) -
          noise * 0.22 * progress -
          connectivityDrop * 0.35 * progress -
          fragmentationLoad * 1.2 -
          temporalLoad * 0.22,
      ),
    );
    record.participation_growth = roundMetric(
      toNumber(record.participation_growth) +
        noise * 0.25 * progress +
        connectivityDrop * 0.30 * progress +
        fragmentationLoad +
        temporalLoad * 0.28,
    );
    record.overlap = roundMetric(
      clamp(
        toNumber(record.overlap) -
          noise * 0.30 * progress -
          connectivityDrop * 0.35 * progress -
          fragmentationLoad * 1.32 -
          temporalLoad * 0.38,
      ),
    );
    record.recovery_score = roundMetric(
      clamp(
        toNumber(record.recovery_score) -
          noise * 0.34 * progress -
          connectivityDrop * 0.42 * progress -
          fragmentationLoad * 1.18 -
          temporalLoad * 0.45,
      ),
    );
  }

  const expectedShells = generated.expected_shells;
  const maxDepth = Math.max(...Object.values(expectedShells).map((value) => toNumber(value)), 1);

  for (let windowIndex = 0; windowIndex < generated.windows.length; windowIndex += 1) {
    const window = generated.windows[windowIndex];
    const probeHistories = window.probe_histories;
    for (const [node, values] of Object.entries(probeHistories)) {
      if (node === generated.source) {
        continue;
      }
      const depth = toNumber(expectedShells[node] ?? 1);
      let shift = Math.round(connectivityDrop * depth * 2 + temporalDistortion * depth * 1.4);
      if (perturbation > 0.55 && (node === "mid" || node === "far")) {
        shift += 1;
      }

      let attenuation =
        noise * (0.10 + 0.05 * depth) +
        connectivityDrop * (0.10 * depth) +
        perturbation * (0.09 + 0.04 * depth) +
        temporalDistortion * (0.05 + 0.03 * depth);
      if (perturbation > 0.35 && node === "mid") {
        attenuation += 0.05 * perturbation;
      }
      if (perturbation > 0.45 && node === "far") {
        attenuation += 0.11 * perturbation;
      }

      const adjusted = [];
      for (let sampleIndex = 0; sampleIndex < values.length; sampleIndex += 1) {
        const progress = (sampleIndex + 1) / Math.max(values.length, 1);
        const wobble = noise * 0.03 * ((windowIndex + sampleIndex + depth) % 3);
        const timingBleed = temporalDistortion * 0.04 * sampleIndex;
        const degraded = clamp(toNumber(values[sampleIndex]) - attenuation * progress - wobble - timingBleed);
        adjusted.push(roundMetric(degraded));
      }

      if (shift > 0) {
        const padded = Array.from({ length: shift }, () => 0);
        const shifted = padded.concat(adjusted).slice(0, adjusted.length);
        adjusted.splice(0, adjusted.length, ...shifted);
      }

      if (perturbation > 0.62 && node === "far" && adjusted.length >= 3) {
        adjusted[1] = roundMetric(Math.max(adjusted[1], clamp(0.70 - 0.06 * windowIndex)));
        adjusted[2] = roundMetric(Math.max(adjusted[2], clamp(0.81 - 0.05 * windowIndex)));
      }

      probeHistories[node] = adjusted;
    }
  }

  const distanceProxy = generated.distance_proxy;
  const nodes = generated.nodes;
  const baseScale = 1 + noise * 0.10 + connectivityDrop * 0.20 + perturbation * 0.08 + temporalDistortion * 0.03;

  for (const [left, middle, right] of combinationsOfThree(nodes)) {
    void left;
    void middle;
    void right;
  }

  for (let i = 0; i < nodes.length; i += 1) {
    for (let j = i + 1; j < nodes.length; j += 1) {
      const left = nodes[i];
      const right = nodes[j];
      const depthGap = Math.abs(toNumber(expectedShells[left] ?? 0) - toNumber(expectedShells[right] ?? 0));
      let current = toNumber(distanceProxy[left][right]) * baseScale;
      current += connectivityDrop * 0.20 * depthGap / Math.max(maxDepth, 1);
      current += perturbation * 0.16 * depthGap / Math.max(maxDepth, 1);
      current += temporalDistortion * 0.05;
      setSymmetricDistance(distanceProxy, left, right, roundMetric(current));
    }
  }

  if (connectivityDrop > 0) {
    setSymmetricDistance(
      distanceProxy,
      "source",
      "far",
      roundMetric(toNumber(distanceProxy.source.far) + 0.80 * connectivityDrop),
    );
    setSymmetricDistance(
      distanceProxy,
      "near",
      "mid",
      roundMetric(Math.max(0.2, toNumber(distanceProxy.near.mid) - 0.30 * connectivityDrop)),
    );
  }

  if (perturbation > 0) {
    setSymmetricDistance(
      distanceProxy,
      "source",
      "far",
      roundMetric(toNumber(distanceProxy.source.far) + 0.90 * perturbation),
    );
    setSymmetricDistance(
      distanceProxy,
      "near",
      "far",
      roundMetric(Math.max(0.2, toNumber(distanceProxy.near.far) - 0.85 * perturbation)),
    );
    setSymmetricDistance(
      distanceProxy,
      "mid",
      "far",
      roundMetric(toNumber(distanceProxy.mid.far) + 0.30 * perturbation),
    );
  }

  if (temporalDistortion > 0) {
    setSymmetricDistance(
      distanceProxy,
      "source",
      "mid",
      roundMetric(toNumber(distanceProxy.source.mid) + 0.45 * temporalDistortion),
    );
    setSymmetricDistance(
      distanceProxy,
      "near",
      "far",
      roundMetric(toNumber(distanceProxy.near.far) + 0.35 * temporalDistortion),
    );
  }

  return generated;
}

function evaluateScenario(scenario) {
  const nodes = [...scenario.nodes];
  const source = String(scenario.source);
  const expectedShells = Object.fromEntries(
    Object.entries(scenario.expected_shells).map(([key, value]) => [key, toNumber(value)]),
  );
  const windows = arrivalWindows(scenario);

  const referenceArrival = windows[0];
  const referenceEdges = reconstructInfluenceEdges(referenceArrival);
  const referenceDepths = causalDepths(referenceEdges, source);

  const orderingScores = [];
  const depthDrifts = [];
  const shellOverlaps = [];

  for (const arrivalOrdering of windows) {
    const currentEdges = reconstructInfluenceEdges(arrivalOrdering);
    const currentDepths = causalDepths(currentEdges, source);
    const mismatch = orderCompatibility(referenceEdges, arrivalOrdering);
    const completeness = completenessFraction(arrivalOrdering, nodes);
    orderingScores.push(acyclicityScore(currentEdges, nodes) * (1 - mismatch) * completeness);
    depthDrifts.push(normalizedDepthDrift(referenceDepths, currentDepths, nodes));
    shellOverlaps.push(shellOverlapFraction(expectedShells, currentDepths));
  }

  const persistence = normalizePersistence(scenario);
  const ordering = average(orderingScores);
  const depthDrift = average(depthDrifts);
  const shellOverlap = average(shellOverlaps);
  const triangleViolation = triangleViolationRate(scenario.distance_proxy, nodes);
  const distanceCoherence = Math.max(0, 1 - 0.5 * (shellOverlap + triangleViolation));
  const classification = classifyRow(persistence, ordering, depthDrift, distanceCoherence);

  return {
    scenario: String(scenario.name),
    description: String(scenario.description),
    structural_retention: persistence,
    temporal_consistency: ordering,
    causal_deformation: depthDrift,
    geometric_integrity: distanceCoherence,
    persistence,
    ordering,
    depth_drift: depthDrift,
    distance_coherence: distanceCoherence,
    shell_overlap: shellOverlap,
    triangle_violation_rate: triangleViolation,
    classification,
  };
}

function normalizePersistence(scenario) {
  const tauGrid = scenario.tau_grid.map((value) => toNumber(value));
  const value = persistenceTime(scenario.persistence_history, tauGrid, THRESHOLDS);
  if (!tauGrid.length || tauGrid[tauGrid.length - 1] <= 0) {
    return 0;
  }
  return value / tauGrid[tauGrid.length - 1];
}

function persistenceTime(history, tauGrid, thresholds) {
  let lastTau = toNumber(tauGrid[0] ?? 0);
  for (let index = 0; index < history.length; index += 1) {
    const record = history[index];
    const tau = toNumber(tauGrid[index]);
    const widthGrowth = toNumber(record.width_growth);
    const concentration = toNumber(record.concentration_retention);
    const participationGrowth = toNumber(record.participation_growth);
    const overlap = toNumber(record.overlap);
    const score = toNumber(record.recovery_score);

    const persistent =
      widthGrowth <= thresholds.maxWidthGrowth &&
      concentration >= thresholds.minConcentration &&
      participationGrowth <= thresholds.maxParticipationGrowth &&
      overlap >= thresholds.minOverlap &&
      score >= thresholds.minRecoveryScore;

    if (!persistent) {
      break;
    }
    lastTau = tau;
  }
  return lastTau;
}

function arrivalWindows(scenario) {
  const threshold = toNumber(scenario.ordering_threshold);
  return scenario.windows.map((window) => frontArrivalOrder(window.probe_histories, threshold));
}

function frontArrivalOrder(probeHistories, threshold) {
  const entries = Object.entries(probeHistories).sort(([left], [right]) => left.localeCompare(right));
  return Object.fromEntries(
    entries.map(([probeName, history]) => [probeName, firstThresholdCrossing(history, threshold)]),
  );
}

function firstThresholdCrossing(signal, threshold) {
  for (let index = 0; index < signal.length; index += 1) {
    if (toNumber(signal[index]) >= threshold - 1.0e-12) {
      return index;
    }
  }
  return null;
}

function reconstructInfluenceEdges(eventTimes) {
  const ordered = Object.entries(eventTimes)
    .sort((left, right) => {
      const leftValue = left[1] === null ? Number.POSITIVE_INFINITY : left[1];
      const rightValue = right[1] === null ? Number.POSITIVE_INFINITY : right[1];
      if (leftValue !== rightValue) {
        return leftValue - rightValue;
      }
      return left[0].localeCompare(right[0]);
    })
    .filter((entry) => entry[1] !== null);

  const edges = [];
  for (let index = 0; index < ordered.length - 1; index += 1) {
    edges.push([ordered[index][0], ordered[index + 1][0]]);
  }
  return edges;
}

function acyclicityScore() {
  return 1;
}

function causalDepths(edges, source) {
  const adjacency = new Map();
  for (const [src, dst] of edges) {
    if (!adjacency.has(src)) {
      adjacency.set(src, new Set());
    }
    adjacency.get(src).add(dst);
  }

  const depths = { [source]: 0 };
  const queue = [source];
  while (queue.length) {
    const node = queue.shift();
    const nextNodes = [...(adjacency.get(node) || [])].sort();
    for (const nextNode of nextNodes) {
      if (!(nextNode in depths)) {
        depths[nextNode] = depths[node] + 1;
        queue.push(nextNode);
      }
    }
  }
  return depths;
}

function orderCompatibility(edgeSet, arrivalOrdering) {
  let comparable = 0;
  let mismatches = 0;
  for (const [src, dst] of edgeSet) {
    const srcTime = arrivalOrdering[src];
    const dstTime = arrivalOrdering[dst];
    if (srcTime === null || dstTime === null) {
      continue;
    }
    comparable += 1;
    if (srcTime > dstTime) {
      mismatches += 1;
    }
  }
  return comparable === 0 ? 0 : mismatches / comparable;
}

function completenessFraction(arrivalOrdering, nodes) {
  const arrived = nodes.filter((node) => arrivalOrdering[node] !== null).length;
  return arrived / Math.max(nodes.length, 1);
}

function normalizedDepthDrift(referenceDepths, currentDepths, nodes) {
  const maxDepth = Math.max(0, ...Object.values(referenceDepths));
  const missingDepth = maxDepth + 1;
  const diffs = nodes.map((node) => {
    const current = currentDepths[node] ?? missingDepth;
    const reference = referenceDepths[node] ?? missingDepth;
    return Math.abs(current - reference);
  });
  return average(diffs) / Math.max(missingDepth, 1);
}

function shellOverlapFraction(expectedShells, currentDepths) {
  const maxShell = Math.max(0, ...Object.values(expectedShells));
  const missingShell = maxShell + 1;
  let mismatches = 0;
  for (const [node, expected] of Object.entries(expectedShells)) {
    if ((currentDepths[node] ?? missingShell) !== expected) {
      mismatches += 1;
    }
  }
  return mismatches / Math.max(Object.keys(expectedShells).length, 1);
}

function triangleViolationRate(distanceProxy, nodes) {
  const triples = combinationsOfThree(nodes);
  if (!triples.length) {
    return 0;
  }

  let violations = 0;
  for (const [left, middle, right] of triples) {
    const distances = [
      toNumber(distanceProxy[left][middle]),
      toNumber(distanceProxy[left][right]),
      toNumber(distanceProxy[middle][right]),
    ].sort((a, b) => a - b);
    if (distances[2] > distances[0] + distances[1] + 1.0e-9) {
      violations += 1;
    }
  }
  return violations / triples.length;
}

function classifyRow(persistence, ordering, depthDrift, distanceCoherence) {
  if (persistence >= 0.9 && ordering >= 0.9 && depthDrift <= 0.05 && distanceCoherence >= 0.9) {
    return "stable";
  }
  if (persistence >= 0.65 && ordering >= 0.75 && depthDrift <= 0.15 && distanceCoherence >= 0.6) {
    return "marginal";
  }
  return "unstable";
}

function deriveCoherenceScore(result) {
  return clamp(
    result.structural_retention * 0.34 +
      result.temporal_consistency * 0.24 +
      (1 - clamp(result.causal_deformation)) * 0.22 +
      result.geometric_integrity * 0.20,
  );
}

function deriveConfidence(result) {
  if (result.classification === "stable") {
    const margin = Math.min(
      result.structural_retention - 0.90,
      result.temporal_consistency - 0.90,
      0.05 - result.causal_deformation,
      result.geometric_integrity - 0.90,
    );
    return clamp(0.78 + margin * 2.3, 0.55, 0.99);
  }

  if (result.classification === "marginal") {
    const margin = Math.min(
      result.structural_retention - 0.65,
      result.temporal_consistency - 0.75,
      0.15 - result.causal_deformation,
      result.geometric_integrity - 0.60,
    );
    return clamp(0.58 + margin * 1.9, 0.42, 0.9);
  }

  const deficit = Math.max(
    0.65 - result.structural_retention,
    0.75 - result.temporal_consistency,
    result.causal_deformation - 0.15,
    0.60 - result.geometric_integrity,
  );
  return clamp(0.55 + deficit * 0.8, 0.42, 0.97);
}

function createStatCard(stat) {
  const card = document.createElement("div");
  card.className = "stat-card";
  card.innerHTML = `<strong>${escapeHtml(stat.label)}</strong><span>${escapeHtml(String(stat.value))}</span>`;
  return card;
}

function createMetricBar(metric) {
  const card = document.createElement("div");
  card.className = "metric-bar";
  card.innerHTML = `
    <div class="metric-text">
      <strong>${escapeHtml(metric.label)}</strong>
      <span>${escapeHtml(metric.short)} ${formatMetric(metric.display)}</span>
    </div>
    <div class="metric-track">
      <div class="metric-fill" style="width: ${clamp(metric.value) * 100}%"></div>
    </div>
  `;
  return card;
}

function createListItem(text) {
  const item = document.createElement("li");
  item.textContent = text;
  return item;
}

function createReferenceLink(ref) {
  const anchor = document.createElement("a");
  anchor.href = ref.href;
  anchor.target = "_blank";
  anchor.rel = "noreferrer";
  anchor.textContent = ref.label;
  return anchor;
}

function setSymmetricDistance(matrix, left, right, value) {
  matrix[left][right] = value;
  matrix[right][left] = value;
}

function createAxisSeries(center, step, count, min, max) {
  const span = step * (count - 1);
  const start = clamp(center - (step * Math.floor(count / 2)), min, Math.max(min, max - span));
  const values = [];
  for (let index = 0; index < count; index += 1) {
    values.push(roundMetric(start + index * step));
  }
  return values;
}

function combinationsOfThree(nodes) {
  const output = [];
  for (let i = 0; i < nodes.length; i += 1) {
    for (let j = i + 1; j < nodes.length; j += 1) {
      for (let k = j + 1; k < nodes.length; k += 1) {
        output.push([nodes[i], nodes[j], nodes[k]]);
      }
    }
  }
  return output;
}

function syncModifierInputs() {
  elements.noiseSlider.value = state.modifiers.noise;
  elements.connectivitySlider.value = state.modifiers.connectivityDrop;
  elements.perturbationSlider.value = state.modifiers.perturbation;
  elements.temporalSlider.value = state.modifiers.temporalDistortion;
  elements.noiseValue.textContent = formatMetric(state.modifiers.noise);
  elements.connectivityValue.textContent = formatMetric(state.modifiers.connectivityDrop);
  elements.perturbationValue.textContent = formatMetric(state.modifiers.perturbation);
  elements.temporalValue.textContent = formatMetric(state.modifiers.temporalDistortion);
}

function isZeroModifier(modifiers) {
  return (
    modifiers.noise <= 0 &&
    modifiers.connectivityDrop <= 0 &&
    modifiers.perturbation <= 0 &&
    modifiers.temporalDistortion <= 0
  );
}

function deepClone(value) {
  return JSON.parse(JSON.stringify(value));
}

function average(values) {
  if (!values.length) {
    return 0;
  }
  return values.reduce((sum, value) => sum + value, 0) / values.length;
}

function clamp(value, min = 0, max = 1) {
  return Math.max(min, Math.min(max, toNumber(value)));
}

function roundMetric(value) {
  return Math.round(toNumber(value) * 1_000_000) / 1_000_000;
}

function formatMetric(value) {
  return toNumber(value).toFixed(2);
}

function formatFixed(value, digits) {
  return toNumber(value).toFixed(digits);
}

function capitalize(value) {
  return value.charAt(0).toUpperCase() + value.slice(1);
}

function almostEqual(left, right) {
  return Math.abs(left - right) < 1.0e-9;
}

function toNumber(value) {
  return Number.parseFloat(String(value)) || 0;
}

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}
