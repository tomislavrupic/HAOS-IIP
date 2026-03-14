import * as THREE from "./vendor/three.module.js";

function clamp01(value) {
  return Math.max(0, Math.min(1, value));
}

function pointVector(point) {
  return new THREE.Vector3(point.x, point.y, point.z);
}

class AtlasThreeGridScene {
  constructor() {
    this.canvas = null;
    this.renderer = null;
    this.scene = null;
    this.camera = null;
    this.root = null;
    this.ghostPoints = null;
    this.activePoints = null;
    this.edgeLines = null;
    this.pulsePoints = null;
    this.trailLine = null;
    this.centerMesh = null;
    this.envelopeMesh = null;
    this.boundsMesh = null;
  }

  ensure(canvas) {
    if (this.renderer && this.canvas === canvas) {
      return;
    }

    this.canvas = canvas;
    this.renderer = new THREE.WebGLRenderer({
      canvas,
      alpha: true,
      antialias: true,
      powerPreference: "high-performance",
    });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.renderer.setClearColor(0x000000, 0);

    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(42, 1, 0.1, 20);
    this.camera.position.set(0, 0, 3.35);
    this.camera.lookAt(0, 0, 0);

    this.root = new THREE.Group();
    this.scene.add(this.root);

    const ambient = new THREE.AmbientLight(0xffffff, 0.9);
    const key = new THREE.DirectionalLight(0xa4bfd1, 0.32);
    key.position.set(1.4, 2.1, 2.8);
    this.scene.add(ambient, key);

    this.boundsMesh = new THREE.LineSegments(
      new THREE.EdgesGeometry(new THREE.BoxGeometry(1, 1, 1)),
      new THREE.LineBasicMaterial({
        color: 0x5f7383,
        transparent: true,
        opacity: 0.2,
      }),
    );
    this.root.add(this.boundsMesh);

    this.ghostPoints = new THREE.Points(
      new THREE.BufferGeometry(),
      new THREE.PointsMaterial({
        color: 0x7d90a2,
        size: 0.018,
        transparent: true,
        opacity: 0.22,
        depthWrite: false,
      }),
    );

    this.activePoints = new THREE.Points(
      new THREE.BufferGeometry(),
      new THREE.PointsMaterial({
        size: 0.046,
        vertexColors: true,
        transparent: true,
        opacity: 0.95,
        depthWrite: false,
        blending: THREE.AdditiveBlending,
      }),
    );

    this.edgeLines = new THREE.LineSegments(
      new THREE.BufferGeometry(),
      new THREE.LineBasicMaterial({
        vertexColors: true,
        transparent: true,
        opacity: 0.64,
      }),
    );

    this.pulsePoints = new THREE.Points(
      new THREE.BufferGeometry(),
      new THREE.PointsMaterial({
        size: 0.062,
        vertexColors: true,
        transparent: true,
        opacity: 0.98,
        depthWrite: false,
        blending: THREE.AdditiveBlending,
      }),
    );

    this.trailLine = new THREE.Line(
      new THREE.BufferGeometry(),
      new THREE.LineBasicMaterial({
        color: 0x7ea2a8,
        transparent: true,
        opacity: 0.48,
      }),
    );

    this.centerMesh = new THREE.Mesh(
      new THREE.SphereGeometry(0.04, 18, 18),
      new THREE.MeshBasicMaterial({
        color: 0xf0f7fb,
      }),
    );

    this.envelopeMesh = new THREE.Mesh(
      new THREE.SphereGeometry(0.08, 18, 18),
      new THREE.MeshBasicMaterial({
        color: 0x7ea2a8,
        transparent: true,
        opacity: 0.08,
        wireframe: true,
      }),
    );

    this.root.add(
      this.ghostPoints,
      this.activePoints,
      this.edgeLines,
      this.pulsePoints,
      this.trailLine,
      this.centerMesh,
      this.envelopeMesh,
    );
  }

  resize() {
    const width = this.canvas.clientWidth || this.canvas.width || 1;
    const height = this.canvas.clientHeight || this.canvas.height || 1;
    this.renderer.setSize(width, height, false);
    this.camera.aspect = width / Math.max(height, 1);
    this.camera.updateProjectionMatrix();
  }

  replaceGeometry(object, positions, colors = null) {
    if (object.geometry) {
      object.geometry.dispose();
    }

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
    if (colors) {
      geometry.setAttribute("color", new THREE.Float32BufferAttribute(colors, 3));
    }

    object.geometry = geometry;
    object.visible = positions.length > 0;
  }

  updateNodes(frame) {
    const activeColor = new THREE.Color(frame.regimeFill);
    const neutralColor = new THREE.Color("#9db1c4");
    const ghostPositions = [];
    const activePositions = [];
    const activeColors = [];

    frame.nodes.forEach((node) => {
      if (!node.sampleVisible) {
        return;
      }

      if (node.intensity > 0.08) {
        const color = neutralColor.clone().lerp(activeColor, clamp01(node.intensity * 1.2));
        color.multiplyScalar(0.76 + node.intensity * 0.48);
        activePositions.push(node.point.x, node.point.y, node.point.z);
        activeColors.push(color.r, color.g, color.b);
        return;
      }

      ghostPositions.push(node.point.x, node.point.y, node.point.z);
    });

    this.ghostPoints.material.size = 0.014 + frame.currentWidth * 0.02;
    this.activePoints.material.size = 0.03 + frame.energyNorm * 0.035;
    this.replaceGeometry(this.ghostPoints, ghostPositions);
    this.replaceGeometry(this.activePoints, activePositions, activeColors);
  }

  updateEdges(frame) {
    const coolColor = new THREE.Color("#7ea2a8");
    const warmColor = new THREE.Color("#b49463");
    const positions = [];
    const colors = [];

    frame.activeEdges.forEach((edge) => {
      const from = frame.nodes[edge.startIndex];
      const to = frame.nodes[edge.endIndex];
      const color = edge.signedAlignment >= 0 ? coolColor.clone() : warmColor.clone();
      color.multiplyScalar(0.34 + edge.flowStrength * 1.16);

      positions.push(from.point.x, from.point.y, from.point.z);
      positions.push(to.point.x, to.point.y, to.point.z);
      colors.push(color.r, color.g, color.b);
      colors.push(color.r, color.g, color.b);
    });

    this.edgeLines.material.opacity = 0.28 + frame.energyNorm * 0.34;
    this.replaceGeometry(this.edgeLines, positions, colors);
  }

  updatePulses(frame) {
    const coolColor = new THREE.Color("#f0f7fb");
    const warmColor = new THREE.Color("#f3d9b0");
    const positions = [];
    const colors = [];

    frame.activeEdges.forEach((edge) => {
      const from = frame.nodes[edge.startIndex].point;
      const to = frame.nodes[edge.endIndex].point;
      const phase = edge.signedAlignment >= 0 ? frame.pulsePhaseSeed : 1 - frame.pulsePhaseSeed;
      const color = edge.signedAlignment >= 0 ? coolColor.clone() : warmColor.clone();
      color.multiplyScalar(0.68 + edge.flowStrength * 0.54);

      positions.push(
        from.x + (to.x - from.x) * phase,
        from.y + (to.y - from.y) * phase,
        from.z + (to.z - from.z) * phase,
      );
      colors.push(color.r, color.g, color.b);
    });

    this.pulsePoints.material.size = 0.04 + frame.currentWidth * 0.05;
    this.replaceGeometry(this.pulsePoints, positions, colors);
  }

  updateTrail(frame) {
    const positions = [];

    frame.centerTrail.forEach((sample) => {
      positions.push(sample.point.x, sample.point.y, sample.point.z);
    });

    this.trailLine.material.color.set(frame.regimeFill);
    this.replaceGeometry(this.trailLine, positions);
  }

  updateCenter(frame) {
    this.centerMesh.position.copy(pointVector(frame.currentCenterPoint));
    this.centerMesh.material.color.set(frame.regimeFill);
    this.centerMesh.scale.setScalar(0.72 + frame.currentWidth * 2.4);

    this.envelopeMesh.position.copy(pointVector(frame.currentCenterPoint));
    this.envelopeMesh.material.color.set(frame.regimeFill);
    this.envelopeMesh.scale.setScalar(1 + frame.currentWidth * 5.2);

    this.boundsMesh.material.opacity = 0.12 + frame.energyNorm * 0.12;
  }

  render({ canvas, frame, viewState }) {
    this.ensure(canvas);
    this.resize();

    this.root.rotation.x = viewState.rotationX;
    this.root.rotation.y = viewState.rotationY;
    this.camera.position.z = 3.8 / Math.max(viewState.zoom || 1, 0.35);
    this.camera.lookAt(0, 0, 0);

    this.updateNodes(frame);
    this.updateEdges(frame);
    this.updatePulses(frame);
    this.updateTrail(frame);
    this.updateCenter(frame);
    this.renderer.render(this.scene, this.camera);
    return true;
  }

  clear(canvas) {
    if (!this.renderer || (canvas && this.canvas !== canvas)) {
      return;
    }

    this.renderer.clear();
  }
}

const scene = new AtlasThreeGridScene();

window.AtlasThreeGridRenderer = {
  available: true,
  render(payload) {
    if (!this.available) {
      return false;
    }

    try {
      return scene.render(payload);
    } catch (error) {
      this.available = false;
      console.warn("Atlas Explorer Three.js renderer disabled. Falling back to canvas.", error);
      return false;
    }
  },
  clear(canvas) {
    scene.clear(canvas);
  },
};

window.dispatchEvent(new Event("atlas-three-ready"));
