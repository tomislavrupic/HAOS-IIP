from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from experiments.cst.cst_io import artifact_hashes, read_json, stable_hash
from experiments.cst.run_cst_benchmark import DEFAULT_OUTPUT_DIR, PHASE_ARTIFACTS, run_cst_benchmark


REQUIRED_OUTPUTS = [
    "cst_runs.json",
    "closure_signatures.json",
    "closure_distance_matrix.csv",
    "closure_distance_components.csv",
    "recoverability_vectors.csv",
    "control_distributions.csv",
    "seed_manifest.json",
    "benchmark_result.json",
    "benchmark_report.md",
]


def fail(reason: str) -> None:
    raise SystemExit(f"CST bundle check failed: {reason}")


def hash_payload(payload: dict[str, Any]) -> str:
    comparable = dict(payload)
    comparable.pop("result_hash", None)
    return stable_hash(comparable, prefix="cst_result_")


def main() -> int:
    output_dir = DEFAULT_OUTPUT_DIR
    missing = [name for name in REQUIRED_OUTPUTS if not (output_dir / name).exists()]
    if missing:
        fail(f"missing output files: {missing}")

    disk_result = read_json(output_dir / "benchmark_result.json")
    expected_hash = hash_payload(disk_result)
    if disk_result.get("result_hash") != expected_hash:
        fail("benchmark_result.json result_hash does not match payload")

    memory = run_cst_benchmark(write_outputs=False)
    memory_payload = memory["benchmark_payload"]
    if memory_payload.get("result_hash") != disk_result.get("result_hash"):
        fail("fresh benchmark result_hash differs from benchmark_result.json")

    seed_manifest = read_json(output_dir / "seed_manifest.json")
    current_hashes = artifact_hashes(PHASE_ARTIFACTS)
    if seed_manifest.get("source_artifact_hashes") != current_hashes:
        fail("seed_manifest source_artifact_hashes differ from current frozen inputs")

    if not memory["target_distances"]:
        fail("no target target-pair distances were computed")
    if not memory["control_distances"]:
        fail("no target control-pair distances were computed")

    observations = memory["observations"]
    target_observations = [item for item in observations if item.observation_kind == "target"]
    if not target_observations:
        fail("no target observations were produced")
    for observation in target_observations:
        if observation.provenance.run_instance_id == observation.closure_signature.branch_signature_id:
            fail("run_instance_id and branch_signature_id must remain separate")

    print(
        json.dumps(
            {
                "status": "ok",
                "verdict": disk_result.get("verdict"),
                "result_hash": disk_result.get("result_hash"),
                "observations": len(observations),
                "target_distances": len(memory["target_distances"]),
                "control_distances": len(memory["control_distances"]),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
