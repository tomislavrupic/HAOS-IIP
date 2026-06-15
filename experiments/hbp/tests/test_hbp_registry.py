from __future__ import annotations

import copy
import unittest

from experiments.hbp.hbp_validation import assess_contract, contract_from_dict
from experiments.hbp.pb01_network_recovery import run_pb01_network_recovery as pb01
from experiments.hbp.run_hbp_registry import build_registry


def contract_payload() -> dict:
    return {
        "bridge_id": "TEST-BRIDGE",
        "domain": "synthetic",
        "classification": "PREDICTIVE_BRIDGE",
        "state_schema": {
            "variables": ["x"],
            "state_space": "[0, 1]",
            "units": {"x": "dimensionless"},
            "valid_ranges": {"x": "0..1"},
            "missing_data_policy": "invalid",
        },
        "units": {"x": "dimensionless"},
        "dynamics": {
            "update_law": "x_(t+1)=x_t",
            "parameters": {"alpha": "dimensionless"},
            "perturbations": ["identity"],
            "time_axis": "normalized_step",
            "boundary_conditions": "fixed synthetic boundary",
        },
        "boundary_conditions": "fixed synthetic boundary",
        "symmetries": ["label permutation is a control"],
        "haos_mapping": [
            {
                "haos_source": "frozen_feature",
                "domain_target": "recovery_quality",
                "mapping_function": "predeclared calibration",
                "units_before": "dimensionless",
                "units_after": "dimensionless",
                "semantic_justification": "test calibrated mapping",
                "mapping_status": "CALIBRATED",
                "uncertainty": "bootstrap",
                "failure_conditions": ["holdout failure"],
            }
        ],
        "observation_map": [
            {
                "observable": "recovery_quality",
                "measurement_rule": "bounded residual",
                "units": "dimensionless",
                "valid_range": "0..1",
                "target_role": "prediction target",
            }
        ],
        "prediction_target": "recovery_quality",
        "calibration_split": "calibration",
        "holdout_split": "holdout",
        "controls": [
            {
                "name": "null_control",
                "purpose": "destroy test signal",
                "preserves": ["sample count"],
                "destroys": ["target relation"],
                "invalid_if": "target relation survives",
            }
        ],
        "baselines": [
            {
                "name": "mean_predictor",
                "family": "mean",
                "prediction_rule": "calibration mean",
                "training_scope": "calibration",
                "leakage_risk": "low",
            }
        ],
        "uncertainty": "bootstrap CI",
        "falsification": {
            "fail_conditions": ["not better than baseline"],
            "invalidation_conditions": ["target leakage"],
            "underpowered_policy": "return UNDERPOWERED",
            "leakage_policy": "return TARGET_LEAKAGE_DETECTED",
        },
        "verdict_logic": "requires replication; no intervention-backed unique generative mechanism",
        "provenance": {
            "source_artifacts": [],
            "code_artifacts": [],
            "precommitment_artifact": None,
            "external_data_status": "synthetic_only",
            "level_history": [
                "FORMAL_CORRESPONDENCE",
                "ANALOGY_BRIDGE",
                "OPERATIONAL_MAPPING",
            ],
            "notes": ["test contract"],
        },
    }


class HBPRegistryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.pb01 = pb01.run_pb01(write_outputs=False)
        cls.registry = build_registry(write_outputs=False)

    def test_contract_field_enforcement(self) -> None:
        payload = contract_payload()
        payload["state_schema"] = None
        assessment = assess_contract(contract_from_dict(payload))
        self.assertEqual(assessment.effective_classification, "FORMAL_CORRESPONDENCE")
        self.assertIn("CONTRACT_INCOMPLETE", assessment.labels)
        self.assertIn("state_schema", assessment.missing_required_fields)

    def test_classification_ceiling_for_heuristic_mapping(self) -> None:
        payload = contract_payload()
        payload["classification"] = "PREDICTIVE_BRIDGE"
        payload["haos_mapping"][0]["mapping_status"] = "HEURISTIC"
        assessment = assess_contract(contract_from_dict(payload))
        self.assertEqual(assessment.effective_classification, "OPERATIONAL_MAPPING")
        self.assertIn("OPERATIONAL_MAPPING_PARTIAL", assessment.labels)

    def test_forced_mapping_failure(self) -> None:
        payload = contract_payload()
        payload["haos_mapping"][0]["mapping_status"] = "POST_HOC"
        assessment = assess_contract(contract_from_dict(payload))
        self.assertEqual(assessment.effective_classification, "FORMAL_CORRESPONDENCE")
        self.assertIn("MAPPING_INVALID", assessment.labels)

    def test_unit_validation_failure(self) -> None:
        payload = contract_payload()
        payload["haos_mapping"][0]["units_after"] = ""
        assessment = assess_contract(contract_from_dict(payload))
        self.assertEqual(assessment.effective_classification, "FORMAL_CORRESPONDENCE")
        self.assertIn("DIMENSIONAL_ANALYSIS_FAIL", assessment.labels)

    def test_no_level_skipping(self) -> None:
        payload = contract_payload()
        payload["classification"] = "PREDICTIVE_BRIDGE"
        payload["provenance"]["level_history"] = ["FORMAL_CORRESPONDENCE"]
        assessment = assess_contract(contract_from_dict(payload))
        self.assertEqual(assessment.effective_classification, "FORMAL_CORRESPONDENCE")
        self.assertIn("requested classification skips a bridge level", " ".join(assessment.reasons))

    def test_pb01_deterministic_hashes(self) -> None:
        rerun = pb01.run_pb01(write_outputs=False)
        self.assertEqual(self.pb01["result"]["result_hash"], rerun["result"]["result_hash"])
        self.assertEqual(self.registry["hbp_result"]["pb01_status"], self.pb01["result"]["verdict"]["status"])

    def test_frozen_holdout_isolation(self) -> None:
        contract = self.pb01["contract"]
        dev = set(contract["seeds"]["development"])
        calibration = set(contract["seeds"]["calibration"])
        holdout = set(contract["seeds"]["holdout"])
        self.assertFalse(dev & holdout)
        self.assertFalse(calibration & holdout)
        self.assertGreaterEqual(self.pb01["result"]["case_counts"]["holdout"], contract["verdict_logic"]["min_holdout_cases"])

    def test_baseline_and_control_routing(self) -> None:
        models = {row["model"] for row in self.pb01["baseline_rows"]}
        required_models = {
            "mean_predictor",
            "random_predictor",
            "graph_size_density",
            "degree_centrality",
            "betweenness_centrality",
            "closeness_centrality",
            "pagerank_centrality",
            "eigenvector_centrality",
            "shortest_path_to_perturbation",
            "algebraic_connectivity",
            "domain_diffusion_early_probe",
            "supervised_graph_spectral_model",
            "haos_ablated_score",
            "haos_calibrated_recovery_score",
        }
        self.assertTrue(required_models.issubset(models))
        controls = {row["control"] for row in self.pb01["control_rows"]}
        self.assertIn("intentional_leakage_positive_control", controls)
        self.assertIn("topology_destroyed_graph", controls)
        self.assertIn("parameter_matched_null", controls)

    def test_metric_correctness(self) -> None:
        self.assertAlmostEqual(pb01.spearman([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]), 1.0)
        self.assertAlmostEqual(pb01.spearman([1.0, 2.0, 3.0], [3.0, 2.0, 1.0]), -1.0)
        self.assertAlmostEqual(pb01.roc_auc([0.1, 0.9], [0, 1]) or 0.0, 1.0)
        self.assertAlmostEqual(pb01.pr_auc([0.1, 0.9], [0, 1]) or 0.0, 1.0)

    def test_forced_baseline_non_separation(self) -> None:
        contract = pb01.precommitment_payload()
        rows = [
            {"split": "calibration", "model": pb01.PRIMARY_HAOS_MODEL, "n": 120, "spearman": 0.51},
            {"split": "holdout", "model": pb01.PRIMARY_HAOS_MODEL, "n": 120, "spearman": 0.50},
            {"split": "holdout", "model": "supervised_graph_spectral_model", "n": 120, "spearman": 0.49},
        ]
        controls = [{"control": "intentional_leakage_positive_control", "status": "TARGET_LEAKAGE_DETECTED"}]
        verdict = pb01.verdict_from_metrics(rows, controls, contract)
        self.assertEqual(verdict["status"], "PREDICTION_NOT_DISTINCT_FROM_BASELINES")
        self.assertIn("PREDICTION_NOT_DISTINCT_FROM_BASELINES", verdict["labels"])

    def test_forced_underpowered_result(self) -> None:
        contract = pb01.precommitment_payload()
        min_cases = int(contract["verdict_logic"]["min_holdout_cases"])
        rows = [
            {"split": "calibration", "model": pb01.PRIMARY_HAOS_MODEL, "n": min_cases, "spearman": 0.70},
            {"split": "holdout", "model": pb01.PRIMARY_HAOS_MODEL, "n": min_cases - 1, "spearman": 0.90},
            {"split": "holdout", "model": "supervised_graph_spectral_model", "n": min_cases - 1, "spearman": 0.20},
        ]
        controls = [{"control": "intentional_leakage_positive_control", "status": "TARGET_LEAKAGE_DETECTED"}]
        verdict = pb01.verdict_from_metrics(rows, controls, contract)
        self.assertIn("UNDERPOWERED", verdict["labels"])
        self.assertEqual(verdict["status"], "PREDICTION_NOT_DISTINCT_FROM_BASELINES")

    def test_leakage_positive_control_detected(self) -> None:
        leakage = next(row for row in self.pb01["control_rows"] if row["control"] == "intentional_leakage_positive_control")
        self.assertEqual(leakage["status"], "TARGET_LEAKAGE_DETECTED")

    def test_ambiguous_control_failure_remains_open(self) -> None:
        labels = self.pb01["result"]["verdict"]["labels"]
        self.assertIn("MIXED_OPEN", labels)
        self.assertIn("PREDICTION_NOT_DISTINCT_FROM_BASELINES", labels)
        self.assertIn("CONTROL_INVALID", labels)

    def test_registry_does_not_upgrade_existing_bridges(self) -> None:
        rows = {row["bridge_id"]: row for row in self.registry["registry_payload"]["registry_rows"]}
        self.assertEqual(rows["BELL-BRIDGE-V0-3-TERMINAL"]["effective_classification"], "FORMAL_CORRESPONDENCE")
        self.assertNotIn(
            rows["SYNTHETIC-RELATIONAL-CALIBRATION"]["effective_classification"],
            {"EMPIRICALLY_SUPPORTED_BRIDGE", "PHYSICAL_MECHANISM_CANDIDATE"},
        )
        self.assertNotIn(
            rows["PB-01-NETWORK-RECOVERY"]["effective_classification"],
            {"EMPIRICALLY_SUPPORTED_BRIDGE", "PHYSICAL_MECHANISM_CANDIDATE"},
        )

    def test_scalar_claims_do_not_override_official_verdict(self) -> None:
        result = copy.deepcopy(self.pb01["result"])
        self.assertGreater(result["verdict"]["haos_holdout_spearman"], result["verdict"]["best_baseline_holdout_spearman"])
        self.assertEqual(result["verdict"]["status"], "PREDICTION_NOT_DISTINCT_FROM_BASELINES")


if __name__ == "__main__":
    unittest.main()
