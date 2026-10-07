import copy
import json
import unittest

from scripts.check_claim_contracts import CONTRACT_DIR, validate_contract


class ClaimContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.source = CONTRACT_DIR / "scale_bridge_66_5.claim.json"
        cls.contract = json.loads(cls.source.read_text(encoding="utf-8"))

    def test_current_contract_is_valid(self) -> None:
        validate_contract(self.contract, self.source)

    def test_missing_perturbation_domain_fails(self) -> None:
        candidate = copy.deepcopy(self.contract)
        candidate["perturbation_domain"]["included"] = []
        with self.assertRaisesRegex(ValueError, "too few items"):
            validate_contract(candidate, self.source)

    def test_external_support_requires_independent_replication(self) -> None:
        candidate = copy.deepcopy(self.contract)
        candidate["bridge_status"] = "EXTERNALLY_SUPPORTED"
        with self.assertRaisesRegex(ValueError, "requires independent replication"):
            validate_contract(candidate, self.source)

    def test_partial_selection_history_caps_scaling_pass(self) -> None:
        candidate = copy.deepcopy(self.contract)
        candidate["overall_status"] = "PASS"
        with self.assertRaisesRegex(ValueError, "caps non-internal claims at OPEN"):
            validate_contract(candidate, self.source)


if __name__ == "__main__":
    unittest.main()
