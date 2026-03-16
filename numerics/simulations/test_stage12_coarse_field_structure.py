import sys
from unittest.mock import MagicMock

# Mock dependencies to allow importing stage12_coarse_field_structure without numpy/scipy
sys.modules['numpy'] = MagicMock()
sys.modules['scipy'] = MagicMock()
sys.modules['scipy.sparse'] = MagicMock()
sys.modules['scipy.sparse.linalg'] = MagicMock()
sys.modules['matplotlib'] = MagicMock()
sys.modules['matplotlib.pyplot'] = MagicMock()
sys.modules['stage10_common'] = MagicMock()
sys.modules['stage11_collective_wave_interaction'] = MagicMock()

import unittest
# Add current directory to path so we can import from numerics
import os
sys.path.append(os.getcwd())

from numerics.simulations.stage12_coarse_field_structure import classify_coarse

class TestClassifyCoarse(unittest.TestCase):
    def test_multi_basin_fluctuating(self):
        # Condition: mean_basin_count_2 > 1.05
        summary = {
            'mean_basin_count_sigma2': 1.1,
            'mean_dominant_area_fraction_sigma2': 0.1,
            'mean_dominant_area_fraction_sigma4': 0.1,
            'basin_lifetime_sigma4': 10.0,
            'coarse_persistence_sigma4': 1.0,
            'mean_envelope_variance_sigma4': 0.1,
            't_final': 100.0
        }
        self.assertEqual(classify_coarse(summary), 'multi-basin fluctuating regime')

    def test_basin_dominated_primary(self):
        # Condition:
        # mean_basin_count_2 <= 1.05
        # AND basin_lifetime >= 0.45 * t_final
        # AND mean_dom_frac_2 <= 0.06
        # AND mean_dom_frac_4 <= 0.35
        # AND persistence_4 >= 0.9
        # AND env_var_4 >= 0.02
        summary = {
            'mean_basin_count_sigma2': 1.0,
            'mean_dominant_area_fraction_sigma2': 0.05,
            'mean_dominant_area_fraction_sigma4': 0.3,
            'basin_lifetime_sigma4': 50.0,
            'coarse_persistence_sigma4': 0.95,
            'mean_envelope_variance_sigma4': 0.03,
            't_final': 100.0
        }
        self.assertEqual(classify_coarse(summary), 'basin-dominated regime')

    def test_diffuse_coarse_field_dom_frac(self):
        # Condition:
        # mean_basin_count_2 <= 1.05
        # AND NOT basin_dominated_primary
        # AND (mean_dom_frac_4 >= 0.35 OR env_var_4 < 0.02)
        summary = {
            'mean_basin_count_sigma2': 1.0,
            'mean_dominant_area_fraction_sigma2': 0.1,
            'mean_dominant_area_fraction_sigma4': 0.4,
            'basin_lifetime_sigma4': 10.0,
            'coarse_persistence_sigma4': 0.5,
            'mean_envelope_variance_sigma4': 0.05,
            't_final': 100.0
        }
        self.assertEqual(classify_coarse(summary), 'diffuse coarse field regime')

    def test_diffuse_coarse_field_env_var(self):
        # Condition:
        # mean_basin_count_2 <= 1.05
        # AND NOT basin_dominated_primary
        # AND env_var_4 < 0.02
        summary = {
            'mean_basin_count_sigma2': 1.0,
            'mean_dominant_area_fraction_sigma2': 0.05,
            'mean_dominant_area_fraction_sigma4': 0.1,
            'basin_lifetime_sigma4': 10.0,
            'coarse_persistence_sigma4': 0.5,
            'mean_envelope_variance_sigma4': 0.01,
            't_final': 100.0
        }
        self.assertEqual(classify_coarse(summary), 'diffuse coarse field regime')

    def test_basin_dominated_fallback(self):
        # Condition:
        # mean_basin_count_2 <= 1.05
        # AND NOT basin_dominated_primary
        # AND NOT diffuse_coarse_field
        summary = {
            'mean_basin_count_sigma2': 1.0,
            'mean_dominant_area_fraction_sigma2': 0.1,
            'mean_dominant_area_fraction_sigma4': 0.1,
            'basin_lifetime_sigma4': 10.0,
            'coarse_persistence_sigma4': 0.5,
            'mean_envelope_variance_sigma4': 0.05,
            't_final': 100.0
        }
        self.assertEqual(classify_coarse(summary), 'basin-dominated regime')

    def test_multi_basin_edge(self):
        # Edge case: mean_basin_count_2 = 1.05 (exactly)
        # Should NOT be multi-basin because it's strictly > 1.05
        summary = {
            'mean_basin_count_sigma2': 1.05,
            'mean_dominant_area_fraction_sigma2': 0.1,
            'mean_dominant_area_fraction_sigma4': 0.4,
            'basin_lifetime_sigma4': 10.0,
            'coarse_persistence_sigma4': 1.0,
            'mean_envelope_variance_sigma4': 0.1,
            't_final': 100.0
        }
        self.assertEqual(classify_coarse(summary), 'diffuse coarse field regime')

    def test_basin_dominated_edge(self):
        # Edge case: exactly at thresholds for basin-dominated
        summary = {
            'mean_basin_count_sigma2': 1.05,
            'mean_dominant_area_fraction_sigma2': 0.06,
            'mean_dominant_area_fraction_sigma4': 0.35,
            'basin_lifetime_sigma4': 45.0,
            'coarse_persistence_sigma4': 0.9,
            'mean_envelope_variance_sigma4': 0.02,
            't_final': 100.0
        }
        self.assertEqual(classify_coarse(summary), 'basin-dominated regime')

if __name__ == '__main__':
    unittest.main()
