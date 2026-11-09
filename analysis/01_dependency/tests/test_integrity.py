"""Tests use included research snapshots; no network or account connection required."""
import sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from common import load_csv,bh,number
class AnalysisTests(unittest.TestCase):
    def test_bh_known_result(self):
        np.testing.assert_allclose(bh([.01,.04,.03,.002]),[.02,.04,.04,.008])
    def test_bh_missing_is_preserved(self):
        q=bh([np.nan,.01,.2]);self.assertTrue(np.isnan(q[0]));np.testing.assert_allclose(q[1:],[.02,.2])
    def test_model_counts(self):
        rows=load_csv(ROOT/'data/metadata/cell_line_models_1165.csv')
        self.assertEqual(len(rows),1165);self.assertEqual(len(set(r['ModelID'] for r in rows)),1165)
        self.assertEqual(sum(r['Analysis_group']=='COREAD' for r in rows),63)
        self.assertEqual(sum(r['Analysis_group']=='PRIMARY_COMPARATOR' for r in rows),1099)
    def test_primary_integration(self):
        rows=load_csv(ROOT/'results/csv/integration_results_1165.csv');self.assertEqual(len(rows),1720)
        for r in rows:
            self.assertEqual(r['Union_candidate']=='YES',r['CL_candidate']=='YES' or r['Organoid_3DO_candidate']=='YES')
            self.assertEqual(r['Both_systems_candidate']=='YES',r['CL_candidate']=='YES' and r['Organoid_3DO_candidate']=='YES')
    def test_missing_not_zero_effect(self):
        rows=load_csv(ROOT/'results/csv/organoid_results_1720.csv')
        missing=[r for r in rows if r['Data_status']!='EVALUABLE'];self.assertEqual(len(missing),84)
        self.assertTrue(all(r['Chronos_mean_3DO']=='' for r in missing))
    def test_model_overlap_not_donor_claim(self):
        cl={r['ModelID'] for r in load_csv(ROOT/'data/metadata/cell_line_models_1165.csv')}
        org={r['ModelID'] for r in load_csv(ROOT/'data/metadata/organoid_coread_models_and_eb_coefficients.csv')}
        self.assertFalse(cl&org)
if __name__=='__main__':unittest.main()
