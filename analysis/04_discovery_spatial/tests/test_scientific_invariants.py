"""Regression tests for gene-key correspondence and biological aggregation."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('discovery_run',ROOT/'run.py')
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ScientificInvariants(unittest.TestCase):
    def test_moran_gene_join_is_invariant_to_input_row_order(self):
        source=pd.read_csv(ROOT/'data/moran_section_466.csv')
        _,ordered=module.keyed_moran(source)
        _,shuffled=module.keyed_moran(source.sample(frac=1,random_state=20261031))
        pd.testing.assert_frame_equal(ordered.sort_values('Gene').reset_index(drop=True),
                                      shuffled.sort_values('Gene').reset_index(drop=True),
                                      check_exact=False,rtol=1e-12,atol=1e-14)
        self.assertEqual(int(ordered.Moran_Wilcoxon_P.notna().sum()),441)

    def test_serial_sections_are_averaged_before_patient_testing(self):
        section=pd.read_csv(ROOT/'data/section_gene_effects.csv.gz')
        reference=pd.read_csv(ROOT/'data/patient_gene_effects.csv')
        metrics=['interface_effect','tumor_vs_stroma_effect','distance_rho']
        regenerated=section.groupby(['patient','Gene'],as_index=False)[metrics].mean()
        report=module.compare_frames(regenerated,reference[['patient','Gene']+metrics],('patient','Gene'))
        self.assertTrue(report['passed'],report)

    def test_current_hierarchy_is_nested_and_c9_members_are_explicit(self):
        data=pd.read_csv(ROOT/'results/discovery_spatial_recomputed.csv')
        sets={c:set(data.loc[data[c],'Gene']) for c in ['Core47','Spatially_substantiated_Core47','Extended13_5pct','Convergent9']}
        self.assertEqual([len(sets[c]) for c in sets],[47,43,13,9])
        self.assertLessEqual(sets['Convergent9'],sets['Extended13_5pct'])
        self.assertLessEqual(sets['Extended13_5pct'],sets['Spatially_substantiated_Core47'])
        self.assertLessEqual(sets['Spatially_substantiated_Core47'],sets['Core47'])
        self.assertEqual(sets['Convergent9'],{'BCL2L1','CTNNB1','FASN','HMGB1','IRS2','MYC','RNF6','SOX9','YAP1'})


if __name__=='__main__':unittest.main()
