"""Tests for graph scaling, complete libraries and unavailable measurements."""
from pathlib import Path
import sys
import unittest
import numpy as np
import pandas as pd
from scipy import sparse

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from counts import build_rook_edges, summarise_gene, transform_counts


class ScientificInvariants(unittest.TestCase):
    def test_binary_rook_moran_matches_checkerboard_closed_form(self):
        rows=np.array([0,0,1,1])
        cols=np.array([0,1,0,1])
        left,right=build_rook_edges(rows,cols)
        self.assertEqual(len(left),4)
        result=summarise_gene(np.array([1.,-1.,-1.,1.]),np.array(['Tissue']*4),left,right)
        self.assertAlmostEqual(result['Moran_I'],-1.)
        flat=summarise_gene(np.ones(4),np.array(['Tissue']*4),left,right)
        self.assertTrue(np.isnan(flat['Moran_I']))

    def test_normalization_uses_external_complete_library_not_selected_counts(self):
        counts=sparse.csr_matrix([[1,0],[0,0]])
        transformed=transform_counts(counts,np.array([1000.,0.])).toarray()
        self.assertAlmostEqual(transformed[0,0],np.log2(1001.))
        self.assertTrue(np.all(transformed[1]==0))

    def test_missing_hmgb1_remains_missing_and_tumour_is_test_unit(self):
        data=pd.read_csv(ROOT/'results/Convergent9_three_tumours.csv')
        hmgb=data.loc[data.Gene.eq('HMGB1')]
        self.assertEqual(len(hmgb),3)
        self.assertTrue(hmgb[['External_tumor_minus_tissue','External_Moran_I','External_periphery_minus_interior']].isna().all().all())
        tests=pd.read_csv(ROOT/'results/patient_level_tests.csv')
        self.assertEqual(len(tests),12)
        self.assertTrue(tests.N_samples.eq(3).all())
        self.assertTrue(tests.Exact_one_sided_P.eq(.125).all())


if __name__=='__main__':unittest.main()
