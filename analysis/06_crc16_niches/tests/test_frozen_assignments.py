"""Frozen CSV assignments must survive reordering, not corruption or relabelling."""
from pathlib import Path
import importlib.util
import tempfile
import unittest
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('crc16_extractor', ROOT/'extract_crc16.py')
extractor = importlib.util.module_from_spec(spec)
spec.loader.exec_module(extractor)

class FrozenAssignmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frozen_path = ROOT/'input/CRC16_frozen_niches.csv'
        cls.frozen = pd.read_csv(cls.frozen_path, index_col='barcode')
        cls.coordinates = cls.frozen[['tissue','row','col','imagerow','imagecol']].copy()
        cls.sizes = pd.read_csv(ROOT/'provenance/bacterial_burden_by_niche.csv')

    def test_exact_source_map(self):
        result = extractor.recover_niches(self.frozen_path, self.coordinates, self.sizes)
        self.assertEqual(len(result),3106)
        self.assertEqual(result.Niche.nunique(),13)
        pd.testing.assert_series_equal(result.Niche,self.frozen.Niche)

    def test_source_reordering_preserves_barcode_labels(self):
        coordinates = self.coordinates.iloc[::-1]
        result = extractor.recover_niches(self.frozen_path,coordinates,self.sizes)
        self.assertTrue(result.index.equals(coordinates.index))
        pd.testing.assert_series_equal(result.Niche,self.frozen.loc[coordinates.index,'Niche'])

    def test_coordinate_disagreement_rejected(self):
        changed = self.coordinates.copy()
        changed.iloc[0,changed.columns.get_loc('imagecol')] += 1
        with self.assertRaisesRegex(ValueError,'coordinates differ'):
            extractor.recover_niches(self.frozen_path,changed,self.sizes)

    def test_missing_and_duplicate_barcodes_rejected(self):
        for changed in (self.coordinates.iloc[1:],pd.concat([self.coordinates,self.coordinates.iloc[:1]])):
            with self.subTest(length=len(changed)),self.assertRaises(ValueError):
                extractor.recover_niches(self.frozen_path,changed,self.sizes)

    def test_changed_niche_sizes_rejected(self):
        sizes = self.sizes.copy()
        sizes.loc[0,'n_spots'] += 1
        with self.assertRaisesRegex(ValueError,'sizes differ'):
            extractor.recover_niches(self.frozen_path,self.coordinates,sizes)

if __name__=='__main__':
    unittest.main()
