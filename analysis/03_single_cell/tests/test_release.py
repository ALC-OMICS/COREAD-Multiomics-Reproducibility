import csv, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

class ReleaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with (ROOT/'results/singlecell_466_results_with_sensitivities.csv').open(encoding='utf-8',newline='') as f:
            cls.rows=list(csv.DictReader(f))
    def test_frozen_universe(self):
        self.assertEqual(len(self.rows),466); self.assertEqual(len({r['Gene'] for r in self.rows}),466)
    def test_tiers(self):
        self.assertEqual(sum(r['SingleCell_priority_tier']=='A' for r in self.rows),65)
        self.assertEqual(sum(r['SingleCell_priority_tier']=='B' for r in self.rows),134)
        self.assertEqual(sum(r['SingleCell_priority_tier']=='C' for r in self.rows),267)
    def test_core(self):
        core=[r for r in self.rows if r['Spatial_followup_core_47']=='YES']
        self.assertEqual(len(core),47)
        self.assertTrue(all(r['SingleCell_priority_tier']=='A' for r in core))
        self.assertTrue(all(r['Composition_balanced_positive_FDR05']=='YES' for r in core))
    def test_gene_matching(self):
        self.assertEqual(sum(r['SC_match_method']=='Exact_symbol' for r in self.rows),455)
        self.assertEqual(sum(r['SC_match_method']!='Exact_symbol' for r in self.rows),11)
    def test_primary_counts(self):
        sig=[r for r in self.rows if float(r['Tumor_vs_normal_epi_BH_FDR'])<.05]
        self.assertEqual(len(sig),334)
        self.assertEqual(sum(float(r['Tumor_vs_normal_epi_median_delta_log2CPM'])>0 for r in sig),219)
        self.assertEqual(sum(float(r['Tumor_vs_normal_epi_median_delta_log2CPM'])<0 for r in sig),112)
    def test_mmr_counts(self):
        sig=[r for r in self.rows if float(r['MMRd_vs_MMRp_epi_BH_FDR'])<.05]
        self.assertEqual(len(sig),139)
        self.assertEqual(sum(float(r['MMRd_vs_MMRp_epi_median_delta_log2CPM'])>0 for r in sig),96)
        self.assertEqual(sum(float(r['MMRd_vs_MMRp_epi_median_delta_log2CPM'])<0 for r in sig),40)

if __name__=='__main__': unittest.main()
