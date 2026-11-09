"""Read-only regression tests of the shipped numerical data. No network is used."""
from pathlib import Path
import csv
import json
import unittest
import numpy as np
from scipy.stats import false_discovery_control

R = Path(__file__).resolve().parents[1]

def read(name):
    with (R/'results'/name).open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))

class ArchivedAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.genes = (R/'input/union_candidates_466.txt').read_text().split()
        cls.primary = read('primary_466_results.csv')
        cls.external = read('external_466_results_wide.csv')
        cls.statistics = read('permutation_statistics.csv')

    def test_frozen_universe_and_cohorts(self):
        self.assertEqual(len(self.genes), 466)
        self.assertEqual(len(set(self.genes)), 466)
        self.assertEqual({r['Gene'] for r in self.primary}, set(self.genes))
        self.assertEqual(len(self.primary), 466)
        self.assertEqual(len(self.external), 1398)
        for co in ['AC_ICAM', 'CPTAC2', 'PANCAN']:
            self.assertEqual(sum(r['Cohort'] == co for r in self.external), 466)

    def test_primary_unique_patient_union(self):
        z = np.load(R/'results/COREAD_GDC_calls.npz', allow_pickle=False)
        joint = z['avail__MUTATIONS'] & z['avail__CNA']
        self.assertEqual(int(joint.sum()), 572)
        self.assertEqual(len(z['patients']), 631)
        event = (z['event__Driver_mutation'] | z['event__Amplification'] | z['event__Deep_deletion']).astype(bool)
        counts = event[:, joint].sum(axis=1)
        rows = {r['Gene']: r for r in self.primary}
        for g, k in zip(z['genes'], counts):
            self.assertEqual(int(rows[g]['DNA_driver_CNA_altered']), int(k))
            self.assertAlmostEqual(float(rows[g]['DNA_driver_CNA_frequency']), k/572, places=13)
        self.assertEqual(rows['APC']['DNA_driver_CNA_altered'], '417')
        self.assertEqual(rows['KRAS']['DNA_driver_CNA_altered'], '244')

    def test_modality_frequencies_and_missingness(self):
        for r in self.external:
            for ep in ['DNA_driver_CNA', 'Driver_mutation', 'Amplification', 'Deep_deletion',
                       'CNA_any', 'Any_exported_mutation', 'DNA_broad', 'RNA_high', 'RNA_low']:
                n = float(r[ep+'_N']) if r[ep+'_N'] else 0
                f = r[ep+'_Frequency']
                if n > 0:
                    self.assertAlmostEqual(float(f), float(r[ep+'_Altered'])/n, places=13)
                else:
                    self.assertEqual(f, '')
        ac = [r for r in self.external if r['Cohort'] == 'AC_ICAM']
        self.assertEqual(sum(float(r['Amplification_N'] or 0)>0 for r in ac), 456)
        self.assertTrue(all(r['Driver_mutation_Frequency']=='' and r['DNA_driver_CNA_Frequency']=='' for r in ac))

    def test_empirical_probabilities_and_null_summaries(self):
        z = np.load(R/'results/permutation_null_draws.npz', allow_pickle=False)
        self.assertEqual(len(self.statistics), 50)
        self.assertEqual(len(z.files), 46)
        unavailable = 0
        for r in self.statistics:
            if not r['Empirical_P_greater']:
                unavailable += 1
                self.assertEqual(r['Permutations'], '0')
                continue
            a = z[r['Test_ID']]
            self.assertEqual(a.shape, (10000,))
            b = int((a >= float(r['Spearman_rho'])-1e-12).sum())
            self.assertEqual(b, int(r['Null_exceedances']))
            self.assertAlmostEqual(float(r['Empirical_P_greater']), (b+1)/10001, places=14)
            self.assertAlmostEqual(float(r['Null_mean']), float(a.mean()), places=13)
            self.assertAlmostEqual(float(r['Null_SD']), float(a.std(ddof=1)), places=13)
            np.testing.assert_allclose(np.quantile(a, [.025,.975]),
                [float(r['Null_percentile_2p5']), float(r['Null_percentile_97p5'])], atol=1e-13)
        self.assertEqual(unavailable, 4)

    def test_BH_families(self):
        families = {'EXTERNAL_COREAD_CONCORDANCE':18,
                    'COLON_ONLY_EXTERNAL_SENSITIVITY':18,
                    'TCGA_REPRESENTATION_ROBUSTNESS':10}
        for family, n in families.items():
            rows = [r for r in self.statistics if r['BH_family']==family and r['Empirical_P_greater']]
            self.assertEqual(len(rows), n)
            expected = false_discovery_control([float(r['Empirical_P_greater']) for r in rows], method='bh')
            np.testing.assert_allclose(expected, [float(r['BH_FDR']) for r in rows], atol=1e-14, rtol=1e-12)

    def test_bootstrap_ranks(self):
        z = np.load(R/'results/primary_rank_bootstrap.npz', allow_pickle=False)
        self.assertEqual(z['rank'].shape, (2000,466))
        self.assertEqual(z['frequencies'].shape, (2000,466))
        med = np.median(z['rank'], axis=0)
        q = np.quantile(z['rank'], [.025,.975], axis=0)
        p = {r['Gene']:r for r in self.primary}
        for i,g in enumerate(z['genes']):
            self.assertAlmostEqual(med[i], float(p[g]['DNA_rank_boot_median']))
            self.assertAlmostEqual(q[0,i], float(p[g]['DNA_rank_boot_CI95_low']))
            self.assertAlmostEqual(q[1,i], float(p[g]['DNA_rank_boot_CI95_high']))

    def test_frozen_top30_and_ACICAM_patient_extract(self):
        top = (R/'results/primary_top30_genes.txt').read_text().split()
        self.assertEqual(top, [r['Gene'] for r in self.primary[:30]])
        z = np.load(R/'results/AC_ICAM_top30_calls.npz', allow_pickle=False)
        self.assertEqual(list(z['genes']), top)
        self.assertEqual(z['codes'].shape, (30,281))
        self.assertEqual(len(set(z['patients'])),281)
        union = ((z['codes'] & 7)>0).sum(axis=1)
        rows = {r['Gene']:r for r in self.external if r['Cohort']=='AC_ICAM'}
        for g,k in zip(top,union):
            self.assertEqual(k,int(rows[g]['DNA_broad_Altered']))

    def test_compact_source_values(self):
        files = [R.parents[1]/'source_data/figure_2'/f'Figure2{letter}_source.csv' for letter in 'ABC']
        self.assertEqual(len(files), 3)
        primary = {r['Gene']:r for r in self.primary}
        external = {(r['Cohort'],r['Gene']):r for r in self.external}
        cats = ['mutation','amp','del','rna_high','rna_low']
        for file in files:
            with file.open(newline='') as f:
                rows = list(csv.DictReader(f))
            self.assertEqual([r['Gene'] for r in rows],[r['Gene'] for r in self.primary[:30]])
            for r in rows:
                marginal = np.array([float(r[c+'_frequency']) for c in cats])
                segments = np.array([float(r[c+'_composition']) for c in cats])
                np.testing.assert_allclose(segments,marginal/marginal.sum(),rtol=1e-13)
                self.assertAlmostEqual(float(r['DNA_frequency']),int(r['DNA_altered'])/int(r['DNA_n']),places=13)
                if r['Cohort']=='TCGA':
                    expected = primary[r['Gene']]['DNA_driver_CNA_frequency']
                else:
                    endpoint = 'DNA_broad_Frequency' if r['Cohort']=='AC_ICAM' else 'DNA_driver_CNA_Frequency'
                    expected = external[(r['Cohort'],r['Gene'])][endpoint]
                self.assertAlmostEqual(float(r['DNA_frequency']),float(expected),places=13)

    def test_no_editorial_graphics_shipped(self):
        self.assertFalse(any((R/'figures').rglob('*.svg')))
        self.assertFalse(any((R/'figures').rglob('*.png')))

if __name__ == '__main__':
    unittest.main()
