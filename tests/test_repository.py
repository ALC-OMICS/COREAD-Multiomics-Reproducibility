"""Integrity-policy regression tests; no editing of the live repository."""
from pathlib import Path
import importlib.util,tempfile,unittest
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('repository_verifier',ROOT/'scripts/verify_repository.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
class ScopeTests(unittest.TestCase):
    def test_current_repository_integrity(self):
        report=v.verify(ROOT)
        self.assertTrue(report['passed'],report['errors'][:20])
    def test_path_traversal_rejected(self):
        for name in ['../secret','/absolute','.git/config','x/../../a','a\\b']:
            with self.subTest(name=name),self.assertRaises(ValueError):v.safe_path(ROOT,name)
    def test_editorial_files_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            for name in ['plot.svg','paper.docx','supplement.xlsx','manuscript/draft.md','figures/test.csv']:
                p=root/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('test')
            self.assertEqual(len(v.scope_errors(root)),5)
    def test_numeric_sources_allowed(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/'source_data/figure_1/counts.csv';p.parent.mkdir(parents=True);p.write_text('gene,count\nA,1\n')
            self.assertEqual(v.scope_errors(root),[])
    def test_missing_inputs_not_blessed_by_refresh(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):v.write_manifest(Path(d))
if __name__=='__main__':unittest.main()
