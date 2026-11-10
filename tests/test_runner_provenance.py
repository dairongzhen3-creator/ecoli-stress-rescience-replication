"""Failure-state tests. Synthetic corrupt inputs are never scientific data."""
import json, sys, tempfile, unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import reviewer_run

class ProvenanceStateTests(unittest.TestCase):
    def test_corrupt_original_archive_does_not_claim_verification(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td); fake=p/'not_original.zip';fake.write_bytes(b'not a valid original archive')
            a=type('Args',(),{'author_zip':str(fake),'stress_feather':None,'rmd':None,'realtrace_csv':None,'realtrace_receipt':None,'out':str(p/'output')})()
            with self.assertRaises(RuntimeError):reviewer_run.run(a)
            v=json.loads((p/'output'/'ANALYSIS_RUN_RECEIPT.json').read_text())
            self.assertFalse(v['success'])
            self.assertFalse(v['full_original_zip_reverified_in_this_run'])
            self.assertFalse(v['original_realtrace_member_reverified_in_this_run'])
            self.assertEqual(v['mode'],'ORIGINAL_ZENODO_ZIP')
    def test_missing_preextracted_source_does_not_claim_zip_verified(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)
            a=type('Args',(),{'author_zip':None,'stress_feather':str(p/'missing.feather'),
              'rmd':str(p/'rmd'),'realtrace_csv':str(p/'selected.csv.gz'),
              'realtrace_receipt':str(p/'receipt.json'),'out':str(p/'output')})()
            with self.assertRaises(RuntimeError):reviewer_run.run(a)
            v=json.loads((p/'output'/'ANALYSIS_RUN_RECEIPT.json').read_text())
            self.assertFalse(v['success'])
            self.assertFalse(v['full_original_zip_reverified_in_this_run'])
            self.assertFalse(v['original_realtrace_member_reverified_in_this_run'])
