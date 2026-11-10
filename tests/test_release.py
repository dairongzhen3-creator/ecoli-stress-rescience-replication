"""Offline tests using published aggregates, NOT biological experiments."""
import csv
import hashlib
import math
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import reviewer_run
import original_zenodo_archive
import realtrace_from_original

class TestOfflineReference(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with (ROOT/'references/fig6b_reference_counts.csv').open(newline='') as f:
            cls.b=list(csv.DictReader(f))
        with (ROOT/'references/fig6c_time_series.csv').open(newline='') as f:
            cls.c=list(csv.DictReader(f))
        with (ROOT/'references/fig6d_group_means.csv').open(newline='') as f:
            cls.d=list(csv.DictReader(f))
    def test_has_expected_source_files(self):
        self.assertEqual(len(self.b),4)
        self.assertEqual(len(self.c),616)
        self.assertEqual(len(self.d),8)
    def test_fig6b_source_counts(self):
        self.assertEqual(sum(int(r['Survive'])+int(r['No Growth'])+int(r['Unknown']) for r in self.b),1153)
    def test_unknown_excluded_from_b_fraction(self):
        for r in self.b:
            s,n,u=(int(r['Survive']),int(r['No Growth']),int(r['Unknown']))
            self.assertAlmostEqual(float(r['fraction_survival']),s/(s+n),places=12)
            self.assertAlmostEqual(float(r['author_error_bar_sqrt_pq_div_N_all']),math.sqrt(s/(s+n)*(1-s/(s+n))/(s+n+u)),places=12)
    def test_rpos_ratio(self):
        r=next(r for r in self.b if r['CAM']=='rpoS_KO')
        self.assertEqual((int(r['Survive']),int(r['No Growth']),int(r['Unknown'])),(6,294,5))
    def test_6c_points_per_condition(self):
        from collections import Counter
        self.assertEqual(Counter(x['CAM'] for x in self.c),{'CAM_after_burst':154,'CAM_throughout':154,'Unperturbed':154,'rpoS_KO':154})
    def test_6d_valid_count(self):
        self.assertEqual(sum(int(r['N_corr_nonmissing']) for r in self.d),32400)
    def test_6d_all_finite(self):
        for r in self.d:
            self.assertGreater(float(r['mean_C']),0)
            self.assertGreaterEqual(float(r['sd_C']),0)
    def test_original_archive_locked(self):
        self.assertEqual(original_zenodo_archive.SOURCE_SIZE,842127932)
        self.assertEqual(original_zenodo_archive.SOURCE_SHA,'0697803959ffa64b8b2bc849d716bc841bcb5afe8a42e92b53f397f8c7d7494b')
    def test_original_realtrace_locked(self):
        self.assertEqual(realtrace_from_original.EXPECTED_MEMBER_SIZE,264629352)
        self.assertEqual(realtrace_from_original.EXPECTED_MEMBER_SHA,'17bdb0bc9e45f0aff5435f0a168d87850a3a29370ccf53449201fbb63f18b217')
    def test_runner_source_reference(self):
        self.assertEqual(reviewer_run.EXPECTED_STRESS,'52c9cf4c7607ecd32a621bb410c967ff5a198d3cdf95f37323008319bb9272ac')
    def test_no_original_author_data_in_bundle(self):
        files=list(ROOT.rglob('*'))
        for f in files:
            if f.is_file():
                self.assertNotEqual(f.suffix.lower(),'.feather')
                self.assertFalse(f.name.endswith('_data.zip'))
                self.assertNotEqual(f.name,'Gervais_et_al_analysis.Rmd')

if __name__=='__main__':unittest.main()
