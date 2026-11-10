"""Finite-record tipping-point tests, including negative controls."""
import csv,sys,tempfile,unittest
from copy import deepcopy
from fractions import Fraction
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from tipping_point import analyze
R=Path(__file__).resolve().parents[1]
SRC=R/'references/source_id_partition_counts.csv'
class TippingPointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.actual=analyze(SRC)
    def altered(self,change):
        with SRC.open(encoding='utf8',newline='') as f:rd=csv.DictReader(f);rows=list(rd);fields=rd.fieldnames
        change(rows)
        with tempfile.TemporaryDirectory() as td:
            fn=Path(td)/'input.csv'
            with fn.open('w',encoding='utf8',newline='') as o:
                w=csv.DictWriter(o,fieldnames=fields);w.writeheader();w.writerows(rows)
            return analyze(fn)
    def test_partition_count(self):self.assertEqual(self.actual['source_id_partition_count'],37)
    def test_unknown_counts(self):self.assertEqual(self.actual['relevant_unresolved_labels'],{'Unperturbed':65,'CAM_after_burst':65})
    def test_exact_negative_endpoint(self):self.assertEqual(Fraction(self.actual['sharp_minimum_difference_fraction']),Fraction(-2539,3533145))
    def test_minimum_one_label_change(self):self.assertEqual(Fraction(self.actual['smallest_single_label_shift_fraction']),Fraction(1,540))
    def test_one_flip_removes_reversal(self):self.assertGreater(self.actual['minimum_after_one_opposite_extreme_flipped_pp'],0)
    def test_common_missing_propensity_positive(self):
        self.assertGreater(self.actual['illustrative_shared_unknown_regrowth_fraction_q0_pp'],0)
        self.assertGreater(self.actual['illustrative_shared_unknown_regrowth_fraction_q1_pp'],0)
    def test_rejects_invalid_counts(self):
        with self.assertRaisesRegex(ValueError,'Invalid source counts'):
            self.altered(lambda rows:rows[0].update({'S':'-1'}))
    def test_rejects_duplicated_source_partition(self):
        with self.assertRaisesRegex(ValueError,'Duplicated date'):
            self.altered(lambda rows:rows[1].update({'date':rows[0]['date'],'position_from_software_id':rows[0]['position_from_software_id']}))
