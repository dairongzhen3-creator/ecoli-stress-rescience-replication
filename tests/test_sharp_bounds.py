"""Independent brute-force verification of sharp unknown-label AUC optimizer."""
import csv
from fractions import Fraction
import itertools
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'src'))
from sharp_bounds import auc_known, sharp_auc_bounds, midranks2


def brute_auc(pairs):
    p=[(score,lab) for score,lab in pairs if lab!='Unknown']
    npos=sum(lab=='Survive' for _,lab in p)
    nneg=sum(lab=='No Growth' for _,lab in p)
    if npos==0 or nneg==0:return None
    win2=0
    for a,b in itertools.product([score for score,lab in p if lab=='Survive'],
                                 [score for score,lab in p if lab=='No Growth']):
        win2+=int(a>b)*2+int(a==b)
    return Fraction(win2,2*npos*nneg)


def brute_bounds(pairs):
    unknown_indices=[i for i,(s,lab) in enumerate(pairs) if lab=='Unknown']
    vals=[]
    for bits in itertools.product(['Survive','No Growth'],repeat=len(unknown_indices)):
        q=list(pairs)
        for i,lab in zip(unknown_indices,bits):q[i]=(q[i][0],lab)
        auc=brute_auc(q)
        if auc is not None:vals.append(auc)
    return min(vals),max(vals)

class TestSharpBounds(unittest.TestCase):
    def test_midrank_ties(self):
        self.assertEqual(midranks2([1.0,2.0,2.0,0.0]),[4,7,7,2])

    def test_exact_known_matches_pairwise(self):
        q=[(0.2,'Survive'),(0.2,'No Growth'),(1.0,'Survive'),(1.5,'No Growth'),(2.0,'Unknown')]
        self.assertEqual(auc_known(q),brute_auc(q))

    def test_unknown_sharp_by_exhaustive_enumeration(self):
        scenarios=[
            [(0,'Survive'),(1,'No Growth'),(2,'Unknown')],
            [(1,'Survive'),(1,'No Growth'),(1,'Unknown'),(2,'Unknown')],
            [(5,'Survive'),(2,'No Growth'),(7,'Unknown'),(3,'Unknown'),(8,'Unknown')],
            [(0,'Survive'),(2,'Survive'),(1,'No Growth'),(3,'No Growth'),(4,'Unknown'),(0,'Unknown'),(2,'Unknown')],
            [(0,'No Growth'),(2,'Unknown'),(1,'Unknown')],
            [(0,'Survive'),(2,'Unknown'),(1,'Unknown')]
        ]
        for pairs in scenarios:
            with self.subTest(pairs=pairs):
                min_b,max_b=brute_bounds(pairs)
                sharp=sharp_auc_bounds(pairs)
                self.assertEqual(Fraction(sharp['sharp_lower']['exact']),min_b)
                self.assertEqual(Fraction(sharp['sharp_upper']['exact']),max_b)
                self.assertGreaterEqual(min_b,0)
                self.assertLessEqual(max_b,1)

    def test_group_robust_order_independent(self):
        p=ROOT/'references/missing_outcome_bounds_receipt.json'
        import json
        o=json.loads(p.read_text())
        self.assertEqual(sum(g['Unknown'] for g in o['group_bounds'].values()),245)
        self.assertEqual(sum(g['Survive']+g['No_Growth']+g['Unknown'] for g in o['group_bounds'].values()),1153)
        self.assertTrue(all(x['ordering_invariant_for_all_binary_unknown_labels'] for x in o['adjacent_group_contrast_certificate']))
        self.assertEqual(Fraction(o['adjacent_group_contrast_certificate'][1]['sharp_min_difference']['exact']),Fraction(5,2378))

    def test_known_only_auc_regression(self):
        import json
        o=json.loads((ROOT/'references/missing_outcome_bounds_receipt.json').read_text())
        # Independent reference values come from the frozen reference table, rounded 5 decimals.
        expected={('CAM_after_burst','bolA'):0.65013,('CAM_after_burst','Med3'):0.59401,
                  ('CAM_throughout','bolA'):0.30628,('CAM_throughout','Med3'):0.71591}
        for (group,reporter),v in expected.items():
            with self.subTest(group=group,reporter=reporter):
                got=o['exploratory_auc_sharp_unknown_bounds'][group][reporter]
                self.assertAlmostEqual(got['known_only']['decimal'],v,places=5)
                self.assertLess(got['sharp_lower']['decimal'],.5)
                self.assertGreater(got['sharp_upper']['decimal'],.5)

if __name__=='__main__':unittest.main()
