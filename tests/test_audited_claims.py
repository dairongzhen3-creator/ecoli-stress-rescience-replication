"""Independent table and article arithmetic checks; software-record data only."""
import csv,math,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1]
class AuditedClaims(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  with (R/'references/fig6b_reference_counts.csv').open(newline='') as f:cls.b=list(csv.DictReader(f))
  with (R/'references/fig6d_group_means.csv').open(newline='') as f:cls.d=list(csv.DictReader(f))
 def test_four_condition_total_1153(self):
  self.assertEqual(sum(int(r['Survive'])+int(r['No Growth'])+int(r['Unknown']) for r in self.b),1153)
 def test_group_rates_with_decimal_rounding(self):
  expectations={'CAM_after_burst':.646706586826, 'CAM_throughout':.192660550459,'Unperturbed':.982062780269,'rpoS_KO':.02}
  for row in self.b:
   self.assertAlmostEqual(int(row['Survive'])/(int(row['Survive'])+int(row['No Growth'])),expectations[row['CAM']],places=11)
 def test_unknown_error_denominator_not_known(self):
  rows=[r for r in self.b if r['CAM'].startswith('CAM_')]
  for r in rows:
   s,n,u=map(int,[r['Survive'],r['No Growth'],r['Unknown']]);p=s/(s+n)
   self.assertLess(math.sqrt(p*(1-p)/(s+n+u)),math.sqrt(p*(1-p)/(s+n)))
 def test_eight_group_realtrace_counts_and_means(self):
  self.assertEqual(len(self.d),8)
  self.assertEqual(sum(int(r['N_corr_nonmissing']) for r in self.d),32400)
  for r in self.d:self.assertTrue(math.isfinite(float(r['mean_C'])))
 def test_cam_bola_and_med3_contrasts(self):
  lookup={(r['prom'],r['CAM']):float(r['mean_C']) for r in self.d}
  self.assertGreater(lookup['bolA','CAM_after_burst'],lookup['bolA','CAM_throughout'])
  self.assertGreater(lookup['E08','CAM_after_burst'],lookup['E08','CAM_throughout'])
  self.assertGreater(lookup['bolA','CAM_after_burst'],lookup['bolA','Unperturbed'])
  self.assertLess(lookup['E08','CAM_after_burst'],lookup['E08','Unperturbed'])
if __name__=='__main__':unittest.main()
