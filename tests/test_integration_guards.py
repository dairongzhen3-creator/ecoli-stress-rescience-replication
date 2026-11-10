"""Independent negative controls for integrated exact-bounds workflow."""
import sys,unittest,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'src'))
from sharp_bounds import run
class IntegrationGuards(unittest.TestCase):
 def test_missing_source_does_not_create_output(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);out=p/'out'
   with self.assertRaises(FileNotFoundError):run(p/'missing1.csv',p/'missing2.csv',out)
   self.assertFalse(out.exists())
 def test_tamper_source_has_hash_gate(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);a=p/'a.csv';b=p/'b.csv';out=p/'out'
   a.write_text('S,N\n1,2\n');b.write_text('x,y\na,b\n')
   with self.assertRaisesRegex(ValueError,'INTEGRITY_FAIL'):run(a,b,out)
   self.assertFalse(out.exists())
 def test_preexisting_output_rejected(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);out=p/'out';out.mkdir()
   with self.assertRaises(FileExistsError):run(p/'a',p/'b',out)
 def test_master_runner_has_missing_outcome_step(self):
  body=(ROOT/'src/reviewer_run.py').read_text()
  self.assertIn('sharp_binary_completion',body)
  self.assertIn('MISSING_OUTCOME_BOUNDS',body)
if __name__=='__main__':unittest.main()
