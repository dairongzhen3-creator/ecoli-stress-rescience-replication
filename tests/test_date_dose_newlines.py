import pathlib,unittest
class DateDosePortableNewlines(unittest.TestCase):
 def test_frozen_reference_uses_lf_and_runner_writes_lf(self):
  root=pathlib.Path(__file__).resolve().parents[1]
  s=(root/'src/date_dose_audit.py').read_text()
  self.assertIn('lineterminator="\\n"',s)
  ref=(root/'references/fig6d_group_means.csv').read_bytes()
  self.assertNotIn(b'\r\n',ref)
