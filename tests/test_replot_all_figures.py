import unittest, tempfile, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'src'))
from replot_all_figures import SOURCE_NAMES,load_csv
class ReviewPlotTests(unittest.TestCase):
 def test_all_source_inputs_included(self):
  root=pathlib.Path(__file__).resolve().parents[1]/'references'
  self.assertTrue(all((root/file).is_file() for file in SOURCE_NAMES.values()))
 def test_computational_sources_have_expected_record_counts(self):
  root=pathlib.Path(__file__).resolve().parents[1]/'references'
  self.assertEqual(len(load_csv(root/SOURCE_NAMES['fig6c'])),616)
  self.assertEqual(len(load_csv(root/SOURCE_NAMES['fig6d'])),8)
 def test_reject_existing_output(self):
  from replot_all_figures import plot_all
  root=pathlib.Path(__file__).resolve().parents[1]/'references'
  with tempfile.TemporaryDirectory() as tmp:
   with self.assertRaises(FileExistsError):plot_all(root,pathlib.Path(tmp))
