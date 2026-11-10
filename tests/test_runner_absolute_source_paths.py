import unittest, ast, pathlib
class SourceInputPathRegression(unittest.TestCase):
 def test_input_modes_resolve_path_to_absolute_before_changing_cwd(self):
  src=(pathlib.Path(__file__).resolve().parents[1]/'src/reviewer_run.py').read_text()
  self.assertIn('stress=Path(a.stress_feather).expanduser().resolve()',src)
  self.assertIn('rmd=Path(a.rmd).expanduser().resolve()',src)
  self.assertIn('extract_original_zip(str(Path(a.author_zip).expanduser().resolve()),private)',src)
  ast.parse(src)
