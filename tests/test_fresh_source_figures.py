"""Regression gate: source-to-figure rendering must never substitute frozen reference results.

Small synthetic fixtures reuse already certified derived output tables, not biological raw data.
"""
import csv
import json
import shutil
import tempfile
import unittest
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from replot_all_figures import plot_all,computed_source_files,COMPUTED_NAMES,SOURCE_NAMES

# Test fixtures are bundled separately; the unit test needs no prior private project directory.
FIXTURE = ROOT/'tests'/'fixtures'/'seven_stage_computed_summary'

class FreshFigureSourceTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='fresh_figures_')
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.computed=self.base/'computed'
        self.computed.mkdir()
        for key,rel in COMPUTED_NAMES.items():
            dest=self.computed/rel
            dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(ROOT/'references'/SOURCE_NAMES[key],dest)
        shutil.copy2(FIXTURE/'ANALYSIS_RUN_RECEIPT.json',self.computed/'ANALYSIS_RUN_RECEIPT.json')

    def test_source_run_produces_identical_unchanged_figures_to_staged_reference(self):
        produced=self.base/'fresh_images'
        frozen=self.base/'frozen_images'
        a=plot_all(self.computed,produced,computed_run=True)
        b=plot_all(ROOT/'references',frozen)
        self.assertEqual(len(a),6)
        self.assertEqual(len(b),6)
        self.assertTrue((produced/'FIGURE_SOURCE_SHA256_RECEIPT.json').is_file())
        receipt=json.loads((produced/'FIGURE_SOURCE_SHA256_RECEIPT.json').read_text())
        self.assertEqual(set(receipt['source_inputs_sha256']),set(COMPUTED_NAMES))
        for p in sorted(produced.glob('*.png')):
            self.assertEqual(p.read_bytes(),(frozen/p.name).read_bytes(),p.name)

    def test_source_run_uses_computed_data_not_stale_repository_reference(self):
        baseline=self.base/'before'
        changed=self.base/'after'
        plot_all(self.computed,baseline,computed_run=True)
        groups=self.computed/COMPUTED_NAMES['groups']
        with groups.open(encoding='utf8',newline='') as f:
            reader=csv.DictReader(f)
            data=list(reader);fields=reader.fieldnames
        d=next(r for r in data if r['group']=='CAM_after_burst')
        # Alter ONLY a synthetic, temporary copy. Maintain S+N+U counts.
        d['Survive']=str(int(d['Survive'])+1)
        d['No_Growth']=str(int(d['No_Growth'])-1)
        with groups.open('w',encoding='utf8',newline='') as f:
            wr=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');wr.writeheader();wr.writerows(data)
        plot_all(self.computed,changed,computed_run=True)
        self.assertNotEqual((baseline/'figure1_recovery_fraction.png').read_bytes(),
                            (changed/'figure1_recovery_fraction.png').read_bytes())
        for name in ('figure2_bola_time_series.png','figure3_bola_reporter_means.png',
                     'figure4_med3_reporter_means.png','figure6_reporter_rank_auc.png',
                     'figure5_matched_date_dose_bounds.png'):
            self.assertEqual((baseline/name).read_bytes(),(changed/name).read_bytes(),name)

    def test_reject_failed_run_receipt(self):
        verdict=self.computed/'ANALYSIS_RUN_RECEIPT.json'
        data=json.loads(verdict.read_text(encoding='utf8'))
        data['success']=False
        verdict.write_text(json.dumps(data),encoding='utf8')
        with self.assertRaisesRegex(ValueError,'has not passed'):
            computed_source_files(self.computed)

    def test_reject_incomplete_run_receipt(self):
        verdict=self.computed/'ANALYSIS_RUN_RECEIPT.json'
        data=json.loads(verdict.read_text(encoding='utf8'))
        data['steps']=data['steps'][:-1]
        verdict.write_text(json.dumps(data),encoding='utf8')
        with self.assertRaisesRegex(ValueError,'has not passed'):
            computed_source_files(self.computed)

    def test_reject_incorrect_source_hash_in_verdict(self):
        verdict=self.computed/'ANALYSIS_RUN_RECEIPT.json'
        data=json.loads(verdict.read_text(encoding='utf8'))
        data['source_sha256']['rmd']='0'*64
        verdict.write_text(json.dumps(data),encoding='utf8')
        with self.assertRaisesRegex(ValueError,'source hashes'):
            computed_source_files(self.computed)

    def test_reject_missing_run_source(self):
        (self.computed/COMPUTED_NAMES['auc']).unlink()
        with self.assertRaises(FileNotFoundError):
            computed_source_files(self.computed)

if __name__=='__main__':unittest.main()
