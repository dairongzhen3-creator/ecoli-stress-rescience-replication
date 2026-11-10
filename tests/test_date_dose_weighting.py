import csv,gzip,sys,tempfile,unittest,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from date_dose_audit import analyze
from reviewer_figures import load_groups,load_strata

class DateDoseWeightingTests(unittest.TestCase):
    def test_date_dose_small_with_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'s.csv.gz'
            with gzip.open(p,'wt',newline='') as f:
                w=csv.DictWriter(f,fieldnames=['CAM','prom','date','STRESS','Corr_selected','cell']);w.writeheader()
                w.writerows([
                dict(CAM='T',prom='F',date='20221212',STRESS='dose1',Corr_selected='10',cell='a'),
                dict(CAM='T',prom='F',date='20221212',STRESS='dose1',Corr_selected='20',cell='a'),
                dict(CAM='T',prom='F',date='20221212',STRESS='dose1',Corr_selected='',cell='b'),
                dict(CAM='T',prom='F',date='20221212',STRESS='dose2',Corr_selected='50',cell='c')])
            out=Path(tmp)/'out';r=analyze(p,out)
            self.assertEqual(r['source_rows'],4);self.assertEqual(r['valid_corrected_rows'],3);self.assertEqual(r['date_dose_strata'],2)
            with (out/'date_dose_equal_stratum_sensitivity.csv').open() as f:v=list(csv.DictReader(f))[0]
            self.assertAlmostEqual(float(v['source_row_weighted_mean_FP_um']),80/3)
            self.assertAlmostEqual(float(v['equal_date_dose_mean_FP_um']),32.5)
            self.assertEqual(int(v['source_software_cells']),3)
            with self.assertRaises(FileExistsError):analyze(p,out)
    def test_date_dose_bad_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'in.gz'
            with gzip.open(p,'wt') as f:f.write('foo,bar\na,b\n')
            with self.assertRaises(ValueError):analyze(p,Path(tmp)/'out')
    def test_date_dose_nonfinite_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'in.gz'
            with gzip.open(p,'wt') as f:
                f.write('CAM,prom,date,STRESS,Corr_selected,cell\nT,F,20220101,D,inf,x\n')
            with self.assertRaises(ValueError):analyze(p,Path(tmp)/'out')
    def test_groups_and_strata_references(self):
        pp=ROOT/'references'
        d=load_groups(pp/'fig6b_group_counts.csv');s=load_strata(pp/'date_dose_binary_completion_bounds.csv')
        self.assertEqual(sum(int(x['Unknown']) for x in d.values()),245)
        self.assertEqual(len(s),8)
        self.assertEqual(sum(x['order_identified']=='False' for x in s),2)
    def test_source_figure_b_sum_rejections(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'figb.csv'
            with p.open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=['group','Survive','No_Growth','Unknown','N_all']);w.writeheader()
                for name in ['Unperturbed','CAM_after_burst','CAM_throughout','rpoS_KO']:
                    w.writerow(dict(group=name,Survive=2,No_Growth=2,Unknown=0,N_all=7))
            with self.assertRaises(AssertionError):load_groups(p)
if __name__=='__main__':unittest.main()
