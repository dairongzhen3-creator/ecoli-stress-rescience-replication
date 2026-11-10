"""Exact partition-bound tests: independent brute-force binary completions on toy ledgers."""
import csv
import itertools
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from partition_weighted_bounds import parse_ledger,summarize,run,digest

FIELDS=['cell','date','STRESS','CAM','fate']

def make_file(tmp, data):
    p=Path(tmp)/'toy.csv'
    with p.open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,lineterminator='\n');w.writeheader();w.writerows(data)
    return p

def r(date,pos,lane,cycle,group,fate):
    return {'cell':f'{date}.{pos}.{lane}.{cycle}','date':date,'STRESS':'H2O2_0.65mM', 'CAM':group,'fate':fate}

class PositionSharpTests(unittest.TestCase):
 def setUp(self):
    self.t=tempfile.TemporaryDirectory()
    self.addCleanup(self.t.cleanup)
 def test_fractional_equal_position_does_not_pool_counts(self):
    toy=[r('20220101','1','1',1,'CAM_after_burst','Survive'),
         r('20220101','1','1',2,'CAM_after_burst','Unknown')]
    toy += [r('20220101','2','1',i,'CAM_after_burst','No Growth') for i in range(1,7)]
    toy += [r('20220101','3','1',1,'CAM_throughout','No Growth'), r('20220101','4','1',1,'Unperturbed','Survive'),r('20220101','5','1',1,'rpoS_KO','No Growth')]
    f=make_file(self.t.name,toy)
    h,pos,g,n=parse_ledger(f,digest(f));per,groups,c,check=summarize(pos)
    group=next(x for x in groups if x['CAM']=='CAM_after_burst')
    self.assertEqual(n,11)
    self.assertEqual(group['lower']['numerator'],1)
    self.assertEqual(group['lower']['denominator'],4) # (1/2+0)/2
    self.assertEqual(group['upper']['decimal'],.5)   # (2/2+0)/2
    self.assertNotEqual(Fraction(1,8),Fraction(1,4))
 def test_global_extrema_bruteforce_for_two_arms(self):
    toy=[r('20220101','1','1',1,'CAM_after_burst','Survive'),
         r('20220101','1','1',2,'CAM_after_burst','Unknown'),
         r('20220101','2','1',1,'CAM_after_burst','Unknown'),
         r('20220101','3','1',1,'CAM_throughout','No Growth'),
         r('20220101','3','1',2,'CAM_throughout','Unknown')]
    f=make_file(self.t.name,toy)
    _,pos,_,n=parse_ledger(f,digest(f))
    # pad missing groups? summarize requires all 4, use extra known cases
    toy +=[r('20220101','4','1',1,'Unperturbed','Survive'),r('20220101','5','1',1,'rpoS_KO','No Growth')]
    f=make_file(self.t.name,toy)
    _,pos,_,_=parse_ledger(f,digest(f))
    pp,groups,contrasts,extra=summarize(pos)
    target=next(x for x in contrasts if x['arm_A']=='CAM_after_burst' and x['arm_B']=='CAM_throughout')
    diffs=[]
    for completion in itertools.product(('Survive','No Growth'),repeat=3):
        tt=[];j=0
        for row in toy:
            row=dict(row)
            if row['fate']=='Unknown':row['fate']=completion[j];j+=1
            tt.append(row)
        p=make_file(self.t.name,tt)
        _,poss,_,_=parse_ledger(p,digest(p))
        _,gg,_,_=summarize(poss)
        x=next(z for z in gg if z['CAM']=='CAM_after_burst')
        y=next(z for z in gg if z['CAM']=='CAM_throughout')
        diffs.append(Fraction(x['lower']['numerator'],x['lower']['denominator'])-Fraction(y['lower']['numerator'],y['lower']['denominator']))
    self.assertEqual(min(diffs),Fraction(target['difference_lower']['numerator'],target['difference_lower']['denominator']))
    self.assertEqual(max(diffs),Fraction(target['difference_upper']['numerator'],target['difference_upper']['denominator']))
 def test_duplicate_id_rejected(self):
    toy=[r('20220101','1','1',1,'CAM_after_burst','Survive')]*2
    f=make_file(self.t.name,toy)
    with self.assertRaisesRegex(ValueError,'DUPLICATE'):parse_ledger(f,digest(f))
 def test_bad_id_and_date_rejected(self):
    toy=[r('20220101','1','1',1,'CAM_after_burst','Unknown')]
    toy[0]['cell']='not_a_software_id';f=make_file(self.t.name,toy)
    with self.assertRaisesRegex(ValueError,'UNRECOGNIZED'):parse_ledger(f,digest(f))
    toy[0]['cell']='20220101.1.1.1';toy[0]['date']='20220202';f=make_file(self.t.name,toy)
    with self.assertRaisesRegex(ValueError,'DATE_KEY_DISAGREEMENT'):parse_ledger(f,digest(f))
 def test_source_sha_and_existing_output_rejected(self):
    toy=[r('20220101','1','1',1,'CAM_after_burst','Survive')]
    f=make_file(self.t.name,toy)
    with self.assertRaisesRegex(ValueError,'SOURCE_SHA_MISMATCH'):parse_ledger(f,'0000')
    dest=Path(self.t.name)/'exists';dest.mkdir()
    with self.assertRaisesRegex(FileExistsError,'OUTPUT_EXISTS'):run(f,dest,expected_hash=digest(f),require_reference=False)
 def test_arm_mixing_in_same_source_position_rejected(self):
    toy=[r('20220101','1','1',1,'CAM_after_burst','Survive'),r('20220101','1','1',2,'CAM_throughout','No Growth')]
    f=make_file(self.t.name,toy)
    _,pos,_,_=parse_ledger(f,digest(f))
    with self.assertRaisesRegex(ValueError,'MULTIPLE_ARM'):summarize(pos)

if __name__=='__main__':unittest.main()
