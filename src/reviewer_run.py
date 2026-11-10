#!/usr/bin/env python3
"""Run all three Gervais 2026 Fig6B/C/D source computations, without copying original data into code package.

This is a numerical software replication. Corrected RealTrace table is derived
from *author* observations, not a new physical experiment. Validation of the
physical cell, exposure time and memory carrier is out of scope.
"""
import argparse,hashlib,json,os,platform,subprocess,sys,traceback
from pathlib import Path
from original_zenodo_archive import extract as extract_original_zip
from realtrace_from_original import extract as extract_realtrace
from fig6bcd_crosspanel import run as crosspanel
SRC=Path(__file__).resolve().parent
REPO=SRC.parent
EXPECTED_STRESS='52c9cf4c7607ecd32a621bb410c967ff5a198d3cdf95f37323008319bb9272ac'
EXPECTED_RT_ZIP='0697803959ffa64b8b2bc849d716bc841bcb5afe8a42e92b53f397f8c7d7494b'
EXPECTED_RT_MEMBER='17bdb0bc9e45f0aff5435f0a168d87850a3a29370ccf53449201fbb63f18b217'

def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for part in iter(lambda:f.read(4<<20),b''):h.update(part)
 return h.hexdigest()

def run(a):
 out=Path(a.out).expanduser().resolve()
 if out.exists():raise FileExistsError('OUTPUT_EXISTS: refusing to overwrite '+str(out))
 if bool(a.author_zip)==bool(a.stress_feather):raise ValueError('Choose either --author-zip OR --stress-feather')
 if a.stress_feather and not(a.rmd and a.realtrace_csv and a.realtrace_receipt):
  raise ValueError('Pre-extracted mode needs --rmd, --realtrace-csv and --realtrace-receipt')
 out.mkdir(parents=True)
 logs=out/'logs';logs.mkdir()
 state={'mode':'ORIGINAL_ZENODO_ZIP' if a.author_zip else 'PRIOR_VERIFIED_COLAB_SELECTED_SOURCE',
        'python':sys.version,'platform':platform.platform(),'steps':[],
        'full_original_zip_reverified_in_this_run':False,
        'original_realtrace_member_reverified_in_this_run':False,
        'source_defined_software_replication_only':True,'physical_biological_replication':False}
 try:
  if a.author_zip:
   private=out/'AUTHOR_INPUT_PRIVATE';ex=extract_original_zip(str(Path(a.author_zip).expanduser().resolve()),private)
   # Set source-verification flags only after the authentic source/member checks complete.
   state['full_original_zip_reverified_in_this_run']=True
   state['original_realtrace_member_reverified_in_this_run']=True
   stress=Path(ex['stress']['file']);rmd=Path(ex['rmd']['file'])
   rt=out/'REALTRACE_DERIVED_PRIVATE.csv.gz';rec=out/'REALTRACE_DERIVED_RECEIPT.json'
   extract_realtrace(ex['realtrace']['file'],rt,rec,ex['source_zip_sha256'],ex['realtrace']['crc32'])
  else:
   stress=Path(a.stress_feather).expanduser().resolve();rmd=Path(a.rmd).expanduser().resolve();rt=Path(a.realtrace_csv).expanduser().resolve();rec=Path(a.realtrace_receipt).expanduser().resolve()
  if sha(stress)!=EXPECTED_STRESS:raise RuntimeError('Stress Feather SHA mismatch')
  state['source_sha256']={'stress':sha(stress),'rmd':sha(rmd),'selected_realtrace':sha(rt)}
  def step(name,script,arguments):
   cmd=[sys.executable,str(SRC/script)]+[str(v) for v in arguments]
   proc=subprocess.run(cmd,cwd=str(REPO),text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
   (logs/(name+'.log')).write_text(proc.stdout,encoding='utf8')
   state['steps'].append({'name':name,'returncode':proc.returncode,'cmd':[Path(x).name if i==1 else x for i,x in enumerate(cmd)]})
   if proc.returncode:raise RuntimeError(f'{name} failed ({proc.returncode}): '+proc.stdout[-1400:])
  cycle_dir=out/'FIG6B';reporter_source_dir=out/'CANDIDATE_SOURCE';reporter_summary_dir=out/'FIG6CD';crosspanel_dir=out/'PANEL_CROSSWALK'
  step('Fig6B','fig6b.py',['--feather',stress,'--rmd',rmd,'--reference',REPO/'references/fig6b_reference_counts.csv','--out',cycle_dir])
  partition_dir=out/'PARTITION_WEIGHT_BOUNDS'
  step('equal_partition_binary_completion','partition_weighted_bounds.py',[
       '--source',cycle_dir/'fig6b_cycle_ledger.csv', '--out',partition_dir])
  partition_certificate=json.loads((partition_dir/'partition_weighted_bounds_receipt.json').read_text())
  if partition_certificate['date_position_units']!=37:raise RuntimeError('source ID partition count drift')
  if not partition_certificate['pooled_four_group_ordering_not_invariant_to_equal_position_weighting']:
   raise RuntimeError('partition-equal contrast no longer crosses zero')
  state['source_id_partition_count']=37
  state['equal_partition_contrast_crosses_zero']=True
  step('Fig6CD_source','fig6cd_source.py',['--feather',stress,'--rmd',rmd,'--out',reporter_source_dir])
  step('Fig6CD_RealTrace','fig6cd_aggregate.py',['--source',reporter_source_dir/'reporter_join_keys.csv.gz',
       '--realtrace',rt,'--receipt',rec,'--out',reporter_summary_dir])
  step('Fig6D_date_dose_audit','date_dose_audit.py',[
       '--source',reporter_summary_dir/'fig6d_candidate_observations.csv.gz',
       '--reference',REPO/'references/fig6d_group_means.csv',
       '--out',out/'FIG6D_DATE_DOSE_AUDIT'])
  cert=crosspanel(stress,cycle_dir/'fig6b_cycle_ledger.csv',reporter_summary_dir/'fig6d_candidate_observations.csv.gz',crosspanel_dir)
  state['steps'].append({'name':'Fig6BCD_crosspanel','returncode':0})
  state['crosspanel']=cert
  import pandas as pd
  calculated=pd.read_csv(reporter_summary_dir/'fig6d_group_means.csv')
  known=pd.read_csv(REPO/'references/fig6d_group_means.csv')
  comp=calculated.merge(known,on=['prom','CAM','condition'],suffixes=('_calculated','_known'),validate='one_to_one')
  if len(comp)!=8:raise RuntimeError('D eight group reference join incomplete')
  for f in ('mean_C','sd_C'):
   if (comp[f+'_calculated']-comp[f+'_known']).abs().max()>5e-8:
    raise RuntimeError('D numeric reference mismatch: '+f)
  if not(comp['N_corr_nonmissing_calculated']==comp['N_corr_nonmissing_known']).all():raise RuntimeError('D effective sample counts changed')
  state['D_reference_8_group_fields_within_5e-8']=True
  c=pd.read_csv(reporter_summary_dir/'fig6c_time_series.csv')
  if len(c)!=616:raise RuntimeError('Figure6C count mismatch')
  if not (crosspanel_dir/'crosspanel_software_id_crosswalk.csv').exists():raise RuntimeError('Bridge absent')
  state['Fig6C_timepoints']=len(c);state['Fig6D_groups']=len(calculated)
  bounds_dir=out/'MISSING_OUTCOME_BOUNDS'
  step('sharp_binary_completion','sharp_bounds.py',[
       '--groups',cycle_dir/'fig6b_group_counts.csv',
       '--crosswalk',crosspanel_dir/'crosspanel_software_id_crosswalk.csv',
       '--out',bounds_dir])
  bounds_certificate=json.loads((bounds_dir/'missing_outcome_bounds_receipt.json').read_text())
  if len(bounds_certificate['date_dose_matched_stratum_contrasts'])!=8:raise RuntimeError('strata mismatch')
  if any(not (v['sharp_lower']['decimal']<.5<v['sharp_upper']['decimal'])
         for group in bounds_certificate['exploratory_auc_sharp_unknown_bounds'].values()
         for v in group.values() if isinstance(v,dict)):
      raise RuntimeError('AUC sharp-bound reference changed')
  state['missingness_computed_from_fresh_source']=True
  state['all_auc_intervals_include_half']=True
  state['success']=True
 except Exception as e:
  state['success']=False;state['error']=str(e);state['traceback']=traceback.format_exc(limit=6)
 finally:
  (out/'ANALYSIS_RUN_RECEIPT.json').write_text(json.dumps(state,indent=2,ensure_ascii=False,sort_keys=True)+'\n')
 if not state['success']:raise RuntimeError(state['error'])
 return state
if __name__=='__main__':
 p=argparse.ArgumentParser()
 g=p.add_mutually_exclusive_group(required=True)
 g.add_argument('--author-zip',help='Gervais2026_data.zip original official 842,127,932-byte ZIP')
 g.add_argument('--stress-feather',help='Already verified original 482,147,704-byte stress Feather')
 p.add_argument('--rmd',help='Original Gervais analysis Rmd (pre-extracted mode)')
 p.add_argument('--realtrace-csv',help='SHA-verified selected RealTrace CSV.gz (pre-extracted mode)')
 p.add_argument('--realtrace-receipt',help='SHA source receipt for RealTrace CSV.gz (pre-extracted mode)')
 p.add_argument('--out',required=True)
 a=p.parse_args();print(json.dumps(run(a),indent=2,ensure_ascii=False))
