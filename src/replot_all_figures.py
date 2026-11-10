#!/usr/bin/env python3
"""Recreate all six *software-level* paper figures from source-derived numerical summaries.

This is an entirely original visualization script; it neither reads raw images nor
recalculates upstream RealTrace inferred concentrations. No claims of physical
cell-ID validity or treatment causality are made by rendering these plots.
"""
from __future__ import annotations
import argparse
import csv
import json
import hashlib
from pathlib import Path
from reviewer_figures import make_plots

COMPUTED_NAMES={
 'groups':'FIG6B/fig6b_group_counts.csv',
 'date_dose':'MISSING_OUTCOME_BOUNDS/date_dose_binary_completion_bounds.csv',
 'fig6c':'FIG6CD/fig6c_time_series.csv',
 'fig6d':'FIG6CD/fig6d_group_means.csv',
 'auc':'MISSING_OUTCOME_BOUNDS/missing_outcome_bounds_receipt.json',
}
SOURCE_NAMES={
 'groups':'fig6b_group_counts.csv',
 'date_dose':'date_dose_binary_completion_bounds.csv',
 'fig6c':'fig6c_time_series.csv',
 'fig6d':'fig6d_group_means.csv',
 'auc':'missing_outcome_bounds_receipt.json',
}
LABELS={'Unperturbed':'Unperturbed','CAM_after_burst':'CAM after 5 h',
        'CAM_throughout':'CAM throughout','rpoS_KO':r'$\Delta\mathit{rpoS}$'}
ORDER=list(LABELS)
EXPECTED_SOURCE_SHA256={
 'stress':'52c9cf4c7607ecd32a621bb410c967ff5a198d3cdf95f37323008319bb9272ac',
 'rmd':'7f915a34db6c3326a0be833bc1c04e57127d90684ce281c9fe4b0cd65a48de70',
 'selected_realtrace':'85623f172105e349556a304a222ae2142c9f7775a3aef649f22f114fd53ffbb3',
}

def sha256_file(path:Path):
 h=hashlib.sha256()
 with Path(path).open('rb') as file:
  for block in iter(lambda:file.read(1<<20),b''):h.update(block)
 return h.hexdigest()

def load_csv(f):
 with Path(f).open(encoding='utf8',newline='') as stream:
  return list(csv.DictReader(stream))

def computed_source_files(run_dir:Path):
 """Resolve source tables from a *successful*, newly computed seven-stage run.

 This path deliberately does not read the repository's frozen reference tables.
 They are used only for separate comparison/smoke tests.
 """
 run_dir=Path(run_dir)
 verdict=run_dir/'ANALYSIS_RUN_RECEIPT.json'
 if not verdict.is_file():raise FileNotFoundError('Missing fresh-run verdict: '+str(verdict))
 state=json.loads(verdict.read_text(encoding='utf8'))
 expected={'Fig6B','equal_partition_binary_completion','Fig6CD_source',
           'Fig6CD_RealTrace','Fig6D_date_dose_audit',
           'Fig6BCD_crosspanel','sharp_binary_completion'}
 actual={s.get('name') for s in state.get('steps',[]) if s.get('returncode')==0}
 if not state.get('success') or not expected.issubset(actual):
  raise ValueError('Fresh seven-stage computed source run has not passed')
 if state.get('source_sha256')!=EXPECTED_SOURCE_SHA256:
  raise ValueError('Fresh-run verdict source hashes do not match certified original inputs')
 resolved={key:run_dir/name for key,name in COMPUTED_NAMES.items()}
 for key,path in resolved.items():
  if not path.is_file():raise FileNotFoundError(f'Fresh computed input absent for {key}: {path}')
 return resolved

def plot_all(refdir:Path,outdir:Path,*,computed_run:bool=False):
 """Render figures from source-computed outputs or explicitly from frozen references.

 In the public one-command workflow, computed_run=True is mandatory. A
 references-only plot is useful for visual smoke tests but is NOT source replay.
 """
 refdir=Path(refdir)
 source=computed_source_files(refdir) if computed_run else {key:refdir/name for key,name in SOURCE_NAMES.items()}
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 import numpy as np
 outdir=Path(outdir)
 if outdir.exists():raise FileExistsError('OUTPUT_EXISTS: refusing to overwrite '+str(outdir))
 for file in source.values():
  if not file.is_file():raise FileNotFoundError(str(file))
 outdir.mkdir(parents=True)
 made=make_plots(source['groups'],source['date_dose'],outdir)
 c=load_csv(source['fig6c']);d=load_csv(source['fig6d'])
 if len(c)!=616 or len(d)!=8:raise ValueError('Expected 616 Fig6C rows and eight Fig6D groups')
 # This x variable is the authors' nominal software time, divided by 3600.
 fig,ax=plt.subplots(figsize=(10.2,5.0),layout='constrained')
 for arm in ORDER:
  rows=sorted((r for r in c if r['CAM']==arm),key=lambda r:float(r['time_rel_sp']))
  if len(rows)!=154:raise ValueError('Expected 154 Figure6C time rows for '+arm)
  x=np.array([float(r['time_rel_sp'])/3600 for r in rows]);mu=np.array([float(r['mean_C']) for r in rows]);sd=np.array([float(r['sd_C']) for r in rows])
  ln,=ax.plot(x,mu,linewidth=1.8,label=LABELS[arm]);ax.fill_between(x,mu-sd,mu+sd,color=ln.get_color(),alpha=.13,linewidth=0)
 ax.axvline(0,linestyle=':',color='.4',linewidth=1)
 ax.set(xlabel='Hours from nominal starvation onset (software clock)',ylabel='Corrected bolA reporter concentration (FP/µm)',title='')
 ax.set_xlim(-2.0,20.0);ax.margins(x=0.03);ax.legend(loc='upper left',fontsize=8,ncol=2);ax.grid(alpha=.12)
 fn=outdir/'figure2_bola_time_series.png';fig.savefig(fn,dpi=180,bbox_inches="tight",pad_inches=0.2);plt.close(fig);made.append(str(fn))
 for reporter,name,fname in [('bolA','bolA','figure3_bola_reporter_means.png'),('E08','Med3','figure4_med3_reporter_means.png')]:
  rows={r['CAM']:r for r in d if r['prom']==reporter}
  if set(rows)!=set(ORDER):raise ValueError('Missing Figure6D group '+name)
  means=[float(rows[arm]['mean_C']) for arm in ORDER];sds=[float(rows[arm]['sd_C']) for arm in ORDER]
  fig,ax=plt.subplots(figsize=(8.8,4.5),layout='constrained')
  xx=list(range(4));ax.barh(xx,means,xerr=sds,color='.66',edgecolor='.26',ecolor='.24',capsize=3)
  ax.set_yticks(xx,[LABELS[x] for x in ORDER]);ax.invert_yaxis();ax.set(xlabel='Corrected reporter concentration (FP/µm)',title='');ax.grid(axis='x',alpha=.13)
  fn=outdir/fname;fig.savefig(fn,dpi=180,bbox_inches="tight",pad_inches=0.2);plt.close(fig);made.append(str(fn))
 cert=json.loads(source['auc'].read_text(encoding='utf8'))
 ac=cert['exploratory_auc_sharp_unknown_bounds'];arms=['CAM_after_burst','CAM_throughout'];proms=['bolA','Med3']
 fig,ax=plt.subplots(figsize=(8.1,4.5),layout='constrained');w=.36
 for k,prom in enumerate(proms):
  xx=np.array([0,1])+((k-.5)*w)
  yy=[float(ac[arm][prom]['known_only']['decimal']) for arm in arms]
  bars=ax.bar(xx,yy,width=w,label=prom)
  for bar,y in zip(bars,yy):ax.text(bar.get_x()+w/2,y+.015,f'{y:.3f}',ha='center',fontsize=9)
 ax.set_xticks([0,1],['CAM after 5 h','CAM throughout']);ax.set_ylim(0,1.0);ax.axhline(.5,color='.35',linestyle='--',linewidth=1)
 ax.set(ylabel='Descriptive known-outcome ROC AUC',title='')
 ax.legend(fontsize=9);ax.grid(axis='y',alpha=.12)
 fn=outdir/'figure6_reporter_rank_auc.png';fig.savefig(fn,dpi=180,bbox_inches="tight",pad_inches=0.2);plt.close(fig);made.append(str(fn))
 if computed_run:
  # Captures EXACT derived source bytes read when figures were rendered. This
  # is auditable provenance, not an independent experiment or formal journal QA.
  manifest={
   'schema':'source-figure-link-v1',
   'mode':'FRESH_COMPUTED_RUN',
   'run_verdict_sha256':sha256_file(refdir/'ANALYSIS_RUN_RECEIPT.json'),
   'source_inputs_sha256':{key:sha256_file(path) for key,path in source.items()},
   'figure_png_sha256':{Path(path).name:sha256_file(Path(path)) for path in made},
   'source_hash_lock':EXPECTED_SOURCE_SHA256,
   'computational_scope':'software-output level; not biological replication',
  }
  (outdir/'FIGURE_SOURCE_SHA256_RECEIPT.json').write_text(
    json.dumps(manifest,indent=2,sort_keys=True)+'\n',encoding='utf8')
 return made

if __name__=='__main__':
 p=argparse.ArgumentParser()
 inputs=p.add_mutually_exclusive_group(required=True)
 inputs.add_argument('--computed-run',help='fresh seven-stage output dir (required for source-to-figure proof)')
 inputs.add_argument('--references',help='small frozen summary tables (visual smoke test only)')
 p.add_argument('--out',required=True)
 args=p.parse_args()
 if args.computed_run:
  images=plot_all(Path(args.computed_run),Path(args.out),computed_run=True)
  mode='FRESH_COMPUTED_RUN'
 else:
  images=plot_all(Path(args.references),Path(args.out))
  mode='FROZEN_REFERENCE_SMOKE_ONLY'
 print(json.dumps({'success':True,'mode':mode,'figure_count':len(images),'figures':images},indent=2))
