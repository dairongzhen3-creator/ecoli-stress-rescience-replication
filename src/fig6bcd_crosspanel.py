#!/usr/bin/env python3
"""One-to-one software-cycle source linkage, NOT same physical cell or causal mediation.

Inputs: independently recomputed Fig6B cycle ledger, Fig6D pre-stress RealTrace
row table, and independently checked author stress Feather for 1153 time_start values.
"""
import argparse, hashlib, json
from pathlib import Path
import pandas as pd
import numpy as np
from feather_v1_readonly import FeatherV1

GROUPS=['Unperturbed','CAM_after_burst','CAM_throughout','rpoS_KO']

def run(feather,b_csv,d_gz,out):
 out=Path(out)
 if out.exists():raise FileExistsError('Refuse to overwrite cross-panel output')
 b=pd.read_csv(b_csv,dtype={'cell':str,'date':str})
 d=pd.read_csv(d_gz,compression='gzip',dtype={'cell':str,'date':str})
 if len(b)!=1153 or b.cell.nunique()!=len(b):raise RuntimeError('Figure6B cycle IDs invalid')
 if len(d)!=32940 or d.duplicated(['cell','prom','time_sec']).any():raise RuntimeError('Figure6D source rows invalid')
 if set(d.prom.unique())!={'bolA','E08'}:raise RuntimeError('Reporter group names changed')
 if d.Corr_selected.count()!=32400:raise RuntimeError('Source-corrected concentration count changed')
 if set(b.fate.unique())!={'Survive','No Growth','Unknown'}:raise RuntimeError('Outcome categories changed')
 fg=FeatherV1(feather)
 if fg.num_rows!=2101521:raise RuntimeError('Wrong stress Feather input')
 ix=b.source_row_1based.to_numpy(dtype=np.int64)-1
 b['time_start']=fg.col('time_start',ix)
 b['software_phase_stress_start']=b.m_start-18000
 identity=['cell','date','STRESS','CAM','condition']
 dm=d[identity].drop_duplicates()
 if dm.cell.nunique()!=len(dm):raise RuntimeError('D cell maps multiple source conditions')
 paired=b[identity].merge(dm,on='cell',validate='one_to_one',suffixes=('_B','_D'))
 for field in ('date','STRESS','CAM','condition'):
  if (paired[field+'_B']!=paired[field+'_D']).any():raise RuntimeError('Cross-panel same-ID source contradiction '+field)
 g=d.groupby(['cell','prom'],sort=True).Corr_selected.agg(['size','count','mean']).reset_index()
 m=g.pivot(index='cell',columns='prom',values=['size','count','mean']);m.columns=[f'{name}_{prom}' for name,prom in m.columns];m=m.reset_index()
 linked=b.merge(m,on='cell',how='left',validate='one_to_one',indicator='_merge_panel')
 linked['has_D_source']=linked._merge_panel.eq('both')
 linked['has_D_both']=linked[['mean_bolA','mean_E08']].notna().all(axis=1)
 linked['born_at_or_after_nominal_stress']=linked.time_start>=linked.software_phase_stress_start
 no_d=linked[~linked.has_D_source]
 if len(paired)!=1145 or linked.has_D_source.sum()!=1145 or linked.has_D_both.sum()!=1139:raise RuntimeError('Cross-panel counts changed')
 if len(no_d)!=8 or set(no_d.CAM)!={'rpoS_KO'} or set(no_d.date)!={'20230712'}:raise RuntimeError('Unmatched-cell provenance differs')
 if not no_d.born_at_or_after_nominal_stress.all():raise RuntimeError('Unmatched B-only birth times not as expected')
 if no_d.fate.value_counts().to_dict()!={'No Growth':5,'Survive':3}:raise RuntimeError('B-only fates differ')
 summary=[]
 for group in GROUPS:
  all_b=linked[linked.CAM==group];complete=all_b[all_b.has_D_both]
  known_b=all_b[all_b.fate!='Unknown'];known_complete=complete[complete.fate!='Unknown']
  summary.append(dict(CAM=group,N_B=len(all_b),N_has_D=int(all_b.has_D_source.sum()),
        N_both_D=int(all_b.has_D_both.sum()),B_known_N=len(known_b),B_known_survive=int(known_b.fate.eq('Survive').sum()),
        B_original_known_survival_frac=float(known_b.fate.eq('Survive').mean()),
        complete_known_N=len(known_complete),complete_known_survive=int(known_complete.fate.eq('Survive').sum()),
        complete_case_known_survival_frac=float(known_complete.fate.eq('Survive').mean()) if len(known_complete) else None))
 certificate={'analysis':'crosspanel_software_identity','B_cycles':len(b),'D_candidates':len(d),'B_D_cell_id_match':len(paired),
  'B_D_both_corr_valid':int(linked.has_D_both.sum()),'B_only_eight_count':len(no_d),
  'B_only_all_July2023_rpoS':True,'B_only_all_born_at_or_after_nominal_stress':True,
  'comparison_unit':'source software cell ID, not certified optical physical UID',
  'publication_conclusion':'software computation and observation-set audit only','source_verification_scope':'software_table_level'}
 out.mkdir(parents=True)
 linked.drop(columns=['_merge_panel']).to_csv(out/'crosspanel_software_id_crosswalk.csv',index=False,float_format='%.12g')
 pd.DataFrame(summary).to_csv(out/'crosspanel_complete_case_summary.csv',index=False,float_format='%.12g')
 (out/'crosspanel_replication_receipt.json').write_text(json.dumps(certificate,indent=2,sort_keys=True)+'\n')
 return certificate
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--feather',required=True);a.add_argument('--B',required=True);a.add_argument('--D',required=True);a.add_argument('--out',required=True)
 p=a.parse_args();print(json.dumps(run(p.feather,p.B,p.D,p.out),indent=2))
