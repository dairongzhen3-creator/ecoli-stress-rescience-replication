#!/usr/bin/env python3
"""Source-locked independent reimplementation of Gervais et al. 2026 Figure 6B aggregation.

Research method scope: author-data software cycle outcomes only; not an independent biological experiment,
not a certified culture-level causal analysis, and not a mechanistic inference.
Uses the tested Feather V1 metadata reader from archived 902 verified source reader only for binary decoding;
selection/aggregation/missingness arithmetic is independently rewritten here.
"""
from pathlib import Path
import argparse,csv,hashlib,json,sys
from fractions import Fraction
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from feather_v1_readonly import FeatherV1

SRC_SHA='52c9cf4c7607ecd32a621bb410c967ff5a198d3cdf95f37323008319bb9272ac'
RMD_SHA='7f915a34db6c3326a0be833bc1c04e57127d90684ce281c9fe4b0cd65a48de70'
EXPECTED={'Unperturbed':(219,4,65),'CAM_after_burst':(108,59,65),'CAM_throughout':(42,176,110),'rpoS_KO':(6,294,5)}
ORDER=['Unperturbed','CAM_after_burst','CAM_throughout','rpoS_KO']

def sha_file(path):
 with open(path,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def extract_binary_valid(v,name,idx):
 c=v.cols[name]
 if not c['null_count']:return np.ones(len(idx),dtype=bool)
 bitmap=np.memmap(v.path,dtype=np.uint8,offset=c['offset'],mode='r',shape=((v.num_rows+7)//8,))
 return ((bitmap[idx>>3]>>(idx&7))&1)!=0

def filter_source(v):
 cats=v.category_levels('c_step')
 gm=np.flatnonzero(v.col('c_step')==cats.index('GM.2'))
 # The same row-level filters as the published Rmd, but implemented via vectorized source arrays.
 gm=gm[(v.col('time_start',gm)<v.col('m_start',gm))&(v.col('cell_rank',gm)<3)]
 st=v.cols['STRESS']; offs=np.memmap(v.path,dtype='<i4',offset=st['offset'],mode='r',shape=(v.num_rows+1,))
 gm=gm[(offs[gm+1]-offs[gm])==11]
 payload_start=((st['offset']+4*(v.num_rows+1)+7)//8)*8
 payload=np.memmap(v.path,dtype=np.uint8,offset=payload_start,mode='r',shape=(int(offs[-1]),))
 pos=offs[gm].astype(np.int64)
 valid_stress=np.zeros(len(pos),dtype=bool)
 for label in [b'H2O2_0.65mM',b'H2O2_0.75mM']:
  assert len(label)==11
  equals=np.ones(len(pos),dtype=bool)
  for i,c in enumerate(label):equals &= (payload[pos+i]==c)
  valid_stress |= equals
 assert valid_stress.all(), 'Source contains nonstandard 11-byte stress label'
 return gm

def run(feather:Path,rmd:Path,out:Path,comparison:Path):
 if out.exists():raise FileExistsError('Refusing to overwrite output directory: '+str(out))
 assert feather.stat().st_size==482147704 and sha_file(feather)==SRC_SHA,'Wrong Feather bytes'
 assert rmd.stat().st_size==207248 and sha_file(rmd)==RMD_SHA,'Wrong author Rmd bytes'
 txt=rmd.read_text('utf-8')
 block=txt.split('```{r fraction_survival_bar_plot}',1)[1].split('```',1)[0]
 required=['filter(c_step == "GM.2"','filter(time_start < m_start)','filter(cell_rank < 3)','group_by(condition, CAM, fate, cell)','frac_survival','error_bar_survival = sqrt((frac_survival*(1-frac_survival))/(sum(N_cells)))']
 assert all(k in block for k in required),'Author Rmd figure code no longer matches archival anchor'
 x=FeatherV1(feather);assert x.num_rows==2101521 and len(x.cols)==26
 ix=filter_source(x)
 assert len(ix)==322154,f'Wrong Figure6B row set size: {len(ix)}'
 strings=x.col('cell',ix)
 # Use first occurrence of each cell cycle in row order, rather than verified source reader's sort-based selection.
 d=pd.DataFrame({'ix':ix,'cell':strings})
 reps=d.loc[~d.cell.duplicated(keep='first')].copy()
 n=len(reps);assert n==1153
 rows=reps.ix.to_numpy(np.int64)
 def col(name):return x.col(name,rows)
 d2=pd.DataFrame({
  'source_row_1based':rows+1,'cell':reps.cell.to_numpy(),
  'date':col('date'),'condition':col('condition'),
  'STRESS':col('STRESS'),'end_type':col('end_type'),
  'm_start':col('m_start'),'time_end':col('time_end'),
  'CAM_C_INIT_SP':col('CAM_C_INIT_SP'),'DrpoS':col('DrpoS')})
 valid=extract_binary_valid(x,'CAM_C_INIT_SP',rows)
 d2.loc[~valid,'CAM_C_INIT_SP']='<null>'
 cond=d2.condition.astype(str)
 a=cond.str.contains('CAM',regex=False)
 dr=d2.DrpoS.astype(bool)
 d2['CAM']=np.select([a&d2.CAM_C_INIT_SP.eq('0µM'),a&d2.CAM_C_INIT_SP.eq('16µM'),~a&~dr,~a&dr],
                       ['CAM_after_burst','CAM_throughout','Unperturbed','rpoS_KO'],default='INVALID')
 assert not d2.CAM.eq('INVALID').any()
 fate=np.select([d2.end_type.eq('div'),d2.end_type.isin(['eod','lysis'])],['Survive','No Growth'],default='Unknown')
 after=(d2.time_end>d2.m_start+72000)
 d2['fate']=np.where((fate=='Unknown')&after,'No Growth',fate)
 assert d2.cell.is_unique
 grp=d2.groupby(['CAM','fate'],dropna=False).size().unstack(fill_value=0).reindex(index=ORDER,columns=['Survive','No Growth','Unknown'],fill_value=0)
 outrows=[]
 for group in ORDER:
  s,nog,u=(int(q) for q in grp.loc[group])
  assert (s,nog,u)==EXPECTED[group]
  known=s+nog;alln=known+u
  p=Fraction(s,known);lo=Fraction(s,alln);hi=Fraction(s+u,alln)
  outrows.append({'group':group,'Survive':s,'No_Growth':nog,'Unknown':u,'N_known':known,'N_all':alln,
       'rate_known_only':float(p),'unknown_fraction':u/alln,'author_errorbar_SE_allN':float((float(p)*(1-float(p))/alln)**.5),
       'known_only_binomial_SE':float((float(p)*(1-float(p))/known)**.5),
       'missingness_lower':float(lo),'missingness_upper':float(hi),
       'rate_fraction_exact':str(p),'missingness_lower_fraction_exact':str(lo),'missingness_upper_fraction_exact':str(hi)})
 assert sum(q['N_all'] for q in outrows)==1153
 # Compare with frozen published-reconstruction, without using it as an input to the new aggregation.
 with open(comparison,encoding='utf-8') as f:old={r['CAM']:r for r in csv.DictReader(f)}
 assert set(old)==set(ORDER)
 for q in outrows:
  z=old[q['group']]
  assert (q['Survive'],q['No_Growth'],q['Unknown']) == (int(z['Survive']),int(z['No Growth']),int(z['Unknown']))
  for oldkey,newkey in [('fraction_survival','rate_known_only'),('author_error_bar_sqrt_pq_div_N_all','author_errorbar_SE_allN'),('known_only_binomial_se','known_only_binomial_SE'),('missing_lower','missingness_lower'),('missing_upper','missingness_upper')]:
   assert abs(float(z[oldkey])-q[newkey])<1e-14,(q['group'],oldkey)
 a1=next(q for q in outrows if q['group']=='CAM_after_burst')
 a0=next(q for q in outrows if q['group']=='CAM_throughout')
 gap=100*(a1['rate_known_only']-a0['rate_known_only'])
 lo=100*(a1['missingness_lower']-a0['missingness_upper'])
 hi=100*(a1['missingness_upper']-a0['missingness_lower'])
 day=d2.groupby(['date','STRESS','CAM','fate']).size().unstack(fill_value=0).reindex(columns=['Survive','No Growth','Unknown'],fill_value=0).reset_index()
 day['N_known']=day['Survive']+day['No Growth'];day['N_all']=day.N_known+day.Unknown
 day['rate_known_only']=day.Survive/day.N_known
 out.mkdir(parents=True,exist_ok=False)
 with (out/'fig6b_group_counts.csv').open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,fieldnames=list(outrows[0]));w.writeheader();w.writerows(outrows)
 d2[['source_row_1based','cell','date','STRESS','condition','CAM','end_type','fate','time_end','m_start']].sort_values('cell').to_csv(out/'fig6b_cycle_ledger.csv',index=False)
 day.to_csv(out/'fig6b_date_dose_counts.csv',index=False)
 summary={'status':'PASS_SOURCE_LOCKED_SOFTWARE_NUMERIC_REPLICATION','source_feather_sha256':SRC_SHA,'source_rmd_sha256':RMD_SHA,
 'source_feather_rows':x.num_rows,'eligible_GM2_H2O2_source_rows':len(ix),'software_cycle_objects':len(d2),'groups':outrows,
 'known_only_CAM_difference_pp':gap,'no_assumption_unknown_only_difference_range_pp':[lo,hi],
 'reference_group_table_sha256':sha_file(comparison),'original_Rmd_chunk_sha256':hashlib.sha256(block.encode()).hexdigest(),
 'data_reanalysis_independent_algorithm':True,'raw_source_independent_experiment':False,'physical_instrument_time_verified':False,
 'cluster_level_error_estimation_verified':False,
 'status_qualifier':'Numerical source-defined software cell-cycle re-implementation, NOT independent experimental replication or causal estimation.'}
 (out/'fig6b_replication_receipt.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True)+'\n',encoding='utf8')
 print(json.dumps({k:summary[k] for k in ['status','source_feather_rows','eligible_GM2_H2O2_source_rows','software_cycle_objects','known_only_CAM_difference_pp','no_assumption_unknown_only_difference_range_pp']},indent=2))
 return summary

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--feather',type=Path,required=True);p.add_argument('--rmd',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
 a=p.parse_args();run(a.feather,a.rmd,a.out,a.reference)
