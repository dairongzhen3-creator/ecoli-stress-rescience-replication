#!/usr/bin/env python3
"""Source-locked upstream eligibility and channel-audit for Gervais Fig. 6C/D.
This program explicitly DOES NOT reproduce Corr_conc or Fig.6C/D; those require the
separate RealTrace stress Feather table which is not among the local parent files.
"""
import argparse,csv,hashlib,json,mmap,os,sys,collections,gzip
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).parent))
from feather_v1_readonly import FeatherV1
SRC_HASH='52c9cf4c7607ecd32a621bb410c967ff5a198d3cdf95f37323008319bb9272ac'
RMD_HASH='7f915a34db6c3326a0be833bc1c04e57127d90684ce281c9fe4b0cd65a48de70'
RMD_C_START='```{r Pre_stress_oncentration_traces}'
RMD_D_START='```{r concentration_bar_plot}'
H2O2={b'H2O2_0.65mM',b'H2O2_0.75mM'}

def sha(path):
 with open(path,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def extract_rmd_chunk(text,start):
 assert text.count(start)==1,start
 return text.split(start,1)[1].split('```',1)[0]

class RawStrings:
 def __init__(self,view):
  self.view=view;self.fd=open(view.path,'rb');self.mm=mmap.mmap(self.fd.fileno(),0,access=mmap.ACCESS_READ)
  self._cols={}
 def close(self):self.mm.close();self.fd.close()
 def offsets(self,name):
  if name not in self._cols:
   c=self.view.cols[name];assert c['type_name']=='UTF8' and c['encoding']==0
   start_off=c['offset']+(((self.view.num_rows+7)//8+7)//8*8 if c['null_count'] else 0)
   arr=np.memmap(self.view.path,dtype='<i4',offset=start_off,shape=(self.view.num_rows+1,),mode='r')
   data=((start_off+(self.view.num_rows+1)*4+7)//8)*8
   assert data+int(arr[-1])<=c['offset']+c['total_bytes']
   self._cols[name]=(arr,data)
  return self._cols[name]
 def get_bytes(self,name,i):
  off,start=self.offsets(name);a=int(off[int(i)]);b=int(off[int(i)+1]);return self.mm[start+a:start+b]
 def col(self,name,indices):
  off,start=self.offsets(name);return [self.mm[start+int(off[int(i)]):start+int(off[int(i)+1])] .decode('utf8') for i in indices]

def valid_numeric(view,name,indices):
 c=view.cols[name];indices=np.asarray(indices,dtype=np.int64)
 if c['null_count']==0:return np.ones(len(indices),dtype=bool)
 b=np.memmap(view.path,dtype=np.uint8,offset=c['offset'],shape=((view.num_rows+7)//8,),mode='r')
 return (b[indices>>3]>>(indices&7)&1)!=0

def run(feather,rmd,out):
 feather=Path(feather);rmd=Path(rmd);out=Path(out)
 if out.exists():raise FileExistsError(f'Refusing to overwrite {out}')
 assert feather.stat().st_size==482147704 and sha(feather)==SRC_HASH,'source FEATHER invalid'
 assert rmd.stat().st_size==207248 and sha(rmd)==RMD_HASH,'source RMD invalid'
 text=rmd.read_text(encoding='utf8')
 c=extract_rmd_chunk(text,RMD_C_START);d=extract_rmd_chunk(text,RMD_D_START)
 # Do not silently relax source-derived code gates.
 anchors_c=['bolA_as_ch_1()', 'left_join(myframes_RealTrace_stress_exp)', 'filter(!is.na(Corr_conc))','Corr_conc_maturation', 'group_by(prom, condition, CAM, time_sec)', 'mean(Corr_conc, na.rm=T)', 'sd(Corr_conc, na.rm=T)', 'filter(prom == "bolA")']
 anchors_d=['bolA_as_ch_1()', 'left_join(myframes_RealTrace_stress_exp)', 'filter(c_step %in% c("GM.1", "SPM.1"))', 'Corr_conc_maturation', 'group_by(prom, condition, CAM, c_step)', 'Corr_conc[time_sec > m_end - 2*3600]','filter(c_step == "SPM.1")']
 assert all(a in c for a in anchors_c),[a for a in anchors_c if a not in c]
 assert all(a in d for a in anchors_d),[a for a in anchors_d if a not in d]
 for k in ['mean_C = mean(Corr_conc[time_sec > m_end - 2*3600], na.rm=T)','sd_C = sd(Corr_conc[time_sec > m_end - 2*3600], na.rm=T)']:
  assert k in d,k
 v=FeatherV1(feather);assert v.num_rows==2101521 and len(v.cols)==26
 cats=v.category_levels('c_step');assert set(['GM.1','SPM.1','GM.2','STRESS.1']).issubset(cats)
 phase_num=v.col('c_step')
 # Figure 6D uses only GM.1/SPM.1. Figure 6C's visible endpoint at 20h
 # also permits the first STRESS.1 software tick (stress start), so preserve it.
 early=np.isin(phase_num,[cats.index('GM.1'),cats.index('SPM.1')])
 stress_first=(phase_num==cats.index('STRESS.1'))&(v.col('time_sec')==v.col('m_start'))
 ix=np.flatnonzero(early|stress_first)
 reader=RawStrings(v)
 try:
  is_stress=np.fromiter((reader.get_bytes('STRESS',i) in H2O2 for i in ix),dtype=bool,count=len(ix))
  ix=ix[is_stress]
  # Original Rmd uses bolA_as_ch_1: duplicate secondary reporter when prom2 is bolA.
  p1=reader.col('prom',ix);p2=reader.col('prom2',ix)
  b1=np.asarray([p in ('bolA','E08') for p in p1]);b2=np.asarray([p=='bolA' for p in p2])
  # primary E08/bolA and secondary bolA are the exact two reporters for stress Fig6D;
  # other primary promoters remain in source; do not misclassify them as reporter replicates.
  i1=ix[b1]; i2=ix[b2]
  desc=pd.DataFrame({
   'source_row_1based':np.r_[i1,i2]+1,
   'reporter_role':['primary']*len(i1)+['secondary_bolA']*len(i2),
   'prom':[p1[j] for j in np.flatnonzero(b1)]+['bolA']*len(i2),
   'Fluo1':reader.col('Fluo1',i1)+reader.col('Fluo2',i2),
   'date':np.r_[v.col('date',i1),v.col('date',i2)],
   'phase':[cats[j] for j in np.r_[v.col('c_step',i1),v.col('c_step',i2)]],
   'condition':reader.col('condition',i1)+reader.col('condition',i2),
   'STRESS':reader.col('STRESS',i1)+reader.col('STRESS',i2),
   'cell':reader.col('cell',i1)+reader.col('cell',i2),
   'ugl':reader.col('ugl',i1)+reader.col('ugl',i2),
   'time_sec':np.r_[v.col('time_sec',i1),v.col('time_sec',i2)],
   'm_start':np.r_[v.col('m_start',i1),v.col('m_start',i2)],
   'm_end':np.r_[v.col('m_end',i1),v.col('m_end',i2)],
   'DrpoS':np.r_[v.col('DrpoS',i1),v.col('DrpoS',i2)],
   'plasmid':np.r_[v.col('plasmid',i1),v.col('plasmid',i2)],
   'fluo_ampl_raw_source':np.r_[v.col('fluo_ampl_ch_1',i1),v.col('fluo_ampl_ch_2',i2)],
   'length_um':np.r_[v.col('length_um',i1),v.col('length_um',i2)]
  })
  for n,s in [('primary',i1),('secondary_bolA',i2)]:
   mask=desc.reporter_role.eq(n).to_numpy(); fl='fluo_ampl_ch_1' if n=='primary' else 'fluo_ampl_ch_2'
   valid=valid_numeric(v,fl,s)
   desc.loc[mask,'fluo_valid']=valid
  desc['fluo_valid']=desc['fluo_valid'].astype(bool)
    # Actual source label mirrors Rmd via first-CAM 0/16 and rpoS flag; do not rely solely on labels.
  ixall=np.r_[i1,i2]
  z=reader.col('CAM_C_INIT_SP',ixall)
  is_null=valid_numeric(v,'CAM_C_INIT_SP',ixall)
  # Source changes label at early starvation only; read exactly as written.
  cam=pd.Series(z,dtype='str').where(is_null,'<NULL>').fillna('<NULL>')
  dr=desc.DrpoS.to_numpy(dtype=bool)
  isc=desc.condition.str.contains('CAM',regex=False).to_numpy()
  desc['CAM']=np.select([isc&(cam=='0µM'),isc&(cam=='16µM'),(~isc)&(~dr),(~isc)&dr],
     ['CAM_after_burst','CAM_throughout','Unperturbed','rpoS_KO'],default='INVALID')
  assert not desc.CAM.eq('INVALID').any(),desc.loc[desc.CAM.eq('INVALID'),['condition','date','CAM']].head()
  # Reproduce declared D window: same source `c_step == SPM.1`, relative time from SPM.1 start,
  # and cutoff raw time_sec > raw m_end - 2h.  Never infer actual local medium arrival.
  desc['d_final_2h']=(desc.phase.eq('SPM.1')) & (desc.time_sec>desc.m_end-7200)
  desc['c_visible_nominal']=((desc.time_sec-desc.m_start)>=-7200)&((desc.time_sec-desc.m_start)<=72000) # pre-RealTrace approximation, only audit
  # For C, anchor relative to the unique SPM.1 start in each condition, not each row's own m_start.
  start=desc.loc[desc.phase.eq('SPM.1')].groupby('condition').m_start.nunique()
  assert (start==1).all(),f'Nonunique SPM.1 starts: {start[start!=1]}'
  mp=desc.loc[desc.phase.eq('SPM.1')].groupby('condition').m_start.first().to_dict()
  desc['time_rel_sp']=desc.time_sec-desc.condition.map(mp)
  desc['c_visible_nominal']=(desc.time_rel_sp>-7200)&(desc.time_rel_sp<=72000)&desc.prom.eq('bolA')
  # Important: no RealTrace source file and no Corr_conc; raw fluorescence/length are only
  # upstream source diagnostics; they are NOT identical with normalized corrected concentration.
  desc['raw_AU_per_um']=np.where(desc.fluo_valid&(desc.length_um>0),desc.fluo_ampl_raw_source/desc.length_um,np.nan)
  all_counts=[]
  grouping=['date','STRESS','CAM','prom','Fluo1','phase','reporter_role']
  for key,grp in desc.groupby(grouping,sort=True,dropna=False):
   dgr=grp[grp.d_final_2h]
   row={k:val for k,val in zip(grouping,key)}
   row.update(source_rows=len(grp),unique_software_cells=grp.cell.nunique(),unique_software_lanes=grp.ugl.nunique(),
              total_valid_raw_signal=int(grp.fluo_valid.sum()),
              d_final2h_rows=int(len(dgr)),d_final2h_valid_raw=int(dgr.fluo_valid.sum()),
              d_final2h_distinct_software_cells=int(dgr.cell.nunique()),
              d_final2h_distinct_software_lanes=int(dgr.ugl.nunique()),
              d_final2h_nominal_ticks=int(dgr.time_sec.nunique()),
              d_final2h_raw_AU_per_um_mean=float(dgr.raw_AU_per_um.mean()) if dgr.raw_AU_per_um.notna().any() else None,
              d_final2h_raw_AU_per_um_sd=float(dgr.raw_AU_per_um.std(ddof=1)) if len(dgr)>1 else None)
   all_counts.append(row)
  # True D grouped by condition, CAM, prom, phase *across both dates and doses*; source codes treat each condition (might imply dose).
  selected=desc[(desc.phase=='SPM.1')&desc.d_final_2h].copy()
  for g,sub in [('all',desc),('D_preinfer',selected)]:
   if not len(sub):raise ValueError('source selection empty')
  summary={'source_feather_sha256':SRC_HASH,'source_rmd_sha256':RMD_HASH,
           'source_rows':v.num_rows,'source_rows_GM1_SPM1_H2O2':int((v.col('c_step',ix)!=cats.index('STRESS.1')).sum()),
           'source_rows_STRESS1_first_software_tick_H2O2':int((v.col('c_step',ix)==cats.index('STRESS.1')).sum()),
           'source_C_including_20h_boundary_candidate_rows':int(len(ix)),
           'source_reporter_rows_primary_E08_or_bolA':int(len(i1)),
           'source_reporter_rows_secondary_bolA':int(len(i2)),
           'source_fig6cd_upstream_expanded_rows':int(len(desc)),
           'source_D_final2h_rows':int(len(selected)),
           'source_C_visible_window_rows_from_two_reporters':int(desc.c_visible_nominal.sum()),
           'source_Rmd_RealTrace_joinkey_duplicate_rows':int(desc.duplicated(subset=['condition','phase','cell','Fluo1','DrpoS','prom','plasmid','time_sec']).sum()),
           'source_D_final2h_nonnull_raw_rows':int(selected.fluo_valid.sum()),
           'source_D_conditions':sorted(selected.CAM.unique().tolist()),
           'source_D_software_cells':int(selected.cell.nunique()),
           'source_D_software_lanes':int(selected.ugl.nunique()),
           'source_D_nominal_sample_times':sorted(float(y) for y in selected.time_rel_sp.unique()),
           'rmd_C_chunk_sha256':hashlib.sha256(c.encode()).hexdigest(),
           'rmd_D_chunk_sha256':hashlib.sha256(d.encode()).hexdigest(),
           'realtrace_table_present':False,'actual_Corr_conc_reproduced':False,'actual_Figure6CD_reproduced':False,
           'author_Corrected_Concentration_Requires':['Corr_conc','Corr_conc_maturation'],
           'current_input':'source-only MoMA-export fluorescence and phase/identity columns',
           'note':'Raw AU/um is a diagnostic only, NOT author RealTrace-derived FP/um, and must not be compared quantitatively with published C/D.','source_verification_scope':'software_table_level'}
  out.mkdir(parents=True,exist_ok=False)
  selected_source_columns=['source_row_1based','reporter_role','condition','phase','cell','Fluo1','DrpoS','prom','plasmid','time_sec','time_rel_sp','m_start','m_end','date','STRESS','CAM','fluo_valid']
  with (out/'reporter_join_keys.csv.gz').open('wb') as raw_output:
   with gzip.GzipFile(fileobj=raw_output,mode='wb',filename='',mtime=0,compresslevel=6) as zipped:
    desc[selected_source_columns].rename(columns={'phase':'c_step'}).to_csv(zipped,index=False)
  pd.DataFrame(all_counts).to_csv(out/'reporter_source_groups.csv',index=False,float_format='%.12g')
  selected[['source_row_1based','date','STRESS','CAM','condition','prom','Fluo1','reporter_role','cell','ugl','time_sec','time_rel_sp','fluo_valid','raw_AU_per_um']].sort_values(['date','STRESS','CAM','prom','cell','time_sec']).to_csv(out/'reporter_last2h_source.csv',index=False,float_format='%.12g')
  # check each condition at the data level independent of original model's RealTrace join
  selected.groupby(['CAM','prom','Fluo1']).agg(source_rows=('cell','size'),valid_raw=('fluo_valid','sum'),software_cells=('cell','nunique'),software_lanes=('ugl','nunique'),dates=('date','nunique'),conditions=('condition','nunique'),time_grid_points=('time_rel_sp','nunique')).reset_index().to_csv(out/'reporter_group_eligibility.csv',index=False)
  (out/'reporter_source_eligibility.json').write_text(json.dumps(summary,indent=2,sort_keys=True,ensure_ascii=False)+'\n',encoding='utf8')
  return summary
 finally:reader.close()

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--feather',required=True);p.add_argument('--rmd',required=True);p.add_argument('--out',required=True);a=p.parse_args();print(json.dumps(run(a.feather,a.rmd,a.out),indent=2))
