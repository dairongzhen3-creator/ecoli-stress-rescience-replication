#!/usr/bin/env python3
"""Independent numeric reproduction of Gervais et al. (2026) Fig. 6C/D.
Reads frozen source software provenance rows and genuine Colab-extracted RealTrace
corrected-concentration rows. Never treats software rows as biological replicates.
No causal identification or physical sampling/valve validation is claimed.
"""
import argparse, csv, gzip, hashlib, io, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

KEY=['condition','c_step','cell','Fluo1','DrpoS','prom','plasmid','time_sec']
REF_SHA={
 'source':'1dbb1c35be05c05c214b805819ca90721fa1c2f25a0891fca968111517a58c53',
 'realtrace':'85623f172105e349556a304a222ae2142c9f7775a3aef649f22f114fd53ffbb3'
}


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for c in iter(lambda: f.read(1<<20),b''): h.update(c)
    return h.hexdigest()


def norm_bool(s):
    v=s.astype(str).str.lower().str.strip().map({'true':'1','false':'0','1':'1','0':'0','t':'1','f':'0'})
    if v.isna().any(): raise ValueError('Bad Boolean key value')
    return v


def export_csv(df,path):
    df.to_csv(path,index=False,float_format='%.12g',na_rep='')


def export_gzip_csv(df,path):
    # Stable compressed outputs: fixed gzip timestamp, no original filename metadata.
    with open(path,'wb') as b:
        with gzip.GzipFile(fileobj=b,mode='wb',filename='',mtime=0,compresslevel=6) as g:
            with io.TextIOWrapper(g,encoding='utf8',newline='') as t:
                df.to_csv(t,index=False,float_format='%.12g',na_rep='')


def make_plot_c(c,op):
    plt.figure(figsize=(9,4.6))
    order=['Unperturbed','CAM_after_burst','CAM_throughout','rpoS_KO']
    labels={'Unperturbed':'Unperturbed','CAM_after_burst':'CAM after 5 h','CAM_throughout':'CAM throughout','rpoS_KO':'rpoS knockout'}
    for cam in order:
        d=c[c.CAM==cam].sort_values('time_rel_sp')
        x=d.time_rel_sp.to_numpy()/3600
        y=d.mean_C.to_numpy(); sd=d.sd_C.to_numpy()
        ln=plt.plot(x,y,label=labels[cam],linewidth=1.55)[0]
        plt.fill_between(x,np.maximum(y-sd,0),y+sd,alpha=.12,color=ln.get_color())
    plt.axvline(0,linestyle=':',linewidth=1)
    plt.xlabel('Hours from nominal starvation onset (author software clock)')
    plt.ylabel('Corrected bolA reporter concentration (FP/µm)')
    plt.xlim(-2,20);plt.ylim(bottom=0)
    plt.legend(loc='upper left',fontsize=8,ncol=2);plt.grid(alpha=.17);plt.tight_layout()
    plt.savefig(op/'fig6c_source_replot.png',dpi=170)
    plt.close()


def make_plot_d(d,op):
    order=['Unperturbed','CAM_after_burst','CAM_throughout','rpoS_KO']
    labels={'Unperturbed':'Unperturbed','CAM_after_burst':'CAM after 5 h','CAM_throughout':'CAM throughout','rpoS_KO':'rpoS knockout'}
    for prom,name in [('bolA','bolA'),('E08','Med3 (E08)')]:
        da=d[d.prom==prom].set_index('CAM').loc[order]
        plt.figure(figsize=(8,3.5))
        yy=np.arange(4)
        plt.barh(yy,da.mean_C.to_numpy(),xerr=da.sd_C.to_numpy(),capsize=3,alpha=.8)
        plt.yticks(yy,[labels[k] for k in order]);plt.gca().invert_yaxis()
        plt.xlabel('Corrected reporter concentration (FP/µm)')
        plt.title('Fig. 6D numeric reproduction: '+name)
        plt.grid(axis='x',alpha=.15);plt.tight_layout()
        plt.savefig(op/f'fig6d_{prom.upper()}_REALTRACE_REPLOT.png',dpi=170)
        plt.close()


def run(source,realtrace,receipt,out,write_rows=True):
    out=Path(out)
    if out.exists(): raise FileExistsError('Refusing existing output directory')
    rec=json.loads(Path(receipt).read_text())
    # The 987,148-row source ledger has a fixed byte-identical checksum.
    # RealTrace output is source-verified through its *actual* run receipt:
    # independent pd/Arrow versions can serialize equivalent CSV bytes differently.
    for role,path in [('source',source),('realtrace',realtrace)]:
        h=sha(path)
        expected=REF_SHA['source'] if role=='source' else rec['selected_gzip_csv_sha256']
        if h!=expected: raise RuntimeError(f'{role} source SHA mismatch {h}; expected={expected}')
    if rec['selected_gzip_csv_bytes']!=Path(realtrace).stat().st_size:
        raise RuntimeError('RealTrace extraction receipt byte-size mismatch')
    if rec['member_name']!='Gervais2026_data/Exported_data/myframes_RealTrace_stress_exp.feather' or rec['member_sha256']!='17bdb0bc9e45f0aff5435f0a168d87850a3a29370ccf53449201fbb63f18b217':
        raise RuntimeError('Unrecognized original RealTrace member in extraction receipt')
    src=pd.read_csv(source,compression='gzip',low_memory=False)
    rt=pd.read_csv(realtrace,compression='gzip',low_memory=False)
    assert len(src)==987148 and len(rt)==1347150
    assert set(KEY+['time_rel_sp','date','STRESS','CAM','m_end','source_row_1based']).issubset(src.columns)
    assert set(KEY+['Corr_conc','Corr_conc_maturation']).issubset(rt.columns)
    for df in [src,rt]:
        for k in ('DrpoS','plasmid'):df[k]=norm_bool(df[k])
        for k in ('condition','c_step','cell','Fluo1','prom'):df[k]=df[k].astype(str)
        df['time_sec']=pd.to_numeric(df.time_sec).astype('float64')
    if src.duplicated(KEY).any() or rt.duplicated(KEY).any(): raise RuntimeError('Duplicate 8-key source or RealTrace row')
    merged=src.merge(rt[KEY+['Corr_conc','Corr_conc_maturation']],on=KEY,how='left',validate='one_to_one',indicator=True)
    if len(merged)!=len(src):raise RuntimeError('Non1to1 join')
    for k in ['Corr_conc','Corr_conc_maturation']:
        merged[k]=pd.to_numeric(merged[k],errors='coerce')
    merged['Corr_selected']=np.where(merged.Fluo1=='GFP',merged.Corr_conc,merged.Corr_conc_maturation)
    # Fig6C is an Rmd post-join Corr_conc nonmissing filter, THEN maturation substitution;
    # plotted range is -2 < t <= 20h. The source includes one first-stress tick.
    cw=(merged.prom=='bolA')&(merged.time_rel_sp > -7200)&(merged.time_rel_sp<=72000)
    c_ok=cw&merged.Corr_conc.notna()
    dfc=merged.loc[c_ok]
    c=dfc.groupby(['prom','condition','CAM','time_rel_sp'],dropna=False,sort=True).Corr_selected.agg(N_corr_nonmissing='count',mean_C='mean',sd_C='std').reset_index()
    # Fig6D requires SPM.1, the last 2h using author m_end relative threshold.
    dw=(merged.c_step=='SPM.1')&(merged.time_sec>merged.m_end-7200)
    dfd=merged.loc[dw]
    d=dfd.groupby(['prom','condition','CAM'],dropna=False,sort=True).Corr_selected.agg(N_corr_nonmissing='count',mean_C='mean',sd_C='std').reset_index()
    if len(c)!=616 or len(d)!=8:raise AssertionError((len(c),len(d)))
    if int((merged._merge=='both').sum())!=966777: raise AssertionError('Unexpected source match count')
    if int((dw&merged.Corr_selected.notna()).sum())!=32400: raise AssertionError('Unexpected D count')
    # Quantify RealTrace join losses; use counts, not imputation or MAR.
    missing=merged.assign(C_candidate=cw,D_candidate=dw,C_accepted=c_ok,match=(merged._merge=='both'),corrected=merged.Corr_selected.notna())
    gates=[]
    for lbl,mask in [('ALL_8KEY_INPUT',np.ones(len(missing),dtype=bool)),('FIG6C_PLOT_WINDOW',cw),('FIG6D_LAST2H',dw)]:
        sub=missing.loc[mask]
        gates.append(dict(window=lbl,candidates=len(sub),matched=int(sub.match.sum()),missing_key=int((~sub.match).sum()),realtrace_nonmissing=int(sub.corrected.sum()),source_cell_ids=sub.cell.nunique(),source_lane_ids=sub.cell.str.extract(r'^(.*)\.\d+$')[0].nunique()))
    gates=pd.DataFrame(gates)
    dgroups=dfd.groupby(['prom','CAM'],sort=True).agg(source_candidates=('cell','size'),matched=('Corr_selected','count'),unique_cell=('cell','nunique'),unique_lane=('cell',lambda s:s.str.extract(r'^(.*)\.\d+$')[0].nunique())).reset_index()
    # Four source strata for each treatment and reporter, not independent culture counts.
    strata=dfd.groupby(['prom','CAM','date','STRESS'],sort=True).Corr_selected.agg(N='count',mean='mean',sd='std').reset_index()
    pooled=d.rename(columns={'N_corr_nonmissing':'pooled_N','mean_C':'pooled_mean','sd_C':'pooled_sd'})
    eq=strata.groupby(['prom','CAM']).agg(four_strata=('mean','size'),equal_strata_mean=('mean','mean'),strata_mean_min=('mean','min'),strata_mean_max=('mean','max')).reset_index()
    weight=pooled.merge(eq,on=['prom','CAM'],validate='one_to_one')
    weight['equal_minus_pooled_pct']=(weight.equal_strata_mean/weight.pooled_mean-1)*100
    # Per-software-cell averaging is an alternative descriptive estimand only.
    cm=dfd.groupby(['prom','condition','CAM','cell'],sort=True).Corr_selected.mean().reset_index(name='cellmean')
    cellmeans=cm.groupby(['prom','CAM']).cellmean.agg(software_cell_N='count',cell_equal_mean='mean').reset_index()
    weight=weight.merge(cellmeans,on=['prom','CAM'],validate='one_to_one')
    weight['cell_equal_minus_pooled_pct']=(weight.cell_equal_mean/weight.pooled_mean-1)*100
    # Baseline / late reference marks for the exact reported source time indices.
    anchor=c[c.time_rel_sp.isin([0,14400,64800,71280])].copy()
    # C versus D prefilter significance on actual values.
    c_na_asym=int((cw & merged.Corr_conc.isna() & merged.Corr_conc_maturation.notna()).sum())
    d_na_asym=int((dw & merged.Corr_conc.isna() & merged.Corr_conc_maturation.notna()).sum())
    # matched RealTrace extraction selected rows originally contain both columns nonmissing
    rt_finite={k:int(np.isfinite(rt[k].to_numpy(dtype=float)).sum()) for k in ['Corr_conc','Corr_conc_maturation']}
    output={
        'status':'FIG6C_AND_FIG6D_REALTRACE_NUMERIC_REPLICATION_PASSED',
        'scientific_scope':'published software fluorescence concentration, descriptive pooled rows only',
        'source_keys_rows':len(src),'realtrace_extracted_rows':len(rt),
        'matched_8keys':int((merged._merge=='both').sum()),'unmatched_8keys':int((merged._merge=='left_only').sum()),
        'C_aggregated_points':len(c),'C_plot_source_rows':int(cw.sum()),'C_valid_corr_rows':int(c_ok.sum()),
        'D_groups':len(d),'D_source_rows':int(dw.sum()),'D_valid_corr_rows':int((dw&merged.Corr_selected.notna()).sum()),
        'C_raw_missing_mature_valid':c_na_asym,'D_raw_missing_mature_valid':d_na_asym,
        'realtrace_all_numeric_values_finite':rt_finite,
        'four_strata_per_condition_reporter':bool((weight.four_strata==4).all()),
        'largest_abs_pct_equal_vs_pooled':float(weight.equal_minus_pooled_pct.abs().max()),
        'source_keys_sha256':sha(source),'extracted_realtrace_sha256':sha(realtrace),
        'original_zip_sha256_colab_reported_not_local_reverified':rec['original_zip_sha256'],
        'original_realtrace_feather_sha256_colab_reported_not_local_reverified':rec['member_sha256'],
        'limitation':'RealTrace extracted source verified locally by gzip SHA, original 842MB ZIP/member SHA only verified in user-run Colab receipt; D observations not independent cultures; no physical U or Z intervention.'
    }
    out.mkdir(parents=True,exist_ok=False)
    export_csv(c,out/'fig6c_time_series.csv')
    export_csv(d,out/'fig6d_group_means.csv')
    export_csv(strata,out/'fig6d_date_dose_strata.csv')
    export_csv(weight,out/'fig6d_weighting_sensitivity.csv')
    export_csv(dgroups,out/'fig6d_denominators.csv')
    export_csv(gates,out/'reporter_join_coverage.csv')
    export_csv(anchor,out/'fig6c_time_anchors.csv')
    if write_rows:
        cols_c=['source_row_1based','date','STRESS','CAM','condition','cell','time_rel_sp','time_sec','Corr_conc','Corr_conc_maturation','Corr_selected']
        cols_d=['source_row_1based','date','STRESS','CAM','condition','cell','prom','Fluo1','time_sec','m_end','Corr_conc','Corr_conc_maturation','Corr_selected']
        export_gzip_csv(merged.loc[cw,cols_c],out/'fig6c_candidate_observations.csv.gz')
        export_gzip_csv(merged.loc[dw,cols_d],out/'fig6d_candidate_observations.csv.gz')
    # The standalone figure builder renders all publication figures from the audited
    # computed tables; do not render duplicative auxiliary plots here.
    (out/'reporter_replication_receipt.json').write_text(json.dumps(output,sort_keys=True,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps(output,sort_keys=True,ensure_ascii=False,indent=2))
    return output

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--realtrace',required=True);p.add_argument('--receipt',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(a.source,a.realtrace,a.receipt,a.out)
