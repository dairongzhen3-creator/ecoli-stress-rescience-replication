#!/usr/bin/env python3
"""Recompute source-observation date×nominal-dose descriptive sensitivity.
These are not independent culture-unit uncertainties. Python standard library only.
"""
import argparse,csv,gzip,hashlib,json,math,statistics
from collections import defaultdict
from pathlib import Path

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def f12(v):return format(v,'.12g')

def analyze(source_csv_gz,out_dir,known_8=None):
    src=Path(source_csv_gz);out=Path(out_dir)
    if out.exists():raise FileExistsError('OUTPUT_EXISTS '+str(out))
    groups=defaultdict(list);software_ids=defaultdict(set)
    rows=0;unique=set()
    with gzip.open(src,'rt',newline='',encoding='utf-8') as f:
        dr=csv.DictReader(f)
        needed={'CAM','prom','date','STRESS','Corr_selected','cell'}
        if not needed.issubset(dr.fieldnames or []):raise ValueError('REQUIRED_COLUMNS_MISSING')
        for row in dr:
            rows+=1;key=(row['CAM'],row['prom'],row['date'],row['STRESS'])
            software_ids[key].add(row['cell']);unique.add(row['cell'])
            s=row['Corr_selected']
            if s and s.lower() not in ('nan','none'):
                v=float(s)
                if not math.isfinite(v):raise ValueError('NONFINITE_CORRECTED_CONCENTRATION')
                groups[key].append(v)
            elif key not in groups:groups[key]=[]
    if rows==0:raise ValueError('NO_SOURCE_ROWS')
    agg=defaultdict(list)
    strata=[]
    for key,vals in sorted(groups.items()):
        n=len(vals)
        if not n:raise ValueError('EMPTY_STRATUM '+str(key))
        mean=math.fsum(vals)/n;sd=statistics.stdev(vals) if n>1 else 0
        agg[key[:2]].append((n,mean,vals))
        strata.append({'CAM':key[0],'reporter':key[1],'date':key[2],
                       'nominal_H2O2':key[3],
                       'source_rows':rows_for_key(key,src) if False else None,
                       'valid_corrected_rows':n,'source_software_cells':len(software_ids[key]),
                       'mean_FP_per_um':mean,'observational_sd_FP_per_um':sd})
    # Distinct total rows per key, including missing corrected signals, from a second direct pass.
    num_rows=defaultdict(int)
    with gzip.open(src,'rt',newline='',encoding='utf-8') as f:
        for row in csv.DictReader(f):num_rows[(row['CAM'],row['prom'],row['date'],row['STRESS'])]+=1
    for row in strata:row['source_rows']=num_rows[(row['CAM'],row['reporter'],row['date'],row['nominal_H2O2'])]
    weights=[]
    for (cam,prom),ls in sorted(agg.items()):
        N=sum(n for n,_,_ in ls)
        pooled=sum(n*mu for n,mu,_ in ls)/N
        equal=sum(mu for _,mu,_ in ls)/len(ls)
        group_ids=set().union(*(software_ids[key] for key in software_ids if key[:2]==(cam,prom)))
        weights.append({'CAM':cam,'reporter':prom,'source_software_cells':len(group_ids),
                        'valid_corrected_rows':N,'source_row_weighted_mean_FP_um':pooled,
                        'equal_date_dose_mean_FP_um':equal,'equal_minus_pooled_percent':100*(equal/pooled-1),
                        'min_date_dose_mean_FP_um':min(mu for _,mu,_ in ls),
                        'max_date_dose_mean_FP_um':max(mu for _,mu,_ in ls),
                        'num_source_date_dose_strata':len(ls)})
    if known_8 is not None:
        with Path(known_8).open(newline='',encoding='utf-8') as f:known={(r['CAM'],r['prom']):r for r in csv.DictReader(f)}
        if set(known)!={(r['CAM'],r['reporter']) for r in weights}:raise AssertionError('D_GROUP_REFERENCE_KEYS_MISMATCH')
        for r in weights:
            k=(r['CAM'],r['reporter']);x=known[k]
            if int(x['N_corr_nonmissing'])!=r['valid_corrected_rows']:raise AssertionError('D_N_MISMATCH:'+str(k))
            if abs(float(x['mean_C'])-r['source_row_weighted_mean_FP_um'])>1e-8:raise AssertionError('D_MEAN_MISMATCH:'+str(k))
    out.mkdir(parents=True)
    def write(filename,rows_,headers):
        with (out/filename).open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=headers,lineterminator="\n");w.writeheader();
            for row in rows_:w.writerow({k:f12(v) if isinstance(v,float) else v for k,v in row.items()})
    write('date_dose_observation_audit.csv',strata,list(strata[0]))
    write('date_dose_equal_stratum_sensitivity.csv',weights,list(weights[0]))
    rec={'stage':'FIG6D_DATE_DOSE_OBSERVATION_SENSITIVITY','source_sha256':sha(src),
         'source_rows':rows,'valid_corrected_rows':sum(r['valid_corrected_rows'] for r in strata),
         'date_dose_strata':len(strata),'reporter_treatment_groups':len(weights),
         'source_software_cells':len(unique),'reference_checked':known_8 is not None,
         'interpretation':'Descriptive source record reweighting, not independently cultured uncertainty'}
    (out/'date_dose_observation_audit_receipt.json').write_text(json.dumps(rec,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    return rec

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--out',required=True);p.add_argument('--reference')
    a=p.parse_args();print(json.dumps(analyze(a.source,a.out,a.reference),indent=2))
if __name__=='__main__':main()
