#!/usr/bin/env python3
""" Exact finite-population missing-outcome sensitivity.

Takes freshly re-derived, hash-locked author-software tables from the common
independent replication runner; does not infer physical UID or causal intervention.
"""
import csv
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import platform
import sys

BASE = Path(__file__).resolve().parent
F1 = 'fig6b_group_counts.csv'
F2 = 'crosspanel_software_id_crosswalk.csv'
FROZEN_SHA = {
    F1: '41c8de5d427e538a218120f4f2f79a60aad737bd7c8ffc18287a7ef74675ce29',
    F2: 'a8f5f738e430789be668c23f2020cfba8bef9cd1ddaf50e04d7a4937c75fa6c9',
}
ORDER = ['Unperturbed','CAM_after_burst','CAM_throughout','rpoS_KO']

def frac(x):
    return dict(exact=f'{x.numerator}/{x.denominator}', decimal=float(x))

def rows(path):
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))

def midranks2(values):
    """Return 2*ascending 1-index midrank per input item, exact under ties."""
    order = sorted(range(len(values)), key=lambda i:values[i])
    ranks = [None]*len(values)
    i=0
    while i<len(order):
        j=i+1
        while j<len(order) and values[order[j]] == values[order[i]]:
            j+=1
        r2=i+j+1
        for idx in order[i:j]:
            ranks[idx]=r2
        i=j
    return ranks

def auc_known(score_label_pairs):
    r=[(float(s), y) for s,y in score_label_pairs if y in ('Survive','No Growth')]
    if len(r)==0: return None
    rr=midranks2([s for s,_ in r]); np=sum(y=='Survive' for _,y in r);nn=len(r)-np
    if not np or not nn: return None
    sr=sum(rank for rank,(_,y) in zip(rr,r) if y=='Survive')
    return Fraction(sr-np*(np+1),2*np*nn)

def sharp_auc_bounds(score_label_pairs):
    """Exact min/max empirical AUC over all binary assignments of Unknown.

    Formula (with ties): AUC = [sum(2*midrank_positive)
                                   - n_positive*(n_positive+1)]
                                   / (2*n_positive*n_negative).
    All complete-case observations are ranked once; for each k unknown positives,
    rank-sum is minimized/maximized by the k smallest/largest unknown midranks.
    Thus scanning k=0..U is globally sharp, not a uniform-assignment scenario.
    """
    p=list(score_label_pairs)
    values=[float(s) for s,_ in p]
    labels=[lab for _,lab in p]
    if any(y not in ('Survive','No Growth','Unknown') for y in labels):
        raise ValueError('unexpected label')
    ranks=midranks2(values)
    known_pos=sum(y=='Survive' for y in labels)
    known_neg=sum(y=='No Growth' for y in labels)
    unknown_ranks=sorted(rank for rank,y in zip(ranks,labels) if y=='Unknown')
    s_pos=sum(rank for rank,y in zip(ranks,labels) if y=='Survive')
    pu=[0]
    for rank in unknown_ranks:pu.append(pu[-1]+rank)
    U=len(unknown_ranks)
    bestlo=None;besthi=None;kle=None;khi=None
    for k in range(U+1):
        np=known_pos+k;nn=known_neg+U-k
        if not np or not nn:continue
        lo=Fraction(s_pos+pu[k]-np*(np+1),2*np*nn)
        hi=Fraction(s_pos+pu[U]-pu[U-k]-np*(np+1),2*np*nn)
        if bestlo is None or lo<bestlo:bestlo=lo;kle=k
        if besthi is None or hi>besthi:besthi=hi;khi=k
    if bestlo is None:raise ValueError('AUC undefined for all assignments')
    return {'n':len(p),'known_positive':known_pos,'known_negative':known_neg,'unknown':U,
            'known_only':frac(auc_known(p)) if auc_known(p) is not None else None,
            'sharp_lower':frac(bestlo),'sharp_upper':frac(besthi),
            'lower_achieved_unknown_positive_count':kle,
            'upper_achieved_unknown_positive_count':khi}

def run(group_csv, crosswalk_csv, out_dir):
    INP = {F1:Path(group_csv), F2:Path(crosswalk_csv)}
    OUT = Path(out_dir)
    if OUT.exists(): raise FileExistsError('Output directory exists; refusing overwrite: '+str(OUT))
    print('Exact finite-record binary-completion audit', flush=True)
    for name,h in FROZEN_SHA.items():
        got=sha256(INP[name].read_bytes()).hexdigest()
        if h!=got:raise ValueError(f'INTEGRITY_FAIL: {name}: {got}')
    OUT.mkdir(parents=True,exist_ok=False)
    source_groups=rows(INP[F1])
    source_crosswalk=rows(INP[F2])
    if len(source_groups)!=4 or len(source_crosswalk)!=1153:
        raise ValueError('source cardinality changed')
    if set(g['group'] for g in source_groups)!=set(ORDER):raise ValueError('group labels changed')
    if len(set(r['cell'] for r in source_crosswalk))!=1153:raise ValueError('software-ID duplicate')
    groupcount=defaultdict(Counter)
    groups=defaultdict(list)
    for r in source_crosswalk:
        groupcount[r['CAM']][r['fate']]+=1
        groups[r['CAM']].append(r)
    groupresults={}
    for g in source_groups:
        name=g['group'];S=int(g['Survive']);N=int(g['No_Growth']);U=int(g['Unknown']);T=S+N+U
        if (S+N,T)!=(int(g['N_known']),int(g['N_all'])):raise ValueError('source denominator mismatch')
        assert (groupcount[name]['Survive'],groupcount[name]['No Growth'],groupcount[name]['Unknown'])==(S,N,U)
        known=Fraction(S,S+N);low=Fraction(S,T);high=Fraction(S+U,T)
        if abs(float(low)-float(g['missingness_lower']))>1e-10 or abs(float(high)-float(g['missingness_upper']))>1e-10:
            raise ValueError('Reference finite-record bounds do not match')
        groupresults[name]={
            'Survive':S,'No_Growth':N,'Unknown':U,'N_all':T,
            'known_only_survival':frac(known),
            'min_all_unknown_assigned_No_Growth':frac(low),
            'max_all_unknown_assigned_Survive':frac(high),
            'unknown_share':frac(Fraction(U,T))
        }
    contrasts=[]
    for a,b in zip(ORDER,ORDER[1:]):
        ga=groupresults[a];gb=groupresults[b]
        amin=Fraction(ga['Survive'],ga['N_all'])
        amax=Fraction(ga['Survive']+ga['Unknown'],ga['N_all'])
        bmin=Fraction(gb['Survive'],gb['N_all'])
        bmax=Fraction(gb['Survive']+gb['Unknown'],gb['N_all'])
        min_diff=amin-bmax
        max_diff=amax-bmin
        contrasts.append({'higher_group':a,'lower_group':b,
                          'sharp_min_difference':frac(min_diff),'sharp_max_difference':frac(max_diff),
                          'ordering_invariant_for_all_binary_unknown_labels':min_diff>0})
    # Date × nominal H2O2 dose restricts comparisons to overlapping observation strata.
    # Across the four pooled groups, CAM and non-CAM do NOT share a date/dose stratum.
    stratified=defaultdict(Counter)
    for r in source_crosswalk:
        stratified[(r['date'],r['STRESS'],r['CAM'])][r['fate']]+=1
    strata=[]
    for day,dose in sorted({(date,dose) for date,dose,cam in stratified}):
        present=[cam for date0,dose0,cam in stratified if date0==day and dose0==dose]
        if set(present)=={'CAM_after_burst','CAM_throughout'}:
            a,b='CAM_after_burst','CAM_throughout'
        elif set(present)=={'Unperturbed','rpoS_KO'}:
            a,b='Unperturbed','rpoS_KO'
        else:
            raise ValueError(f'unpaired or unexpected strata: {day} {dose} {present}')
        da=stratified[(day,dose,a)];db=stratified[(day,dose,b)]
        Sa,Na,Ua=[da[label] for label in ('Survive','No Growth','Unknown')]
        Sb,Nb,Ub=[db[label] for label in ('Survive','No Growth','Unknown')]
        Ta=Sa+Na+Ua;Tb=Sb+Nb+Ub
        mind=Fraction(Sa,Ta)-Fraction(Sb+Ub,Tb)
        maxd=Fraction(Sa+Ua,Ta)-Fraction(Sb,Tb)
        strata.append({'date':day,'nominal_dose':dose,'higher_group':a,'lower_group':b,
            'higher_S_N_U':[Sa,Na,Ua],'lower_S_N_U':[Sb,Nb,Ub],
            'lower_and_upper_fraction_difference':[frac(mind),frac(maxd)],
            'strict_order_identified_under_arbitrary_unknown':mind>0})
    assert len(strata)==8
    for cam in ORDER:
        for fate in ['Survive','No Growth','Unknown']:
            if sum(c[fate] for (date,dose,c_cam),c in stratified.items() if c_cam==cam)!=groupcount[cam][fate]:
                raise ValueError('stratum aggregates inconsistent')
    # Input identity guard: 8 B source objects without D, 1145 matched, 1139 both reporter valid.
    matched=sum(r['has_D_source']=='True' for r in source_crosswalk)
    both=sum(r['has_D_both']=='True' for r in source_crosswalk)
    if (matched,both)!=(1145,1139):raise ValueError('cross-panel inclusion changed')
    by_cam={}
    for cam in ['CAM_after_burst','CAM_throughout']:
        filtered=[r for r in groups[cam] if r['has_D_both']=='True']
        d={'N_complete_reporters':len(filtered)}
        for prom,col in [('bolA','mean_bolA'),('Med3','mean_E08')]:
            pair=[(float(r[col]),r['fate']) for r in filtered]
            if not all(math.isfinite(s) for s,_ in pair):raise ValueError('bad concentration')
            d[prom]=sharp_auc_bounds(pair)
        by_cam[cam]=d
    result={
        'scope':'AUTHOR_SOFTWARE_CYCLE_LABEL_FINITE_POPULATION_ONLY',
        'status':'EXACT_BINARY_UNKNOWN_ASSIGNMENT_SENSITIVITY_NOT_CAUSAL',
        'provenance':'Source-derived tables from independently executed current runner; exact SHA matches validated frozen source-derived tables; full 842MB ZIP not verified in this run',
        'source_sha256':FROZEN_SHA,
        'source_counts':{'fig6b_software_cycles':1153,'shared_B_D_software_ids':matched,'both_reporters_valid':both},
        'group_bounds':groupresults,
        'adjacent_group_contrast_certificate':contrasts,
        'date_dose_matched_stratum_contrasts':strata,
        'exploratory_auc_sharp_unknown_bounds':by_cam,
        'nonclaims':['No new experiment, physical UID, survival assay or independent culture replicate',
                     'Known-outcome contrast does not identify clinical/biological causal effect',
                     'AUC extrema condition on having both reporters and on arbitrary binary annotation of software Unknown',
                     'No confidence interval, generalization or held-out prediction assertion'],
        'review_status':'not a journal peer review',
        'python':platform.python_version()
    }
    (OUT/'missing_outcome_bounds_receipt.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    with (OUT/'group_binary_completion_bounds.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['group','S','No_Growth','Unknown','N_all','known_only','sharp_low','sharp_high'])
        for name in ORDER:
            g=groupresults[name];w.writerow([name,g['Survive'],g['No_Growth'],g['Unknown'],g['N_all'],
                g['known_only_survival']['exact'],g['min_all_unknown_assigned_No_Growth']['exact'],g['max_all_unknown_assigned_Survive']['exact']])
    with (OUT/'auc_binary_completion_bounds.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['CAM','reporter','N_both','known_pos','known_neg','Unknown','AUC_known','AUC_sharp_low','AUC_sharp_high','k_at_low','k_at_high'])
        for cam,d in by_cam.items():
            for prom in ['bolA','Med3']:
                a=d[prom];w.writerow([cam,prom,a['n'],a['known_positive'],a['known_negative'],a['unknown'],
                    a['known_only']['decimal'],a['sharp_lower']['decimal'],a['sharp_upper']['decimal'],
                    a['lower_achieved_unknown_positive_count'],a['upper_achieved_unknown_positive_count']])
    with (OUT/'date_dose_binary_completion_bounds.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['date','STRESS','group_A','group_B','S_A','NoGrowth_A','Unknown_A','S_B','NoGrowth_B','Unknown_B','diff_low','diff_high','order_identified'])
        for st in strata:
            x=st['higher_S_N_U'];y=st['lower_S_N_U']
            w.writerow([st['date'],st['nominal_dose'],st['higher_group'],st['lower_group'],*x,*y,
                        st['lower_and_upper_fraction_difference'][0]['exact'],
                        st['lower_and_upper_fraction_difference'][1]['exact'],
                        st['strict_order_identified_under_arbitrary_unknown']])
    print(json.dumps({'groups':groupresults,'contrasts':contrasts,'matched_strata':strata,'AUC':by_cam},ensure_ascii=False,indent=2))
    return result

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(description='Exact missing-outcome bounds on freshly reproduced source tables')
    ap.add_argument('--groups',required=True,help='Fresh Figure 6B four group table')
    ap.add_argument('--crosswalk',required=True,help='Fresh Figure 6B/6D software cell crosswalk')
    ap.add_argument('--out',required=True,help='New non-existing results directory')
    a=ap.parse_args(); run(a.groups,a.crosswalk,a.out)

