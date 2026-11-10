#!/usr/bin/env python3
"""Exact, source-partition-weighted binary-completion tipping audit.

An algebraic sensitivity calculation over observed software source keys; it is
not an estimated dropout process or a model of treatment assignment.
"""
from __future__ import annotations
import argparse,csv,json
from fractions import Fraction
from pathlib import Path

ARMS=('Unperturbed','CAM_after_burst','CAM_throughout','rpoS_KO')
EXPECTED_COUNTS={'Unperturbed':9,'CAM_after_burst':7,'CAM_throughout':10,'rpoS_KO':11}

def analyze(path):
    with Path(path).open(encoding='utf8',newline='') as f:rows=list(csv.DictReader(f))
    if len(rows)!=37:raise ValueError('Expected 37 source-ID partitions')
    by={k:[] for k in ARMS};keys=set()
    for r in rows:
        arm=r['CAM']
        if arm not in by:raise ValueError('Unknown treatment label')
        key=(r['date'],r['position_from_software_id'])
        if key in keys:raise ValueError('Duplicated date/second-ID component')
        keys.add(key)
        s,n,u,total=(int(r[k]) for k in ('S','N','U','total'))
        if min(s,n,u,total)<0 or s+n+u!=total or total==0:raise ValueError('Invalid source counts')
        by[arm].append((Fraction(s,total),Fraction(u,total),total,u,r['date']))
    if {k:len(v) for k,v in by.items()}!=EXPECTED_COUNTS:raise ValueError('Unexpected treatment partition counts')
    def avg(arm,j):return sum((r[j] for r in by[arm]),Fraction(0))/len(by[arm])
    a,b='Unperturbed','CAM_after_burst'
    if sum(x[3] for x in by[a])!=65 or sum(x[3] for x in by[b])!=65:raise ValueError('Unexpected unresolved outcomes')
    low=avg(a,0)-avg(b,0)-avg(b,1)
    upper=avg(a,0)+avg(a,1)-avg(b,0)
    if low != Fraction(-2539,3533145):raise ValueError('Extremal comparison has drifted')
    # Starting at the sharp minimum, reversing exactly one adversarial binary
    # completion adds a positive partition-dependent coefficient.
    steps=[Fraction(1,len(by[arm])*row[2]) for arm in (a,b) for row in by[arm] if row[3]]
    smallest=min(steps)
    if low+smallest<=0:raise ValueError('One-bit tipping property drifted')
    common_q0=avg(a,0)-avg(b,0)
    common_q1=common_q0+avg(a,1)-avg(b,1)
    if common_q0<=0 or common_q1<=0:raise ValueError('Common-propensity check drifted')
    if set(r[4] for r in by[a]) & set(r[4] for r in by[b]):raise ValueError('Expected cross-date comparison')
    d={
        'source_id_partition_count':len(rows),
        'partitions_by_arm':{k:len(v) for k,v in by.items()},
        'cross_date_comparison':True,
        'relevant_unresolved_labels':{'Unperturbed':65,'CAM_after_burst':65},
        'sharp_minimum_difference_fraction':str(low),
        'sharp_minimum_percentage_points':100*float(low),
        'sharp_maximum_percentage_points':100*float(upper),
        'smallest_single_label_shift_fraction':str(smallest),
        'smallest_single_label_shift_percentage_points':100*float(smallest),
        'minimum_after_one_opposite_extreme_flipped_pp':100*float(low+smallest),
        'reversal_requires_all_130_oppositely_assigned_unknown_labels':True,
        'illustrative_shared_unknown_regrowth_fraction_q0_pp':100*float(common_q0),
        'illustrative_shared_unknown_regrowth_fraction_q1_pp':100*float(common_q1),
        'shared_unknown_propensity_avoids_negative_expected_contrast':True,
        'scope':'finite-record descriptive scenario with equal source-ID partition weights',
        'assumption_note':'Shared q is an illustrative sensitivity model and not inferred from observed tracking.',
    }
    return d

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--partitions',required=True)
    p.add_argument('--out',required=True)
    a=p.parse_args();out=Path(a.out)
    if out.exists():raise FileExistsError('Refusing to overwrite existing audit output')
    result=analyze(a.partitions)
    out.mkdir(parents=True)
    (out/'tipping_point_summary.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
