#!/usr/bin/env python3
"""Exact sensitivity: exact binary-completion intervals at equal exported date-position weight.

This audits SOURCE software IDs. A dot-delimited cell ID is NOT a physical-cell barcode,
chip/culture ID, intervention assignment log, or randomized biological replicate.
The analytical population is the already audited Figure 6B GM.2 cycle ledger.
"""
from __future__ import annotations
import argparse
import csv
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re

# Derived from the prior original-source replay. This is NOT the raw Feather hash.
EXPECTED_LEDGER_SHA256 = 'f78729215c62d1685bb0af8cc8c8a9111cb3ba2935f186975b0e469e79835828'
EXPECTED_COUNTS = {
    'Unperturbed': (219, 4, 65),
    'CAM_after_burst': (108, 59, 65),
    'CAM_throughout': (42, 176, 110),
    'rpoS_KO': (6, 294, 5),
}
ORDER = ('Unperturbed','CAM_after_burst','CAM_throughout','rpoS_KO')
FATES = ('Survive','No Growth','Unknown')
CELL_RE = re.compile(r'^(\d{8})\.(\d+)\.(\d+)\.(\d+)$')


def digest(file: Path) -> str:
    h = hashlib.sha256()
    with file.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def fmtfrac(x: Fraction) -> dict:
    return {'numerator': x.numerator, 'denominator': x.denominator, 'decimal': float(x)}


def parse_ledger(file: Path, expected_hash: str | None=EXPECTED_LEDGER_SHA256):
    got = digest(file)
    if expected_hash is not None and got != expected_hash:
        raise ValueError(f'SOURCE_SHA_MISMATCH: expected {expected_hash}; found {got}')
    positions = defaultdict(lambda: {'S': 0, 'N': 0, 'U': 0, 'doses': set(), 'groups': set()})
    seen = set()
    group_counts = defaultdict(lambda: [0,0,0])
    with file.open('r', newline='', encoding='utf8') as f:
        rows = csv.DictReader(f)
        required={'cell','date','STRESS','CAM','fate'}
        if not rows.fieldnames or not required.issubset(set(rows.fieldnames)):
            raise ValueError(f'MISSING_COLUMNS: {sorted(required - set(rows.fieldnames or []))}')
        for row in rows:
            cell=row['cell']; m=CELL_RE.fullmatch(cell)
            if not m: raise ValueError('UNRECOGNIZED_SOFTWARE_CELL_KEY: '+repr(cell))
            if m.group(1) != row['date']: raise ValueError('DATE_KEY_DISAGREEMENT: '+cell)
            if cell in seen: raise ValueError('DUPLICATE_SOFTWARE_CELL_ID: '+cell)
            seen.add(cell)
            cam=row['CAM']; fate=row['fate']; stress=row['STRESS']
            if cam not in ORDER: raise ValueError('UNEXPECTED_ARM: '+repr(cam))
            if fate not in FATES: raise ValueError('UNEXPECTED_FATE: '+repr(fate))
            if stress not in ('H2O2_0.65mM','H2O2_0.75mM'):
                raise ValueError('UNEXPECTED_DOSE: '+repr(stress))
            key=(row['date'],m.group(2))
            pos=positions[key]; pos['groups'].add(cam); pos['doses'].add(stress)
            idx=FATES.index(fate)
            group_counts[cam][idx] += 1
            pos[('S','N','U')[idx]] += 1
    if not seen:raise ValueError('NO_SOURCE_CYCLES')
    return got, positions, dict(group_counts), len(seen)


def summarize(positions):
    # Exact Fraction endpoint for each exported date-position stratum.
    per_position=[]; per_group=defaultdict(list)
    arm_overlap=0; dose_overlap=0
    for (date,pos),v in sorted(positions.items(), key=lambda item:(item[0][0],int(item[0][1]))):
        if len(v['groups'])!=1: arm_overlap+=1
        if len(v['doses'])!=1: dose_overlap+=1
        if len(v['groups'])!=1 or len(v['doses'])!=1:
            raise ValueError('SOURCE_POSITION_HAS_MULTIPLE_ARM_OR_DOSE: '+date+'.'+pos)
        cam=next(iter(v['groups']));stress=next(iter(v['doses']))
        s,n,u=(v['S'],v['N'],v['U']);total=s+n+u
        if not total:raise ValueError('ZERO_POSITION_DENOMINATOR')
        lo=Fraction(s,total);hi=Fraction(s+u,total)
        per_group[cam].append((lo,hi,s,n,u))
        per_position.append({'date':date,'position_from_software_id':pos,'STRESS':stress,'CAM':cam,
           'S':s,'N':n,'U':u,'total':total,'lower':fmtfrac(lo),'upper':fmtfrac(hi)})
    group=[];intervals={}
    for cam in ORDER:
        pp=per_group[cam];k=len(pp)
        if not k:raise ValueError('ARM_HAS_NO_POSITIONS: '+cam)
        lo=sum((p[0] for p in pp),Fraction(0,1))/k
        hi=sum((p[1] for p in pp),Fraction(0,1))/k
        vals=tuple(sum(p[i+2] for p in pp) for i in range(3))
        group.append({'CAM':cam,'n_positions':k,'S':vals[0],'N':vals[1],'U':vals[2],
          'n_source_cycles':sum(vals),'lower':fmtfrac(lo),'upper':fmtfrac(hi)})
        intervals[cam]=(lo,hi)
    comp=[]
    for before,after in zip(ORDER,ORDER[1:]):
        l=intervals[before][0]-intervals[after][1]
        h=intervals[before][1]-intervals[after][0]
        comp.append({'arm_A':before,'arm_B':after,'difference_lower':fmtfrac(l),
         'difference_upper':fmtfrac(h),'all_binary_completions_preserve_A_gt_B':l>0})
    return per_position,group,comp,{'position_overlap_arms':arm_overlap,'position_overlap_doses':dose_overlap}


def write_csv(path, fields, rows):
    with path.open('x', newline='',encoding='utf8') as f:
        w=csv.DictWriter(f, fieldnames=fields,lineterminator='\n');w.writeheader()
        for r in rows:
            w.writerow({k:v if not isinstance(v,dict) else v['decimal'] for k,v in r.items() if k in fields})


def run(source,out,expected_hash=EXPECTED_LEDGER_SHA256,require_reference=True):
    source=Path(source);out=Path(out)
    if out.exists():raise FileExistsError('OUTPUT_EXISTS: '+str(out))
    received, positions, group_counts, n=parse_ledger(source,expected_hash)
    if require_reference:
        if n!=1153:raise ValueError(f'COUNT_DRIFT: 1153 expected; got {n}')
        for cam in ORDER:
            if tuple(group_counts.get(cam,())) != EXPECTED_COUNTS[cam]:
                raise ValueError('SOURCE_GROUP_COUNTS_DRIFT: '+cam)
    per_position,groups,contrasts,check=summarize(positions)
    if require_reference:
        if len(per_position)!=37:raise ValueError('SOURCE_POSITION_COUNT_DRIFT')
        if sum(r['n_positions'] for r in groups)!=37:raise ValueError('POSITION_COUNT_SUM_DRIFT')
        if sum(r['n_positions'] for r in groups if r['CAM'].startswith('CAM_'))!=17:
            raise ValueError('CAM_POSITION_COUNT_DRIFT')
        if check['position_overlap_arms'] or check['position_overlap_doses']:
            raise ValueError('SOURCE_POSITION_OVERLAP_DRIFT')
    out.mkdir(parents=True)  # after all source/estimator assertions passed
    write_csv(out/'source_id_partition_counts.csv',
       ['date','position_from_software_id','STRESS','CAM','S','N','U','total','lower','upper'],per_position)
    write_csv(out/'partition_weighted_group_bounds.csv',
       ['CAM','n_positions','n_source_cycles','S','N','U','lower','upper'],groups)
    write_csv(out/'partition_weighted_contrast_bounds.csv',
       ['arm_A','arm_B','difference_lower','difference_upper','all_binary_completions_preserve_A_gt_B'],contrasts)
    cert={
       'schema':'equal-partition-binary-completion-v1',
       'source_ledger_sha256':received, 'source_ledger_is_source_derived_not_original_feather':True,
       'source_ledger_record_count':n,'date_position_units':len(per_position),
       'group_position_units':{g['CAM']:g['n_positions'] for g in groups},
       'position_source_structure':check,
       'method':'Position equal-weighted mean of software-source recovery rates; each unresolved software label independently completed as Survive or No Growth',
       'mathematical_proof':'All coefficients on each missing binary indicator are strictly nonnegative within a group; endpoints are attained by completing all Unknown as 0 or all as 1, respectively. Group contrast endpoints combine the two independently attainable arm endpoints.',
       'sharp_intervals':{g['CAM']:{'lower':g['lower'],'upper':g['upper']} for g in groups},
       'sharp_adjacent_differences':contrasts,
       'pooled_four_group_ordering_not_invariant_to_equal_position_weighting':not contrasts[0]['all_binary_completions_preserve_A_gt_B'],
       'not_a_biological_randomized_experiment':True,
       'not_independent_culture_replicates':True,
       'not_a_statistical_confidence_interval':True,
       'not_identified_causal_effect':True,
       'scope':'finite-record software-output sensitivity'
    }
    (out/'partition_weighted_bounds_receipt.json').write_text(json.dumps(cert,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf8')
    return cert

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',required=True,help='Fig6B source GM2 cycle ledger made by independent Python reimplementation')
    p.add_argument('--out',required=True,help='New, non-existing output directory')
    a=p.parse_args()
    c=run(a.source,a.out)
    print(json.dumps({'success':True,'source_ledger_sha256':c['source_ledger_sha256'],
     'n_positions':c['date_position_units'], 'n_cycles':c['source_ledger_record_count'],
     'pooled_order_sensitivity':c['pooled_four_group_ordering_not_invariant_to_equal_position_weighting']},indent=2))
