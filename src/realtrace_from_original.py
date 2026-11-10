#!/usr/bin/env python3
"""Extract the *author's* complete RealTrace Feather into a selected table.

Run only on the original Feather extracted from the verified Zenodo archive.
The resulting CSV is derived, not independent experimental measurements. Input
requires pyarrow. Original source and attribution are maintained separately.
"""
from pathlib import Path
import argparse,csv,gzip,hashlib,io,json

KEYS=['condition','c_step','cell','Fluo1','DrpoS','prom','plasmid','time_sec']
FIELDS=KEYS+['Corr_conc','Corr_conc_maturation']
EXPECTED_MEMBER_SHA='17bdb0bc9e45f0aff5435f0a168d87850a3a29370ccf53449201fbb63f18b217'
EXPECTED_MEMBER_SIZE=264629352
EXPECTED_MEMBER_NAME='Gervais2026_data/Exported_data/myframes_RealTrace_stress_exp.feather'

def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(4<<20),b''):h.update(b)
 return h.hexdigest()

def extract(feather,output,receipt,original_zip_sha,member_crc='3bcd5583'):
 from importlib.util import find_spec
 if find_spec('pyarrow') is None:
  raise RuntimeError('pyarrow is required to read the original RealTrace Feather; install requirements-full.txt')
 import pandas as pd
 feather,output,receipt=map(Path,[feather,output,receipt])
 if output.exists() or receipt.exists():raise FileExistsError('Refuse to overwrite selected original derived files')
 if feather.stat().st_size!=EXPECTED_MEMBER_SIZE or sha(feather)!=EXPECTED_MEMBER_SHA:
  raise RuntimeError('Original RealTrace Feather source size/SHA mismatch')
 frame=pd.read_feather(feather,columns=FIELDS)
 select=frame.c_step.astype(str).isin(['GM.1','SPM.1','STRESS.1']) & (
     ((frame.prom.astype(str)=='E08')&(frame.Fluo1.astype(str)=='YFP'))|
     ((frame.prom.astype(str)=='bolA')&(frame.Fluo1.astype(str)=='RFP')))
 sub=frame.loc[select,FIELDS].copy()
 if len(sub)!=1347150:raise RuntimeError(f'Unexpected selected row count: {len(sub)}')
 if sub.duplicated(KEYS).any():raise RuntimeError('Duplicate original RealTrace 8-key, no Cartesian join permitted')
 output.parent.mkdir(parents=True,exist_ok=True)
 with output.open('wb') as f:
  with gzip.GzipFile(fileobj=f,mode='wb',filename='',mtime=0,compresslevel=6) as g:
   sub.to_csv(g,index=False,float_format='%.17g')
 receipt_dict={'status':'AUTHOR_SOURCE_FEATHER_VERIFIED','member_name':EXPECTED_MEMBER_NAME,
    'member_original_size_bytes':EXPECTED_MEMBER_SIZE,'member_zip_crc32':member_crc,
    'member_sha256':EXPECTED_MEMBER_SHA,'original_zip_sha256':original_zip_sha,
    'selected_rows':len(sub),'selected_gzip_csv_bytes':output.stat().st_size,
    'selected_gzip_csv_sha256':sha(output),
    'local_source_verified':True,'source_kind':'original author RealTrace stress Feather'}
 receipt.write_text(json.dumps(receipt_dict,sort_keys=True,indent=2)+'\n')
 return receipt_dict

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--feather',required=True);p.add_argument('--output',required=True)
 p.add_argument('--receipt',required=True);p.add_argument('--original-zip-sha',required=True)
 a=p.parse_args();print(json.dumps(extract(a.feather,a.output,a.receipt,a.original_zip_sha),indent=2))
