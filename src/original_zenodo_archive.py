#!/usr/bin/env python3
"""Fail-closed streaming extraction from exactly verified official Zenodo ZIP.

Original data and R Markdown remain OUTSIDE the public code repository.
"""
from pathlib import Path
import argparse,hashlib,json,shutil,zipfile
SOURCE_SHA='0697803959ffa64b8b2bc849d716bc841bcb5afe8a42e92b53f397f8c7d7494b'
SOURCE_SIZE=842127932
MEMBERS={
 'stress':('Gervais2026_data/Exported_data/myframes_stress_exp.feather',482147704,'52c9cf4c7607ecd32a621bb410c967ff5a198d3cdf95f37323008319bb9272ac'),
 'realtrace':('Gervais2026_data/Exported_data/myframes_RealTrace_stress_exp.feather',264629352,'17bdb0bc9e45f0aff5435f0a168d87850a3a29370ccf53449201fbb63f18b217'),
 'rmd':('Gervais2026_data/Gervais_et_al_analysis.Rmd',207248,'7f915a34db6c3326a0be833bc1c04e57127d90684ce281c9fe4b0cd65a48de70')}
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for chunk in iter(lambda:f.read(4<<20),b''):h.update(chunk)
 return h.hexdigest()
def extract(archive,destination):
 archive=Path(archive);dest=Path(destination)
 if dest.exists():raise FileExistsError('Original source extraction directory already exists')
 if archive.stat().st_size!=SOURCE_SIZE or sha(archive)!=SOURCE_SHA:
  raise RuntimeError('Original Zenodo Gervais2026_data.zip byte-size or SHA mismatch')
 with zipfile.ZipFile(archive) as z:
  # z.open(...).read() checks CRC for each selected member when fully consumed.
  entries={k:z.getinfo(name) for k,(name,_,_) in MEMBERS.items()}
  for key,member in entries.items():
   expected_name,expected_size,_=MEMBERS[key]
   if member.file_size!=expected_size:raise RuntimeError(f'Unexpected member size {expected_name}')
  dest.mkdir(parents=True)
  receipt={}
  for key,member in entries.items():
   name,nbytes,expect=MEMBERS[key];dst=dest/Path(name).name
   h=hashlib.sha256();total=0
   with z.open(member) as src,dst.open('wb') as out:
    for piece in iter(lambda:src.read(4<<20),b''):
     out.write(piece);h.update(piece);total+=len(piece)
   if total!=nbytes or h.hexdigest()!=expect:raise RuntimeError(f'Unexpected original member {key} SHA/size')
   receipt[key]={'file':str(dst),'bytes':total,'sha256':h.hexdigest(),'crc32':f'{member.CRC:08x}'}
 receipt['source_zip_sha256']=SOURCE_SHA
 (dest/'EXTRACTED_AUTHOR_SOURCE_MANIFEST.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
 return receipt
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--archive',required=True);p.add_argument('--out',required=True)
 a=p.parse_args();print(json.dumps(extract(a.archive,a.out),indent=2))
