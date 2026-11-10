"""Minimal bounds-checked Feather V1 metadata/column inspector.
Uses publicly documented Feather V1 flatbuffer schema; no arbitrary code from file.
Not a general-purpose Feather reader. Supports plain numeric, UTF8, binary, dictionary;
metadata and bounds checks fail closed. Intended for source verification only.
"""
import os,struct, numpy as np
from pathlib import Path
FMTS={0:('BOOL','u1'),1:('INT8','i1'),2:('INT16','<i2'),3:('INT32','<i4'),4:('INT64','<i8'),5:('UINT8','u1'),6:('UINT16','<u2'),7:('UINT32','<u4'),8:('UINT64','<u8'),9:('FLOAT','<f4'),10:('DOUBLE','<f8'),11:('UTF8',None),12:('BINARY',None),13:('CATEGORY',None),14:('TIMESTAMP','<i8'),15:('DATE','<i4'),16:('TIME','<i8'),17:('LARGE_UTF8',None),18:('LARGE_BINARY',None)}
class BoundsError(ValueError):pass
class FB:
 def __init__(self,b):self.b=memoryview(b);self.n=len(self.b)
 def field(self,p,i):
  self.bounds(p,4);offset=self.i32(p);v=p-offset;self.bounds(v,4);size=self.u16(v);self.bounds(v,size)
  if size < 4 or 4+i*2+2>size:return None
  f=self.u16(v+4+i*2)
  return p+f if f else None
 def bounds(self,p,size):
  if p<0 or size<0 or p+size>self.n:raise BoundsError(f'out of bounds {p=} {size=} {self.n=}')
 def unpack(self,fmt,p):
  self.bounds(p,struct.calcsize(fmt));return struct.unpack_from(fmt,self.b,p)[0]
 def i32(self,p):return self.unpack('<i',p)
 def u32(self,p):return self.unpack('<I',p)
 def i64(self,p):return self.unpack('<q',p)
 def u16(self,p):return self.unpack('<H',p)
 def u8(self,p):return self.unpack('B',p)
 def indirect(self,p):self.bounds(p,4);new=p+self.u32(p);self.bounds(new,4);return new
 def get(self,table,slot,typ,default=None):
  q=self.field(table,slot)
  if q is None:return default
  return {'b':self.u8,'u16':self.u16,'u32':self.u32,'i32':self.i32,'i64':self.i64,'link':self.indirect}[typ](q)
 def string(self,table,slot):
  pos=self.get(table,slot,'link')
  if pos is None:return None
  n=self.u32(pos);self.bounds(pos+4,n)
  return bytes(self.b[pos+4:pos+4+n]).decode('utf8','replace')
 def vector(self,table,slot):
  p=self.get(table,slot,'link')
  if p is None:return []
  n=self.u32(p);self.bounds(p+4,n*4)
  return [self.indirect(p+4+i*4) for i in range(n)]

class FeatherV1:
 def __init__(self,path):
  self.path=Path(path);self.size=self.path.stat().st_size
  with open(path,'rb') as f:
   if f.read(4)!=b'FEA1':raise ValueError('Not Feather V1 magic')
   f.seek(-8,2);end=f.read(8)
   if end[-4:]!=b'FEA1':raise ValueError('Bad footer magic')
   mlen=struct.unpack_from('<I',end,0)[0]
   if mlen>10_000_000 or mlen>self.size-12:raise ValueError('Unreasonable metadata length')
   self.metadata_start=self.size-8-mlen
   f.seek(self.metadata_start);self.metadata=f.read(mlen)
  fb=FB(self.metadata);root=fb.indirect(0)
  self.num_rows=fb.get(root,1,'i64');self.version=fb.get(root,3,'i32')
  self.description=fb.string(root,0)
  self.columns=[]
  for c in fb.vector(root,2):
   name=fb.string(c,0);a=fb.get(c,1,'link')
   if a is None:raise ValueError('column without values')
   dtype=fb.get(a,0,'b',0);encoding=fb.get(a,1,'b',0)
   offset=fb.get(a,2,'i64',0);length=fb.get(a,3,'i64',0)
   null_count=fb.get(a,4,'i64',0);total_bytes=fb.get(a,5,'i64',0)
   metadata_type=fb.get(c,2,'b',0);meta=fb.get(c,3,'link')
   if length!=self.num_rows:raise ValueError(f'column {name} length mismatch {length} != {self.num_rows}')
   if offset<4 or total_bytes<0 or offset+total_bytes>self.metadata_start:raise BoundsError(f'column data out-of-range {name} offset={offset} nbytes={total_bytes}')
   self.columns.append(dict(name=name,type=dtype,type_name=FMTS.get(dtype,('UNKNOWN',None))[0],encoding=encoding,offset=offset,length=length,null_count=null_count,total_bytes=total_bytes,metadata_type=metadata_type,meta=meta))
  self.cols={x['name']:x for x in self.columns}
 def category_levels(self,name):
  c=self.cols[name]
  if c['metadata_type']!=1 or c['meta'] is None:raise ValueError('not a dictionary encoded categorical column: '+name)
  fb=FB(self.metadata);array=fb.get(c['meta'],0,'link')
  typ=fb.get(array,0,'b',0);length=fb.get(array,3,'i64',0)
  off=fb.get(array,2,'i64',0);nbytes=fb.get(array,5,'i64',0)
  if typ!=11 or length<0 or length>100000 or off<4 or off+nbytes>self.metadata_start:raise BoundsError('bad categorical dictionary data')
  offsets=np.memmap(self.path,dtype='<i4',mode='r',offset=off,shape=(length+1,))
  payload_start=((off+(length+1)*4+7)//8)*8
  if payload_start+int(offsets[-1])>off+nbytes:raise BoundsError('bad categorical strings')
  with open(self.path,'rb') as f:
   f.seek(payload_start);payload=f.read(int(offsets[-1]))
  return [payload[int(offsets[i]):int(offsets[i+1])].decode('utf8') for i in range(length)]
 def col(self,name,indices=None):
  c=self.cols[name]
  if c['encoding']!=0:raise NotImplementedError(f'encoded column {name}, {c["encoding"]}')
  typename,dtype=FMTS[c['type']]
  n=c['length'];o=c['offset'];null=c['null_count'];data_start=o
  if null:
   data_start=o+((n+7)//8+7)//8*8
  if typename=='BOOL':
   packed=np.memmap(self.path,dtype='u1',mode='r',offset=data_start,shape=((n+7)//8,))
   if data_start+((n+7)//8)>o+c['total_bytes']:raise BoundsError('boolean bytes beyond range '+name)
   idx=np.arange(n,dtype=np.int64) if indices is None else np.asarray(indices,dtype=np.int64)
   if (idx<0).any() or (idx>=n).any():raise IndexError(name)
   return ((packed[idx>>3]>>(idx&7))&1).astype(bool)
  if dtype is not None:
   sz=np.dtype(dtype).itemsize
   if data_start+n*sz>o+c['total_bytes']:raise BoundsError('numeric bytes beyond range '+name)
   a=np.memmap(self.path,dtype=dtype,mode='r',offset=data_start,shape=(n,))
   return a if indices is None else a[indices]
  if typename in ('UTF8','BINARY','LARGE_UTF8','LARGE_BINARY'):
   index_type='<i8' if typename.startswith('LARGE_') else '<i4'
   width=np.dtype(index_type).itemsize
   offsets=np.memmap(self.path,dtype=index_type,mode='r',offset=data_start,shape=(n+1,))
   # Writer aligns the variable-length value-data buffer to 8 bytes.
   string_start=((data_start+(n+1)*width+7)//8)*8
   if string_start+int(offsets[-1])>o+c['total_bytes']:raise BoundsError('string bytes beyond column')
   with open(self.path,'rb') as f:
    if indices is None:indices=np.arange(n)
    out=[]
    for i in indices:
     i=int(i)
     lo=int(offsets[i]);hi=int(offsets[i+1]);f.seek(string_start+lo);z=f.read(hi-lo)
     out.append(z.decode('utf8','replace') if typename.endswith('UTF8') else z)
    return out
  raise NotImplementedError((typename,name))

if __name__=='__main__':
 import sys,json
 x=FeatherV1(sys.argv[1]);print('rows',x.num_rows,'version',x.version,'metadata_len',len(x.metadata),'description',x.description)
 for c in x.columns:
  print(c)
  try:print('  first3',x.col(c['name'],np.arange(3)))
  except Exception as e:print('  ERROR',type(e).__name__,str(e))
