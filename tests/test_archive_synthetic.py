"""Synthetic source-integrity tests for original archive extractor: NOT an authentic biology rerun."""
import hashlib, tempfile, unittest, zipfile
from pathlib import Path
from unittest import mock
from sys import path
path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import original_zenodo_archive as o

class SourceArchiveIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.t=tempfile.TemporaryDirectory(); self.root=Path(self.t.name)
        self.contents={'stress':b'FEA-original-test', 'realtrace':b'RT-original-test', 'rmd':b'RMD-original-test'}
        self.members={k:('source/'+k+('.Rmd' if k=='rmd' else '.feather'),len(b),hashlib.sha256(b).hexdigest()) for k,b in self.contents.items()}
        self.source=self.root/'mock.zip'
        with zipfile.ZipFile(self.source,'w',zipfile.ZIP_DEFLATED) as z:
            for k,(name,_,_) in self.members.items():z.writestr(name,self.contents[k])
        self.patches=[mock.patch.object(o,'SOURCE_SIZE',self.source.stat().st_size),mock.patch.object(o,'SOURCE_SHA',o.sha(self.source)),mock.patch.object(o,'MEMBERS',self.members)]
        for p in self.patches:p.start()
    def tearDown(self):
        for p in reversed(self.patches):p.stop()
        self.t.cleanup()
    def test_tiny_source_roundtrip_and_receipt(self):
        d=o.extract(self.source,self.root/'extracted')
        self.assertEqual(set(d)-{'source_zip_sha256'},set(self.contents))
        for k,v in d.items():
            if k=='source_zip_sha256':continue
            self.assertEqual(Path(v['file']).read_bytes(),self.contents[k])
        self.assertTrue((self.root/'extracted/EXTRACTED_AUTHOR_SOURCE_MANIFEST.json').exists())
    def test_corrupt_zip_hash_rejected_before_extract(self):
        with self.source.open('ab') as f:f.write(b'corrupt')
        with self.assertRaises(RuntimeError):o.extract(self.source,self.root/'rejected')
        self.assertFalse((self.root/'rejected').exists())
    def test_member_size_mismatch_rejected(self):
        k='realtrace';m=dict(self.members);name,nbytes,sha=m[k];m[k]=(name,nbytes+1,sha)
        with mock.patch.object(o,'MEMBERS',m):
            with self.assertRaises(RuntimeError):o.extract(self.source,self.root/'rejected')
    def test_member_sha_mismatch_rejected(self):
        m=dict(self.members);name,nbytes,sha=m['stress'];m['stress']=(name,nbytes,'0'*64)
        with mock.patch.object(o,'MEMBERS',m):
            with self.assertRaises(RuntimeError):o.extract(self.source,self.root/'bad_member')
    def test_preexisting_destination_rejected(self):
        d=self.root/'existing';d.mkdir()
        with self.assertRaises(FileExistsError):o.extract(self.source,d)
    def test_missing_member_rejected(self):
        m=dict(self.members);m['stress']=('source/absent.feather',len(self.contents['stress']),hashlib.sha256(self.contents['stress']).hexdigest())
        with mock.patch.object(o,'MEMBERS',m):
            with self.assertRaises(KeyError):o.extract(self.source,self.root/'bad_missing')

if __name__=='__main__':unittest.main()
