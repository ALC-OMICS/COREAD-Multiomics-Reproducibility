#!/usr/bin/env python3
from __future__ import annotations
import hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
manifest=ROOT/'MANIFEST_SHA256.tsv'

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
 return h.hexdigest()

bad=[]
for i,line in enumerate(manifest.read_text().splitlines()):
 if i==0 or not line.strip(): continue
 rel,size,digest=line.split('\t')
 p=ROOT/rel
 if not p.exists() or p.stat().st_size!=int(size) or sha(p)!=digest: bad.append(rel)
if bad: raise SystemExit('Manifest mismatch: '+', '.join(bad))
print('MANIFEST: PASS')
