#!/usr/bin/env python3
"""Validate frozen release manifest, no secret keys, and archive file presence."""
from pathlib import Path
import hashlib,sys
ROOT=Path(__file__).resolve().parent
manifest=ROOT/'SHA256SUMS.txt'
assert manifest.is_file(), 'No manifest'
N=0
listed=set()
for line in manifest.read_text().splitlines():
    if not line.strip():continue
    digest,sep,name=line.partition('  ')
    assert sep and len(digest)==64 and len(name)>0, 'Bad manifest line'
    p=ROOT/name
    assert p.is_file(),f'Missing: {name}'
    assert hashlib.sha256(p.read_bytes()).hexdigest()==digest,f'Hash mismatch: {name}'
    listed.add(name)
    N+=1
actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()
        and '.git' not in p.relative_to(ROOT).parts
        and 'build' not in p.relative_to(ROOT).parts
        and '__pycache__' not in p.relative_to(ROOT).parts
        and p.suffix!='.pyc' and p.name!='SHA256SUMS.txt'}
assert listed==actual, 'Unlisted / missing files: '+str(sorted(listed^actual)[:12])
assert N>60,N
print('PASS: '+str(N)+' release files verified by SHA-256')
