#!/usr/bin/env python3
"""Download the externally hosted GTDB/SILVA paired QIIME2 QZV and verify its hash."""
import hashlib, pathlib, urllib.request, sys
from reproduce_independent_validation import SOURCE_URL, SOURCE_SHA256
out=pathlib.Path(__file__).parent/'input'/'asv-seqs-gtdb-silva-fixed.qzv'
out.parent.mkdir(parents=True,exist_ok=True)
tmp=out.with_suffix('.qzv.part')
h=hashlib.sha256()
try:
    with urllib.request.urlopen(SOURCE_URL,timeout=90) as response, tmp.open('wb') as dst:
        while True:
            chunk=response.read(1<<20)
            if not chunk:break
            h.update(chunk);dst.write(chunk)
    if h.hexdigest()!=SOURCE_SHA256:
        tmp.unlink(missing_ok=True)
        raise SystemExit('ERROR: source SHA-256 mismatch; source may have changed')
    tmp.replace(out)
    print('Verified external QZV:',out,'SHA-256:',h.hexdigest())
except Exception:
    tmp.unlink(missing_ok=True)
    raise
