#!/usr/bin/env python3
"""Explicit maintenance: refresh current manifests after reviewed code changes.

Historical *_original manifest records and immutable pseudobulk input checksums
are not rewritten. This command does not recompute scientific reference values.
"""
from pathlib import Path
import csv,hashlib,io,json
from verify_repository import ROOT,files,scope_errors,write_manifest,atomic_text

def main():
    errors=scope_errors(ROOT)
    if errors:raise SystemExit('\n'.join(errors))
    for module in ['01_dependency','02_genomics']:
        root=ROOT/'analysis'/module
        active=root/'provenance/sha256_manifest.json'
        # Active checkpoint manifests protect analytical inputs, frozen numeric
        # outputs and executable code. Historical provenance is covered globally.
        kept={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
              for p in files(root)
              if p != active and p.name != 'MANIFEST_SHA256.tsv' and
              p.relative_to(root).parts[0] != 'provenance'}
        atomic_text(active,json.dumps(kept,indent=2)+'\n')
    for module in ['02_genomics','03_single_cell']:
        root=ROOT/'analysis'/module
        manifest=root/'MANIFEST_SHA256.tsv'
        buffer=io.StringIO(newline='');w=csv.writer(buffer,delimiter='\t',lineterminator='\n')
        w.writerow(['relative_path','bytes','sha256'])
        for p in files(root):
            if p!=manifest:
                w.writerow([p.relative_to(root).as_posix(),p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest()])
        atomic_text(manifest,buffer.getvalue())
    write_manifest(ROOT)
    print('Refreshed active manifests; frozen numerical references were not changed.')
if __name__=='__main__':main()
