#!/usr/bin/env python3
"""Check the code/data-only release, file hashes, source contracts and gene sets.

Refreshing a manifest is a deliberate maintenance action; ordinary verification
never updates an expected checksum or treats a missing measurement as zero.
"""
from __future__ import annotations
import argparse
import ast
import csv
import hashlib
import io
import json
import os
import sys
import tempfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
IGNORED = {'.git', '.venv', '__pycache__', 'run_outputs', '_runs', 'build',
           'default_runs', 'full_runs', '.pytest_cache'}
EDITORIAL_DIRS = {'figures', 'supplementary_tables', 'manuscript', 'submission'}
EDITORIAL_EXTENSIONS = {'.png', '.svg', '.jpg', '.jpeg', '.tif', '.tiff', '.pdf',
                        '.doc', '.docx', '.ppt', '.pptx', '.xlsx', '.xls'}
# These five archived cBioPortal exports are original measurement inputs, not
# manuscript supplementary workbooks. No other spreadsheet attachment is allowed.
SOURCE_SPREADSHEETS = {
    'analysis/02_genomics/input/PATIENT_DATA_oncoprint-COLON-TCGA-GDC.xlsx',
    'analysis/02_genomics/input/PATIENT_DATA_oncoprint-CPTAC2.xlsx',
    'analysis/02_genomics/input/PATIENT_DATA_oncoprint-PANCANCER-PARTE -2.xlsx',
    'analysis/02_genomics/input/PATIENT_DATA_oncoprint-PANCANCER-PARTE 1-.xlsx',
    'analysis/02_genomics/input/PATIENT_DATA_oncoprint-RECTAL-TCGA-GDC.xlsx',
}

MODULES = ['01_dependency', '02_genomics', '03_single_cell',
           '04_discovery_spatial', '05_external_visium_hd',
           '06_crc16_niches', '07_translational']
REQUIRED = [
    'LICENSE', 'LICENSE_NOTICE.md', 'CITATION.cff', 'README.md', 'requirements.txt',
    'run_all.py', 'config/candidate_sets.json', 'source_data/panel_index.csv',
    'analysis/01_dependency/data/inputs/Texto pegado(1).txt',
    'analysis/01_dependency/data/organoids/query_gene_effect_108_screens.csv',
    'analysis/01_dependency/results/csv/organoid_results_1720.csv',
    'analysis/01_dependency/results/csv/organoid_differential_context.csv',
    'analysis/02_genomics/input/union_candidates_466.txt',
    'analysis/02_genomics/results/COREAD_GDC_calls.npz',
    'analysis/03_single_cell/input/pseudobulk/pseudobulk_counts.npz',
    'analysis/03_single_cell/input/pseudobulk/manifest_sha256.json',
    'analysis/04_discovery_spatial/reference/discovery_spatial_expected.csv',
    'analysis/05_external_visium_hd/data/frozen_discovery_annotations.csv',
    'analysis/06_crc16_niches/input/CRC16_host_counts.npz',
    'analysis/06_crc16_niches/input/CRC16_frozen_niches.csv',
] + ['analysis/'+m+'/run.py' for m in MODULES]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def files(root):
    for path in sorted(root.rglob('*')):
        relative = path.relative_to(root)
        if any(x in IGNORED for x in relative.parts):
            continue
        if path.is_file() or path.is_symlink():
            if not path.name.endswith(('.pyc', '.tmp')):
                yield path


def safe_path(root, value):
    relative = PurePosixPath(value)
    if (not value or relative.is_absolute() or '..' in relative.parts or
            '.git' in relative.parts or '\\' in value or '\x00' in value):
        raise ValueError('Unsafe path: '+value)
    path = root.joinpath(*relative.parts)
    if path.is_symlink() or not path.resolve().is_relative_to(root.resolve()):
        raise ValueError('Symlink or escaped path: '+value)
    return path


def scope_errors(root):
    errors = []
    for path in files(root):
        relative = path.relative_to(root)
        if path.is_symlink():
            errors.append('Symlink not allowed in release: '+relative.as_posix())
        if (EDITORIAL_DIRS.intersection(relative.parts) or
                (path.suffix.lower() in EDITORIAL_EXTENSIONS and relative.as_posix() not in SOURCE_SPREADSHEETS) or
                path.name == '.DS_Store' or
                path.name.lower().startswith(('supplementary_table', 'methods_original', 'methods_archived'))):
            errors.append('Editorial/temporary file excluded by release policy: '+relative.as_posix())
    return errors


def atomic_text(path, text):
    fd, name = tempfile.mkstemp(dir=path.parent, prefix='.'+path.name+'.')
    try:
        with os.fdopen(fd, 'w', encoding='utf-8', newline='') as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(name, path)
    finally:
        if Path(name).exists():
            Path(name).unlink()


def write_manifest(root):
    errors = scope_errors(root)
    errors.extend('Missing required input: '+name for name in REQUIRED if not (root/name).is_file())
    if errors:
        raise ValueError('\n'.join(errors))
    manifest = root/'MANIFEST_SHA256.tsv'
    buffer = io.StringIO(newline='')
    writer = csv.writer(buffer, delimiter='\t', lineterminator='\n')
    writer.writerow(['sha256', 'bytes', 'path'])
    for path in files(root):
        if path != manifest:
            writer.writerow([digest(path), path.stat().st_size, path.relative_to(root).as_posix()])
    atomic_text(manifest, buffer.getvalue())


def verify(root=ROOT):
    root = root.resolve()
    errors = scope_errors(root)
    errors.extend('Missing required input: '+name for name in REQUIRED if not (root/name).is_file())
    manifest = root/'MANIFEST_SHA256.tsv'
    records = []
    if not manifest.is_file():
        errors.append('Missing MANIFEST_SHA256.tsv')
    else:
        with manifest.open(newline='', encoding='utf-8') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if reader.fieldnames != ['sha256','bytes','path']:
                errors.append('Invalid manifest schema')
            else:
                records = list(reader)
    expected = set()
    for row in records:
        name = row['path']
        if name in expected:
            errors.append('Duplicate manifest path: '+name)
        expected.add(name)
        try:
            path = safe_path(root, name)
            if not path.is_file():
                errors.append('Missing: '+name)
            elif path.stat().st_size != int(row['bytes']) or digest(path) != row['sha256']:
                errors.append('Changed: '+name)
        except (ValueError, OSError) as exc:
            errors.append(str(exc))
    actual = {p.relative_to(root).as_posix() for p in files(root) if p != manifest}
    errors.extend('Unmanifested: '+name for name in sorted(actual-expected))
    try:
        sets = json.loads((root/'config/candidate_sets.json').read_text())
        for key,n in [('Core47',47),('Spatial43',43),('Convergent9',9),('Extended13',13)]:
            if len(sets[key]) != n or len(set(sets[key])) != n:
                errors.append(key+': incorrect size or duplicate genes')
        c9 = {'BCL2L1','CTNNB1','FASN','HMGB1','IRS2','MYC','RNF6','SOX9','YAP1'}
        if set(sets['Convergent9']) != c9:
            errors.append('Convergent9 membership mismatch')
        for small,big in [('Convergent9','Spatial43'),('Spatial43','Core47'),('Convergent9','Extended13')]:
            if not set(sets[small]) <= set(sets[big]):
                errors.append(small+' not a subset of '+big)
        with (root/'source_data/panel_index.csv').open(newline='',encoding='utf-8') as handle:
            panels = list(csv.DictReader(handle))
        if len(panels) != 31 or len({(p['figure'],p['panel']) for p in panels}) != 31:
            errors.append('Expected 31 unique numeric panel-source records')
        for panel in panels:
            for name in panel['source_files'].split(';'):
                name = name.strip()
                if name and not safe_path(root,name).is_file():
                    errors.append('Missing numeric panel source: '+name)
    except (ValueError, KeyError, OSError) as exc:
        errors.append('Invalid source contract: '+str(exc))
    scripts = 0
    for path in files(root):
        if path.suffix == '.py':
            scripts += 1
            try:
                ast.parse(path.read_text(encoding='utf-8-sig'), filename=path.relative_to(root).as_posix())
            except (SyntaxError, UnicodeError) as exc:
                errors.append('Python syntax/encoding: '+str(exc))
    return {'passed':not errors, 'scope':'code-and-numeric-data-only',
            'files_checked':len(records), 'python_scripts_checked':scripts, 'errors':errors}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-manifest',action='store_true',help='Maintenance only: explicitly refresh checksum inventory')
    args = parser.parse_args()
    if args.write_manifest:
        write_manifest(ROOT)
    report = verify()
    print(json.dumps(report,indent=2))
    return 0 if report['passed'] else 1

if __name__ == '__main__':
    sys.exit(main())
