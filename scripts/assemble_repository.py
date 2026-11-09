#!/usr/bin/env python3
"""Extract the complete COREAD ZIP or assemble all smaller ZIP parts, then verify every file."""
from __future__ import annotations
import argparse,csv,hashlib,json,shutil,sys,zipfile
from pathlib import Path,PurePosixPath
PREFIX='coread-multiaxis-reproducibility'
FULL='COREAD_GitHub_COMPLETO_Convergent9_v0_2_1.zip'
PATTERN='COREAD_GitHub_v0_2_1_parte_*.zip'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--downloads',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--output-dir',type=Path,help='Parent directory for the restored repository; default: downloads directory')
    a=ap.parse_args();downloads=a.downloads.resolve();out=(a.output_dir or downloads).resolve()
    parts=sorted(downloads.glob(PATTERN));archives=[downloads/FULL]if (downloads/FULL).is_file()else parts
    if not archives:raise SystemExit('No se encuentra el ZIP completo ni las partes. Coloca este script y los ZIP en la misma carpeta.')
    package_manifest=downloads/'DOWNLOADS_SHA256.tsv'
    if package_manifest.is_file():
        with package_manifest.open(newline='')as f:metadata={r['file']:r for r in csv.DictReader(f,delimiter='\t')}
        if len(archives)>1 or archives[0].name!=FULL:
            expected={n for n in metadata if n.startswith('COREAD_GitHub_v0_2_1_parte_')and n.endswith('.zip')}
            missing=expected-{p.name for p in archives}
            if missing:raise SystemExit('Faltan partes: '+', '.join(sorted(missing)))
        for p in archives:
            if p.name in metadata and sha(p)!=metadata[p.name]['sha256']:raise SystemExit('ZIP incompleto o alterado: '+p.name)
    target=out/PREFIX
    if target.exists()and any(target.iterdir()):raise SystemExit(f'La carpeta {target} ya contiene archivos. Elige --output-dir con otra carpeta.')
    out.mkdir(parents=True,exist_ok=True);seen=set()
    for archive in archives:
        print('Extrayendo',archive.name,flush=True)
        with zipfile.ZipFile(archive)as z:
            bad=z.testzip()
            if bad:raise SystemExit('Error CRC en '+archive.name+': '+bad)
            for info in z.infolist():
                rel=PurePosixPath(info.filename)
                if rel.is_absolute()or '..'in rel.parts or not rel.parts or rel.parts[0]!=PREFIX:
                    raise SystemExit('Ruta inesperada en ZIP: '+info.filename)
                if info.is_dir():continue
                if info.filename in seen:raise SystemExit('Archivo duplicado entre partes: '+info.filename)
                seen.add(info.filename);p=out.joinpath(*rel.parts);p.parent.mkdir(parents=True,exist_ok=True)
                with z.open(info)as src,p.open('wb')as dst:shutil.copyfileobj(src,dst)
    manifest=target/'MANIFEST_SHA256.tsv'
    if not manifest.exists():raise SystemExit('Falta el manifiesto: probablemente no descargaste todas las partes.')
    with manifest.open(newline='')as f:records=list(csv.DictReader(f,delimiter='\t'))
    errors=[]
    for r in records:
        p=target/r['path']
        if not p.is_file():errors.append('Falta '+r['path'])
        elif p.stat().st_size!=int(r['bytes'])or sha(p)!=r['sha256']:errors.append('No coincide '+r['path'])
    if errors:raise SystemExit('\n'.join(errors)+'\nDescarga todas las partes y vuelve a extraer en una carpeta nueva.')
    print(f'CORRECTO: {len(records)} archivos verificados.\nRepositorio completo: {target}')
if __name__=='__main__':main()
