#!/usr/bin/env python3
"""Package a verified repository as one complete ZIP and smaller normal ZIP parts."""
from __future__ import annotations
import argparse,csv,hashlib,io,json,os,shutil,subprocess,sys,tempfile,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PREFIX='coread-multiaxis-reproducibility'
FULL='COREAD_GitHub_COMPLETO_Convergent9_v0_2_1.zip'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()
def atomic(p,data):
    fd,tmp=tempfile.mkstemp(dir=p.parent,prefix='.'+p.name+'.')
    with os.fdopen(fd,'wb')as f:f.write(data);f.flush();os.fsync(f.fileno())
    os.replace(tmp,p)
def makezip(path,files):
    fd,tmp=tempfile.mkstemp(dir=path.parent,prefix='.'+path.name+'.');os.close(fd)
    with zipfile.ZipFile(tmp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True)as z:
        for rel in files:
            info=zipfile.ZipInfo(f'{PREFIX}/{rel}',date_time=(2026,10,31,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o100644<<16
            with (ROOT/rel).open('rb')as src,z.open(info,'w',force_zip64=True)as dst:shutil.copyfileobj(src,dst)
    with open(tmp,'rb')as f:os.fsync(f.fileno())
    os.replace(tmp,path)
def verify_zip(path,expected):
    with zipfile.ZipFile(path)as z:
        if z.testzip():raise AssertionError('ZIP CRC failure: '+str(path))
        for info in z.infolist():
            rel=info.filename.split('/',1)[1];h=hashlib.sha256()
            with z.open(info)as f:
                for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
            if h.hexdigest()!=expected[rel]:raise AssertionError('ZIP member mismatch: '+rel)
        return z.infolist()
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,default=ROOT.parent/'COREAD_GitHub_delivery_v0_2_1');ap.add_argument('--part-mib',type=int,default=22);a=ap.parse_args()
    subprocess.run([sys.executable,str(ROOT/'scripts/verify_repository.py')],check=True)
    out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    with (ROOT/'MANIFEST_SHA256.tsv').open(newline='')as f:records=list(csv.DictReader(f,delimiter='\t'))
    expected={r['path']:r['sha256']for r in records};expected['MANIFEST_SHA256.tsv']=sha(ROOT/'MANIFEST_SHA256.tsv')
    files=['MANIFEST_SHA256.tsv']+sorted(r['path']for r in records)
    full=out/FULL;makezip(full,files);infos=verify_zip(full,expected)
    groups=[];current=[];size=0;limit=a.part_mib*1024**2
    for info in infos:
        rel=info.filename.split('/',1)[1];weight=info.compress_size+len(info.filename.encode())*2+512
        if current and size+weight>limit:groups.append(current);current=[];size=0
        current.append(rel);size+=weight
    if current:groups.append(current)
    delivered=[full]
    for i,group in enumerate(groups,1):
        p=out/f'COREAD_GitHub_v0_2_1_parte_{i:02d}.zip';makezip(p,group);verify_zip(p,expected);delivered.append(p)
    helper=out/'extraer_repositorio.py';shutil.copyfile(ROOT/'scripts/assemble_repository.py',helper);delivered.append(helper)
    text=f'''REPOSITORIO COREAD COMPLETO — versión 0.2.1\n\nOPCIÓN 1: descarga {FULL} y descomprímelo.\n\nOPCIÓN 2: si prefieres archivos más pequeños, descarga las {len(groups)} partes ZIP,\nextraer_repositorio.py y DOWNLOADS_SHA256.tsv en una misma carpeta.\nAbre una terminal en esa carpeta y ejecuta:\n\n    python3 extraer_repositorio.py\n\nEn Windows puede ser:\n\n    py -3 extraer_repositorio.py\n\nNo es necesario descargar el ZIP completo si tienes todas las partes.\nLas partes son ZIP normales que conservan la misma carpeta raíz; el script\nlas reúne y comprueba todos los archivos contra su huella SHA-256.\nNo requiere instalar paquetes adicionales de Python.\nEl script también puede verificar y extraer el ZIP completo.\nUsa --output-dir OTRA_CARPETA si el destino ya existe.\n\nLa carpeta resultante contiene los siete módulos, datos procesados, figuras,\nSource Data para31paneles, ST_1–ST_14, documentación y validaciones.\nConsulta LEEME.md y README.md dentro del repositorio.\n'''
    instructions=out/'LEEME_DESCARGA.txt';atomic(instructions,text.encode('utf-8'));delivered.append(instructions)
    inventory=[]
    for p in delivered:inventory.append({'file':p.name,'bytes':p.stat().st_size,'sha256':sha(p)})
    b=io.StringIO();w=csv.DictWriter(b,fieldnames=['sha256','bytes','file'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(inventory)
    atomic(out/'DOWNLOADS_SHA256.tsv',b.getvalue().encode())
    report={'release':'0.2.1','root_folder':PREFIX,'repository_files':len(expected),'all_zip_members_hash_verified':True,'parts':len(groups),'files':inventory}
    atomic(out/'DELIVERY_MANIFEST.json',(json.dumps(report,indent=2)+'\n').encode());print(json.dumps(report,indent=2))
if __name__=='__main__':main()
