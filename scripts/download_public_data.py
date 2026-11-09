#!/usr/bin/env python3
"""Download optional public source objects and accept them only after SHA-256 verification."""
import argparse,hashlib,json,os,tempfile,urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(path):
    h=hashlib.sha256()
    with path.open('rb')as f:
        for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--groups',nargs='+',choices=['single_cell','external','crc16'],default=['single_cell','external','crc16'])
    ap.add_argument('--output-dir',type=Path,default=ROOT/'data/raw');a=ap.parse_args()
    manifest=json.loads((ROOT/'provenance/public_sources.json').read_text())
    for row in manifest['files']:
        if row['group']not in a.groups:continue
        p=a.output_dir/row['path'];p.parent.mkdir(parents=True,exist_ok=True)
        if p.is_file()and sha(p)==row['sha256']:
            print('Verified existing',p,flush=True);continue
        failure=None
        for attempt in range(4):
            fd,tmp=tempfile.mkstemp(dir=p.parent,prefix='.'+p.name+'.')
            try:
                url=row['url']+(''if attempt==0 else '?download=1&attempt='+str(attempt))
                req=urllib.request.Request(url,headers={'User-Agent':'COREAD-reproducibility/0.2.0','Accept-Encoding':'identity'})
                with os.fdopen(fd,'wb')as f,urllib.request.urlopen(req,timeout=120)as response:
                    for block in iter(lambda:response.read(8*1024**2),b''):f.write(block)
                    f.flush();os.fsync(f.fileno())
                if Path(tmp).stat().st_size!=row['bytes'] or sha(Path(tmp))!=row['sha256']:
                    raise ValueError('Downloaded bytes do not match the recorded source identity')
                os.replace(tmp,p);print('Downloaded and verified',p,flush=True);failure=None;break
            except Exception as exc:
                failure=exc
                try:os.close(fd)
                except OSError:pass
                Path(tmp).unlink(missing_ok=True)
        if failure:raise RuntimeError(f'Could not verify {row["path"]}: {failure}')
if __name__=='__main__':main()
