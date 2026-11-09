#!/usr/bin/env python3
"""Reconstructed 2026-10-31: derive patient pseudobulks from the public UMI H5.

The original statistical inference script was not supplied. This implementation
uses the archived feature mapping and published cellular annotations, retaining
complete-feature library denominators. It is not an upstream sequencing pipeline.
"""
from pathlib import Path
import argparse, hashlib, io, json, os, time, zlib
import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = 'f435bb2651ff5297d0c24a99daf58850ed67ae1ed6c5ef05fad48fa3f0186670'
C9 = ['BCL2L1','CTNNB1','FASN','HMGB1','IRS2','MYC','RNF6','SOX9','YAP1']

def atomic(path, payload):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.tmp')
    with temp.open('wb') as f: f.write(payload); f.flush(); os.fsync(f.fileno())
    os.replace(temp,path)
    assert path.stat().st_size==len(payload)

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''): h.update(b)
    return h.hexdigest()

def expand_dataset(h5path, dataset, out):
    """Disk-backed decompression and unshuffle of oversized HDF5 chunks.

    Process each physical chunk separately. The public source uses shuffle then
    deflate; byte planes are restored in blocks, so RAM does not scale with the
    multi-gigabyte HDF5 chunk size. Expanded cache files are scratch artifacts.
    """
    out=Path(out); out.parent.mkdir(parents=True,exist_ok=True)
    with h5py.File(h5path,'r') as h:
        d=h[dataset]; shape=tuple(d.shape); dtype=d.dtype; expected=int(d.size*dtype.itemsize)
        if out.exists() and out.stat().st_size==expected:
            return np.memmap(out,mode='r',dtype=dtype,shape=shape)
        temp=out.with_name(out.name+'.tmp')
        filters=[d.id.get_create_plist().get_filter(i)[0] for i in range(d.id.get_create_plist().get_nfilters())]
        if not d.chunks or filters not in ([1],[2,1]):
            raise ValueError(f'Unsupported streaming layout: {dataset}, {filters}; inspect before allocating')
        with open(h5path,'rb') as src,temp.open('wb') as dst:
            for chunk_id in range(d.id.get_num_chunks()):
                info=d.id.get_chunk_info(chunk_id); assert info.filter_mask==0
                decoded=out.with_name(out.name+'.shuffle.tmp');decoder=zlib.decompressobj();written=0
                with decoded.open('wb') as expansion:
                    src.seek(info.byte_offset);remaining=info.size
                    while remaining:
                        chunk=src.read(min(2*1024**2,remaining)); remaining-=len(chunk)
                        if not chunk: raise IOError('Truncated HDF5 chunk')
                        while chunk:
                            block=decoder.decompress(chunk,8*1024**2);chunk=decoder.unconsumed_tail
                            expansion.write(block);written+=len(block)
                    block=decoder.flush();expansion.write(block);written+=len(block)
                    expansion.flush();os.fsync(expansion.fileno())
                assert decoder.eof and written%dtype.itemsize==0
                full_n=written//dtype.itemsize;valid=min(full_n,shape[0]-int(info.chunk_offset[0]))
                if filters==[2,1]:
                    byteplanes=np.memmap(decoded,mode='r',dtype=np.uint8,shape=(dtype.itemsize,full_n))
                    for first in range(0,valid,1_000_000):
                        dst.write(np.ascontiguousarray(byteplanes[:,first:min(first+1_000_000,valid)].T).tobytes())
                    del byteplanes
                else:
                    with decoded.open('rb') as raw:
                        left=valid*dtype.itemsize
                        while left:
                            block=raw.read(min(left,8*1024**2));dst.write(block);left-=len(block)
                decoded.unlink();print(f'Expanded {dataset} chunk {chunk_id+1}/{d.id.get_num_chunks()}',flush=True)
            dst.flush();os.fsync(dst.fileno())
        assert temp.stat().st_size==expected,(dataset,temp.stat().st_size,expected)
        os.replace(temp,out)
        atomic(out.with_suffix(out.suffix+'.json'),json.dumps({'dataset':dataset,'shape':shape,'dtype':str(dtype),'bytes':expected,'source_sha256':EXPECTED},indent=2).encode())
        return np.memmap(out,mode='r',dtype=dtype,shape=shape)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--h5',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--expanded-dir',type=Path,required=True)
    ap.add_argument('--c9-cell-output',type=Path);ap.add_argument('--block-cells',type=int,default=5000)
    a=ap.parse_args(); start=time.time(); a.output.mkdir(parents=True,exist_ok=True)
    sha=digest(a.h5); assert sha==EXPECTED,'Unexpected public source checksum'
    mp=pd.read_csv(ROOT/'input/gene_feature_matches_full466.csv'); genes=mp.Gene.to_numpy(str)
    with h5py.File(a.h5,'r') as h:
        shape=tuple(h['matrix/shape'][:]); assert shape==(43113,370115)
        barcodes=np.array([x.decode() for x in h['matrix/barcodes'][:]])
        names=np.array([x.decode() for x in h['matrix/features/name'][:]])
        ids=np.array([x.decode() for x in h['matrix/features/id'][:]])
        indptr=h['matrix/indptr'][:]
    assert np.array_equal(names[mp.feature_row],mp.H5_feature_name)
    assert np.array_equal(ids[mp.feature_row],mp.H5_feature_id)
    cluster=pd.read_csv(ROOT/'input/GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz')
    meta=pd.read_csv(ROOT/'input/GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz')
    ann=cluster.merge(meta,left_on='sampleID',right_on='cellID',validate='one_to_one').set_index('sampleID').loc[barcodes].reset_index()
    assert len(ann)==shape[1] and not ann.PID.isna().any()
    key=['PID','SPECIMEN_TYPE','MMRStatus','clTopLevel']
    top=ann[key].drop_duplicates().reset_index(drop=True); top.insert(0,'group_id',np.arange(len(top)))
    top_id=ann[key].merge(top,on=key,how='left',sort=False).group_id.to_numpy()
    epmask=ann.clTopLevel.eq('Epi').to_numpy();finekey=['PID','SPECIMEN_TYPE','MMRStatus','cl295v11SubFull']
    fine=ann.loc[epmask,finekey].drop_duplicates().reset_index(drop=True);fine.insert(0,'group_id',np.arange(len(fine)))
    fine_id=np.full(len(ann),-1,dtype=int);fine_id[epmask]=ann.loc[epmask,finekey].merge(fine,on=finekey,how='left',sort=False).group_id.to_numpy()
    assert len(top)==671 and len(fine)==969,(len(top),len(fine))
    feature_map=np.full(shape[0],-1,dtype=np.int32);feature_map[mp.feature_row]=np.arange(len(mp))
    top_counts=np.zeros((len(top),len(mp)),dtype=np.int64);top_positive=np.zeros_like(top_counts)
    fine_counts=np.zeros((len(fine),len(mp)),dtype=np.int64);fine_positive=np.zeros_like(fine_counts)
    libs=np.zeros(shape[1],dtype=np.int64); c9counts=np.zeros((shape[1],len(C9)),dtype=np.int32)
    c9map=np.full(len(mp),-1,dtype=np.int32)
    for i,g in enumerate(C9):c9map[np.where(genes==g)[0][0]]=i
    data=expand_dataset(a.h5,'matrix/data',a.expanded_dir/'data.bin')
    indices=expand_dataset(a.h5,'matrix/indices',a.expanded_dir/'indices.bin')
    b=io.BytesIO();np.save(b,indptr);atomic(a.expanded_dir/'indptr.npy',b.getvalue())
    for s in range(0,shape[1],a.block_cells):
        e=min(s+a.block_cells,shape[1]);p=indptr[s:e+1];lo,hi=int(p[0]),int(p[-1])
        v=np.asarray(data[lo:hi]);idx=np.asarray(indices[lo:hi]);n=np.diff(p)
        if not np.all((v>=0)&(v==np.floor(v))):raise ValueError('Source UMI entries are not non-negative integers')
        cell=np.repeat(np.arange(s,e,dtype=np.int32),n)
        libs[s:e]=np.bincount(cell-s,weights=v,minlength=e-s).astype(np.int64)
        j=feature_map[idx];ok=j>=0;ci=cell[ok];gj=j[ok];val=v[ok]
        np.add.at(top_counts,(top_id[ci],gj),val);np.add.at(top_positive,(top_id[ci],gj),val>0)
        sel=fine_id[ci]>=0
        np.add.at(fine_counts,(fine_id[ci[sel]],gj[sel]),val[sel]);np.add.at(fine_positive,(fine_id[ci[sel]],gj[sel]),val[sel]>0)
        keep=c9map[gj]>=0;c9counts[ci[keep],c9map[gj[keep]]]=val[keep]
        if s%(a.block_cells*10)==0:print(f'Aggregated {e:,}/{shape[1]:,} cells',flush=True)
    top['n_cells']=np.bincount(top_id,minlength=len(top));top['library_size']=np.bincount(top_id,weights=libs,minlength=len(top)).astype(np.int64)
    fine['n_cells']=np.bincount(fine_id[epmask],minlength=len(fine));fine['library_size']=np.bincount(fine_id[epmask],weights=libs[epmask],minlength=len(fine)).astype(np.int64)
    fine=fine.rename(columns={'cl295v11SubFull':'epiFine'})
    atomic(a.output/'top_groups.csv',top.to_csv(index=False).encode());atomic(a.output/'fine_groups.csv',fine.to_csv(index=False).encode())
    b=io.BytesIO();np.savez_compressed(b,genes=genes,top_counts=top_counts,top_positive=top_positive,fine_counts=fine_counts,fine_positive=fine_positive)
    atomic(a.output/'pseudobulk_counts.npz',b.getvalue())
    if a.c9_cell_output:
        b=io.BytesIO();np.savez_compressed(b,genes=np.array(C9),counts=c9counts,library_size=libs,barcodes=barcodes)
        atomic(a.c9_cell_output,b.getvalue())
    report={'passed':True,'implementation':'reconstructed 2026-10-31','source_sha256':sha,'features':int(shape[0]),'cells':int(shape[1]),'top_groups':len(top),'epithelial_fine_groups':len(fine),'complete_feature_library_total':int(libs.sum()),'candidate_genes':len(mp),'seconds':time.time()-start}
    atomic(a.output/'preparation_report.json',json.dumps(report,indent=2).encode())
    manifest={p.name:digest(p) for p in a.output.iterdir() if p.is_file() and p.name!='manifest_sha256.json'}
    atomic(a.output/'manifest_sha256.json',json.dumps(manifest,indent=2).encode());print(json.dumps(report,indent=2))

if __name__=='__main__':main()
