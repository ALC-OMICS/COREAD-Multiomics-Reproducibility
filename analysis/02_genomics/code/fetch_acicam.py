"""Optional network-enabled reconstruction of AC-ICAM archived summaries.
Downloads actual pinned public DataHub files and checks original LFS SHA-256.
Not needed for offline reanalysis of the provided summary inputs.
"""
from pathlib import Path
import urllib.request,hashlib,csv,json,re,math,io
from collections import defaultdict
import numpy as np
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'downloads';RAW=D/'raw_acicam';RAW.mkdir(exist_ok=True)
META=json.loads((ROOT/'provenance/acicam_public_source.json').read_text());COMMIT=META['source_commit']
BASE=f'https://media.githubusercontent.com/media/cBioPortal/datahub/{COMMIT}/public/coad_silu_2022/'
G=(ROOT/'input/union_candidates_466.txt').read_text().split();GS=set(G)

def download(name):
 p=RAW/name
 if not p.exists():
  req=urllib.request.Request(BASE+name,headers={'User-Agent':'COREAD-patient-reanalysis/1.0'})
  with urllib.request.urlopen(req,timeout=120) as r,p.open('wb') as f:
   while data:=r.read(1024*1024):f.write(data)
 expected=META['original_LFS_objects'][name]
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1024*1024):h.update(b)
 if h.hexdigest()!=expected['sha256'] or p.stat().st_size!=expected['bytes']:
  raise ValueError('Source checksum/size mismatch: '+name)
 return p

def matrix(path):
 with path.open(encoding='utf-8-sig',newline='') as f:
  reader=csv.reader(f,delimiter='\t');header=next(reader)
  ix=[j for j,x in enumerate(header) if x.startswith('SER-SILU-CC-')];samples=[header[j] for j in ix];out={}
  for row in reader:
   if row and row[0] in GS:
    if row[0] in out:raise ValueError('Duplicate symbol: '+row[0])
    vals=[]
    for j in ix:
     try:v=float(row[j]);v=v if math.isfinite(v) else np.nan
     except (ValueError,IndexError):v=np.nan
     vals.append(v)
    out[row[0]]=np.array(vals)
 return samples,out

cnasamples,cna=matrix(download('data_cna.txt'));assert len(cnasamples)==281
patients=[re.sub(r'-PT-\d+$','',s) for s in cnasamples];assert len(set(patients))==281
rawcase=f'https://raw.githubusercontent.com/cBioPortal/datahub/{COMMIT}/public/coad_silu_2022/case_lists/cases_sequenced.txt'
with urllib.request.urlopen(rawcase,timeout=60) as r:case=r.read().decode()
seq=case.split('case_list_ids:',1)[1].strip().split('\t');assert set(seq)==set(cnasamples)
mut=defaultdict(set);classes={'Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Splice_Site','Translation_Start_Site'}
with download('data_mutations.txt').open(encoding='utf-8-sig',newline='') as f:
 reader=csv.DictReader((line for line in f if not line.startswith('#')),delimiter='\t')
 assert not any('driver' in h.lower() or 'oncogenic' in h.lower() for h in reader.fieldnames),'Annotation schema changed; review before classifying.'
 for r in reader:
  if r['Hugo_Symbol'] in GS and r['Variant_Classification'] in classes and r['Tumor_Sample_Barcode'] in set(seq):mut[r['Hugo_Symbol']].add(r['Tumor_Sample_Barcode'])
rnasamples,rna=matrix(download('data_mrna_seq_expression_zscores_ref_all_samples.txt'));assert len(rnasamples)==348
cyto={r['Gene']:r['Cytoband'] for r in csv.DictReader((D/'acicam_rna_cytoband.csv').open())}

for filename,headers in [('acicam_cna_counts.csv',['Gene','N','AMP','DEL']),('acicam_mutation_counts.csv',['Gene','N_mut','AnyCodingSplice_mut','N_joint','DNA_anycoding_or_CNA']),('acicam_rna_cytoband.csv',['Gene','N_RNA','RNA_high','RNA_low','Cytoband'])]:
 rows=[]
 for g in G:
  if filename.startswith('acicam_cna'):
   row=[g,int(np.isfinite(cna[g]).sum()),int((cna[g]==2).sum()),int((cna[g]==-2).sum())] if g in cna else [g,'NA','NA','NA']
  elif filename.startswith('acicam_mutation'):
   if g in cna:
    covered={s for j,s in enumerate(cnasamples) if np.isfinite(cna[g][j])};cnaset={s for j,s in enumerate(cnasamples) if abs(cna[g][j])==2};joint=covered&set(seq)
    row=[g,len(seq),len(mut[g]),len(joint),len((mut[g]|cnaset)&joint)]
   else:row=[g,'NA','NA','NA','NA']
  else:
   row=[g,int(np.isfinite(rna[g]).sum()),int((rna[g]>2).sum()),int((rna[g]<-2).sum()),cyto.get(g,'NA')] if g in rna else [g,'NA','NA','NA',cyto.get(g,'NA')]
  rows.append(row)
 # Verify counts before replacing existing archived values.
 out=D/filename
 if out.exists():
  old=list(csv.reader(out.open()));assert old[0]==headers
  assert [[str(v) for v in row] for row in rows]==old[1:],f'Refetched values differ: {filename}'
 with out.open('w',newline='') as f:w=csv.writer(f);w.writerow(headers);w.writerows(rows)
 print('Verified/refetched',filename)
print('AC-ICAM source summaries reproduced without driver inference.')
