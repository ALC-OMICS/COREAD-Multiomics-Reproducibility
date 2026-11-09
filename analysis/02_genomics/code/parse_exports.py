"""Read supplied categorical OncoPrint exports without guessing missing profiles.
Original workbooks remain immutable. OOXML is inspected with Python stdlib;
all final workbook creation uses artifact_tool.
"""
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET
from collections import Counter
import json, hashlib, csv
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'input'
NS={'x':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
FILES={'COAD_GDC':'PATIENT_DATA_oncoprint-COLON-TCGA-GDC.xlsx','READ_GDC':'PATIENT_DATA_oncoprint-RECTAL-TCGA-GDC.xlsx','CPTAC2':'PATIENT_DATA_oncoprint-CPTAC2.xlsx','PANCAN_P1':'PATIENT_DATA_oncoprint-PANCANCER-PARTE 1-.xlsx','PANCAN_P2':'PATIENT_DATA_oncoprint-PANCANCER-PARTE -2.xlsx'}
GENES=(ROOT/'input/union_candidates_466.txt').read_text().split()
assert len(GENES)==len(set(GENES))==466
GI={g:i for i,g in enumerate(GENES)}
MODS=['Driver_mutation','Driver_explicit','Splice_driver_code','Any_exported_mutation','Amplification','Deep_deletion','RNA_high','RNA_low','Protein_high','Protein_low','Structural_variant']

def colindex(s):
 x=0
 for c in s:
  if not c.isalpha():break
  x=x*26+ord(c.upper())-64
 return x-1

def rows(path):
 with ZipFile(path) as z:
  ss=[]
  if 'xl/sharedStrings.xml' in z.namelist():
   ss=[''.join(t.itertext()) for t in ET.fromstring(z.read('xl/sharedStrings.xml'))]
  root=ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
  for row in root.findall('x:sheetData/x:row',NS):
   d={}
   for c in row:
    v=c.find('x:v',NS)
    val=v.text if v is not None else ''.join(c.itertext())
    if c.get('t')=='s':val=ss[int(val)]
    d[colindex(c.get('r'))]=val
   yield int(row.get('r')),d

def parse(key,path):
 rr=iter(rows(path));_,h=next(rr)
 pats=[h[i] for i in range(2,max(h)+1)]
 assert len(pats)==len(set(pats))
 n=len(pats);events={m:np.zeros((466,n),dtype=np.int8) for m in MODS}
 tracks={t:np.zeros(466,dtype=bool) for t in ['CNA','MUTATIONS','MRNA','PROTEIN','STRUCTURAL_VARIANT']}
 all_gene_set=set();labels={};clinical={};out=[]
 for rn,d in rr:
  gene=d.get(0,'');typ=d.get(1,'')
  if typ=='CLINICAL':
   clinical[gene]=[d.get(i+2,'') for i in range(n)];continue
  if typ=='MUTATIONS':all_gene_set.add(gene)
  if gene not in GI:continue
  i=GI[gene]
  if typ in tracks:tracks[typ][i]=True
  for j0,v in d.items():
   if j0<2 or not v:continue
   j=j0-2;labels[v]=labels.get(v,0)+1
   hit=[]
   low=v.lower()
   if typ=='MUTATIONS':
    hit.append('Any_exported_mutation')
    if '(putative driver)' in low:hit += ['Driver_mutation','Driver_explicit']
    if 'splice_rec' in low:hit += ['Driver_mutation','Splice_driver_code']
   elif typ=='CNA':
    if 'amplification' in low or 'amp_rec' in low:hit.append('Amplification')
    if 'deep deletion' in low or 'homdel_rec' in low:hit.append('Deep_deletion')
   elif typ=='MRNA':
    if 'mrna high' in low:hit.append('RNA_high')
    if 'mrna low' in low:hit.append('RNA_low')
   elif typ=='PROTEIN':
    if 'protein high' in low:hit.append('Protein_high')
    if 'protein low' in low:hit.append('Protein_low')
   elif typ=='STRUCTURAL_VARIANT':hit.append('Structural_variant')
   if not hit:raise ValueError((key,gene,typ,v))
   for m in set(hit):events[m][i,j]=1
   out.append([key,gene,pats[j],typ,v])
 availability={}
 names={'MUTATIONS':'Profiled in Mutations','CNA':'Profiled in Putative copy-number alterations from GISTIC'}
 for typ,field in names.items():
  if field in clinical:availability[typ]=np.array([v=='Yes' for v in clinical[field]])
 for typ,fragment in [('MRNA','mRNA'),('PROTEIN','Protein')]:
  fs=[f for f in clinical if f.startswith('Profiled in') and fragment in f]
  if len(fs)==1: availability[typ]=np.array([v=='Yes' for v in clinical[fs[0]]])
 no_samples=np.array([int(v) if v else 1 for v in clinical.get('# Samples per Patient',['1']*n)])
 cohort={'key':key,'patients':pats,'genes':GENES,'events':events,'tracks':tracks,'availability':availability,'sample_counts':no_samples,'clinical':clinical}
 audit={'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'patients':n,'gene_rows':len(all_gene_set),'candidate_tracks':{k:int(v.sum()) for k,v in tracks.items()},'profile_yes':{k:int(v.sum()) for k,v in availability.items()},'multiple_samples':{pats[j]:int(v) for j,v in enumerate(no_samples) if v>1},'candidate_event_labels':labels}
 return cohort,audit,out

def save_npz(cohort):
 key=cohort['key'];kw={f'event__{m}':a for m,a in cohort['events'].items()};kw.update({f'avail__{m}':a for m,a in cohort['availability'].items()});kw.update({f'track__{m}':a for m,a in cohort['tracks'].items()});kw.update(patients=np.array(cohort['patients']),genes=np.array(GENES),sample_counts=cohort['sample_counts'])
 np.savez_compressed(ROOT/f'results/{key}_calls.npz',**kw)
 (ROOT/f'provenance/{key}_clinical.json').write_text(json.dumps(cohort['clinical'],indent=2))

if __name__=='__main__':
 cohorts={};audits=[];allout=[]
 for key,fn in FILES.items():
  c,a,o=parse(key,INPUT/fn);cohorts[key]=c;audits.append(a);allout+=o;save_npz(c);print(key,a['patients'],a['profile_yes'],a['multiple_samples'])
 c1=cohorts['COAD_GDC'];c2=cohorts['READ_GDC'];assert not(set(c1['patients'])&set(c2['patients']))
 merged={'key':'COREAD_GDC','patients':c1['patients']+c2['patients'],'events':{m:np.concatenate([c1['events'][m],c2['events'][m]],axis=1) for m in MODS},'availability':{m:np.r_[c1['availability'][m],c2['availability'][m]] for m in c1['availability']},'tracks':{m:c1['tracks'][m]&c2['tracks'][m] for m in c1['tracks']},'sample_counts':np.r_[c1['sample_counts'],c2['sample_counts']],'clinical':{}}
 save_npz(merged)
 # Merge PanCancer DNA only. Expression in P1 is absent, and P2 lacks profile flags.
 p1=cohorts['PANCAN_P1'];p2=cohorts['PANCAN_P2'];pidx={p:j for j,p in enumerate(p1['patients'])};j2=[pidx[p] for p in p2['patients']]
 p={'key':'PANCAN','patients':p1['patients'],'events':{m:p1['events'][m].copy() for m in MODS},'availability':p1['availability'],'tracks':{m:p1['tracks'][m]|p2['tracks'][m] for m in p1['tracks']},'sample_counts':p1['sample_counts'],'clinical':{}}
 # Each gene retains the patients present in its own export; missing 70 patients for P2-only genes are NOT zeros.
 coverage=np.zeros((466,len(p1['patients'])),dtype=bool)
 coverage[p1['tracks']['MUTATIONS'],:]=True
 for i,g in enumerate(GENES):
  if p2['tracks']['MUTATIONS'][i]:
   if p1['tracks']['MUTATIONS'][i]:
    for m in ['Driver_mutation','Any_exported_mutation','Amplification','Deep_deletion']:
     assert np.array_equal(p1['events'][m][i,j2],p2['events'][m][i]),(g,m,'Conflicting PanCancer overlapping values')
   for m in MODS:p['events'][m][i,j2]=p2['events'][m][i]
   coverage[i,j2]=True
 p['availability']={k:coverage&v[None,:] for k,v in p1['availability'].items() if k in ['CNA','MUTATIONS']}
 # Unverified protein gene coverage and RNA profile availability: unavailable as inferential endpoints.
 save_npz(p)
 with (ROOT/'results/exported_events_long.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['Cohort','Gene','Patient','Track','Exported_call']);w.writerows(allout)
 (ROOT/'provenance/input_audit.json').write_text(json.dumps(audits,indent=2))
 print('CANDIDATES',len(GENES),'COREAD DNA joint',int((merged['availability']['CNA']&merged['availability']['MUTATIONS']).sum()))
 for k,c in {**cohorts,'COREAD_GDC':merged}.items():
  for m,t in [('Driver_mutation','MUTATIONS'),('Amplification','CNA'),('Deep_deletion','CNA'),('RNA_high','MRNA'),('RNA_low','MRNA')]:
   if t in c['availability']:
    outside=int(c['events'][m][:,~c['availability'][t]].sum())
    if outside:print('WARNING_POSITIVE_OUTSIDE_PROFILE',k,m,outside)
