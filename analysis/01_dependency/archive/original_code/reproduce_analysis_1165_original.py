from pathlib import Path
import csv, re, json, zipfile, xml.etree.ElementTree as ET
import numpy as np
from collections import Counter

ROOT = Path(__file__).resolve().parent / "inputs"
OUT = Path(__file__).resolve().parent / "outputs"
OUT.mkdir(exist_ok=True)

def read_text_matrix(path):
    text=path.read_text(encoding="utf-8-sig")
    rows=list(csv.reader(text.splitlines(),delimiter="\t"))
    idrow=next(i for i,r in enumerate(rows) if "DepMap_ID" in r)
    header=next(i for i,r in enumerate(rows) if r and r[0]=="Dep_Map_Gene")
    positions=[j for j,x in enumerate(rows[idrow]) if re.fullmatch(r"ACH-\d{6}", x.strip())]
    records=[r for r in rows[header+1:] if r and r[0].strip()]
    shape_counts=Counter(map(len,records))
    bad=[]
    arr=np.full((len(records),len(positions)),np.nan)
    for i,r in enumerate(records):
        for k,j in enumerate(positions):
            v=r[j].strip() if j<len(r) else ""
            if v and v.lower() not in {"na","nan","n/a","#n/a","null"}:
                try: arr[i,k]=float(v.replace(",","."))
                except ValueError: bad.append((i+header+2,j+1,v))
    return {
        "rows":rows,"idrow":idrow,"header":header,"pos":positions,
        "ids":[rows[idrow][j].strip() for j in positions],
        "names":[rows[idrow+1][j].strip() for j in positions],
        "diagnoses":[rows[idrow+2][j].strip() for j in positions],
        "lineages":[rows[header][j].strip() for j in positions],
        "records":records,"symbols":[r[0].strip() for r in records],
        "values":arr,"bad":bad,"shape_counts":dict(shape_counts)}
new=read_text_matrix(ROOT/"Texto pegado(1).txt")
old=read_text_matrix(ROOT/"Texto pegado.txt")
print("New rows/cols",new["values"].shape,"column positions",new["pos"][:5],new["pos"][-5:],
      "widths",new["shape_counts"],"invalid",new["bad"][:5],
      "range",np.nanmin(new["values"]),np.nanmax(new["values"]))
print("Unique IDs",len(set(new["ids"])), "lineages",Counter(new["lineages"]))
print("Non-Cancerous",[(i,n,d) for i,n,d in zip(new["ids"],new["names"],new["diagnoses"]) if d=="Non-Cancerous"])
print("COREAD",Counter(new["diagnoses"])["Colorectal Adenocarcinoma"])
print("Missing",np.isnan(new["values"]).sum(),"duplicate genes",len(new["symbols"])-len(set(new["symbols"])))
print("Gene nonmatching raw",set(new["symbols"])-set(old["symbols"]))
print("all source prefixes", [r[:5] for r in new["records"][:6]])
print("md lists",[(p.name,p.read_text()[:55]) for p in ROOT.glob("*.md")])

# Parse stored Excel values with the standard library, to avoid reevaluating source formulas.
NS={"s":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
def read_xlsx_cells(path, sheet_number=1):
    with zipfile.ZipFile(path) as z:
        strings=[]
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("s:si",NS):
                strings.append("".join(t.text or "" for t in si.findall(".//s:t",NS)))
        root=ET.fromstring(z.read(f"xl/worksheets/sheet{sheet_number}.xml"))
        cells={}
        for c in root.findall(".//s:sheetData/s:row/s:c",NS):
            v=c.find("s:v",NS)
            if c.attrib.get("t")=="inlineStr":
                value="".join(t.text or "" for t in c.findall(".//s:t",NS))
            elif v is not None:
                value=v.text
                if c.attrib.get("t")=="s": value=strings[int(value)]
                elif c.attrib.get("t") not in {"str","e"}:
                    try: value=float(value); value=int(value) if value.is_integer() else value
                    except (ValueError,TypeError): pass
            else: value=None
            cells[c.attrib["r"]]=value
    return cells
def colname(n):
    out=""
    while n: n,a=divmod(n-1,26);out=chr(65+a)+out
    return out
def read_gene_table(path):
    cells=read_xlsx_cells(path)
    header_row=next(int(re.search(r"\d+",k)[0]) for k,v in cells.items() if v=="Input_order")
    headers=[]
    for j in range(1,300):
        h=cells.get(f"{colname(j)}{header_row}")
        if h is None: break
        headers.append(h)
    data=[]
    for i in range(header_row+1,header_row+1+1720):
        data.append({h:cells.get(f"{colname(j+1)}{i}") for j,h in enumerate(headers)})
    return headers,data,cells,header_row
org_h,orgrows,orgcells,org_header = read_gene_table(ROOT/"Supplementary_Table_COREAD_1720_Dependency_Candidates.xlsx")
prev_h,prevrows,prevcells,prev_header=read_gene_table(ROOT/"Supplementary_Table_COREAD_Cell_Lines.xlsx")
ih,irows,icells,iheader=read_gene_table(ROOT/"Supplementary_Table_COREAD_Cell_Lines_Organoids_Integration.xlsx")
genes=[r["Gene_symbol"] for r in orgrows]
assert len(genes)==len(set(genes))==1720
print("Org header",org_header,"headers",len(org_h),"genes",len(genes),genes[:3],genes[-3:])
print("Previous notes",[(k,v) for k,v in prevcells.items() if k.startswith("A") and isinstance(v,str) and ('bootstrap' in v.lower() or 'pair' in v.lower() or 'Cellosaurus' in v)][:20])
print("New unmatched",set(new["symbols"])-set(genes))
print("New matched",len(set(new["symbols"]) & set(genes)))
print("New trailing lines",new["symbols"][-5:])
json.dump({"new_shape":new["values"].shape,"lineages":dict(Counter(new["lineages"]))},open(OUT/"initial_audit.json","w"))

# --- Next execution step ---

# Extract provenance and exact procedures recorded in the earlier workbooks.
def physical_rows(cells, start, end, cols=16):
    return [[cells.get(f"{colname(j)}{i}") for j in range(1,cols+1)] for i in range(start,end+1)]
for i, r in enumerate(physical_rows(prevcells,1730,1776,16),1730):
    vals=[str(v) for v in r if v is not None]
    if vals: print(i," | ".join(vals)[:1600])
print("\nIntegration notes")
for i, r in enumerate(physical_rows(icells,1730,1780,16),1730):
    vals=[str(v) for v in r if v is not None]
    if vals: print(i," | ".join(vals)[:1500])
# Compare actual supplied values by normalized row symbol and ACH ID.
norm=lambda s: re.sub(r"\s+","",s)
oldsym={norm(s):i for i,s in enumerate(old["symbols"])}
newid={s:j for j,s in enumerate(new["ids"])}
old_indices=[newid[s] for s in old["ids"]]
matched_old=[s for s in new["symbols"] if s in oldsym]
aa=np.array([new["values"][new["symbols"].index(s),old_indices] for s in matched_old])
bb=np.array([old["values"][oldsym[s]] for s in matched_old])
same=(aa==bb)|(np.isnan(aa)&np.isnan(bb))
print("Shared genes",len(matched_old),"shared measurements",same.size,"different",np.count_nonzero(~same))
print("Metadata diagnoses changed",[(i,old["diagnoses"][j],new["diagnoses"][newid[i]]) for j,i in enumerate(old["ids"]) if old["diagnoses"][j]!=new["diagnoses"][newid[i]]])
print("Original source ID name rows",[(n,i) for n,i in zip(new["names"],new["ids"]) if n in {"DLD-1","HCT-15","COLO 201","COLO 205","SW 620","SW480"}])
# Organotypic source methods/models can be reused exactly as source annotations.
print("Org model annotations")
for k,v in orgcells.items():
    if isinstance(v,str) and (v in {"ModelID","Model_ID","SC-002839.CD01"} or "Source"==v):
        print(k, v)
print("Stored citations")
for k,v in prevcells.items():
    if isinstance(v,str) and ("cellosaurus" in v.lower() or "doi.org" in v or "scipy.org" in v):
        print(k,v[:450])

# --- Next execution step ---

import scipy, scipy.stats as ss
from scipy.stats.contingency import odds_ratio
from functools import lru_cache
import platform, hashlib, time
from threadpoolctl import threadpool_limits
threadpool_limits(limits=2)
# Keep source text symbols unchanged; this new file has exact query matches.
gix={g:i for i,g in enumerate(genes)}
sourceix={g:i for i,g in enumerate(new["symbols"])}
order=[g for g in genes if g in sourceix]
X=new["values"][[sourceix[g] for g in order]].copy()
D=np.array(new["diagnoses"])
L=np.array(new["lineages"])
I=np.array(new["ids"])
C=D=="Colorectal Adenocarcinoma"
NC=D=="Non-Cancerous"
O=~C & ~NC
OA=~C
assert C.sum()==63 and O.sum()==1099 and NC.sum()==3 and len(I)==1165
print("scipy",scipy.__version__,"numpy",np.__version__,"python",platform.python_version())
print("Eligible other lineages",[(s,int(np.sum(O&(L==s)))) for s in sorted(set(L[O]))])
print("clinical matched",len(order),"core gene counts",Counter(np.isfinite(X[:,C]).sum(axis=1)))
print("CYP2E1 counts",np.sum(X[order.index("CYP2E1"),C]>.5))

def bh(p):
    p=np.asarray(p,float);q=np.full_like(p,np.nan); ok=np.isfinite(p)
    m=ok.sum()
    if m:
        v=p[ok];idx=np.argsort(v,kind="stable")
        adj=np.minimum.accumulate((v[idx]*m/np.arange(1,m+1))[::-1])[::-1]
        tmp=np.empty(m);tmp[idx]=np.minimum(adj,1)
        q[ok]=tmp
    return q

def wilson(k,n):
    k=np.asarray(k,float);n=np.asarray(n,float)
    p=np.divide(k,n,out=np.full(np.broadcast(k,n).shape,np.nan),where=n>0)
    z=ss.norm.ppf(.975)
    denom=1+z*z/n
    center=(p+z*z/(2*n))/denom
    margin=z*np.sqrt((p*(1-p)+z*z/(4*n))/n)/denom
    return np.maximum(0,center-margin),np.minimum(1,center+margin)

def summarize(x,mask,cut=.5):
    vals=x[:,mask]; valid=np.isfinite(vals);n=valid.sum(1);k=((vals>cut)&valid).sum(1)
    return {"n":n,"k":k,"freq":k/n, "mean":np.nanmean(vals,1),
            "median":np.nanmedian(vals,1), "q25":np.nanquantile(vals,.25,axis=1),
            "q75":np.nanquantile(vals,.75,axis=1), "sd":np.nanstd(vals,axis=1,ddof=1),
            "min":np.nanmin(vals,1),"max":np.nanmax(vals,1),"missing":mask.sum()-n,
            "low":wilson(k,n)[0],"high":wilson(k,n)[1]}

@lru_cache(maxsize=None)
def fp(a,b,c,d):
    return float(ss.fisher_exact([[int(a),int(b)],[int(c),int(d)]],alternative="two-sided").pvalue)

def compare(n1,k1,n0,k0):
    eligible=(n1>=5)&(n0>=5)
    p=np.full(len(n1),np.nan)
    for j in np.flatnonzero(eligible):
        p[j]=fp(int(k1[j]),int(n1[j]-k1[j]),int(k0[j]),int(n0[j]-k0[j]))
    f1=np.divide(k1,n1);f0=np.divide(k0,n0); delta=f1-f0
    lo1,hi1=wilson(k1,n1);lo0,hi0=wilson(k0,n0)
    return {"p":p,"q":bh(p),"delta":delta,
            "low":delta-np.sqrt((f1-lo1)**2+(hi0-f0)**2),
            "high":delta+np.sqrt((hi1-f1)**2+(f0-lo0)**2),"eligible":eligible}
t=time.time()
S=summarize(X,C); T=summarize(X,O); Tall=summarize(X,OA); S9=summarize(X,C,.9);T9=summarize(X,O,.9)
test=compare(S["n"],S["k"],T["n"],T["k"])
test9=compare(S9["n"],S9["k"],T9["n"],T9["k"])
testall=compare(S["n"],S["k"],Tall["n"],Tall["k"])
other_lineages=sorted(set(L[O])); eligible_lineages=[s for s in other_lineages if np.sum(O&(L==s))>=5]
linstats={};lintests={}
for s in other_lineages:
    lst=summarize(X,O&(L==s));linstats[s]=lst
    if s in eligible_lineages:lintests[s]=compare(S["n"],S["k"],lst["n"],lst["k"])
joint=np.stack([lintests[s]["p"] for s in eligible_lineages],axis=1)
jointq=bh(joint.ravel()).reshape(joint.shape)
for j,s in enumerate(eligible_lineages):lintests[s]["jointq"]=jointq[:,j]
print("stats elapsed",time.time()-t,"Fisher family",sum(test["eligible"]),"joint",np.isfinite(joint).sum(),
      "eligible lineages",len(eligible_lineages))
print("Candidate",int((S["k"]>=2).sum()),"highp",int((S9["k"]>=2).sum()),"primary higher",sum((test["q"]<.05)&(test["delta"]>0)),"lower",sum((test["q"]<.05)&(test["delta"]<0)))
prefer=[(g,int(S["k"][j]),int(T["k"][j]),float(test["q"][j]),float(test["delta"][j])) for j,g in enumerate(order) if test["q"][j]<.05 and test["delta"][j]>0]
print("Higher", sorted(prefer,key=lambda x:x[3]))
print("All other1102 higher",sum((testall["q"]<.05)&(testall["delta"]>0)))
print("candidate+counts",[(cut,int((S["freq"]>=cut).sum())) for cut in [.125,.25,.5,.8,1]])
np.savez_compressed(OUT/"input_numeric.npz",X=X,genes=np.array(order),ids=I,diagnoses=D,lineages=L,coread_mask=C)

# --- Next execution step ---

t=time.time()
@lru_cache(maxsize=None)
def conditional_or(a,b,c,d):
    if a+c==0 or b+d==0: return (np.nan,0.0,np.inf,"NOT_IDENTIFIABLE_NO_VARIATION")
    result=odds_ratio([[a,b],[c,d]],kind="conditional")
    ci=result.confidence_interval(confidence_level=.95)
    s="INFINITE" if np.isinf(result.statistic) else ("ZERO" if result.statistic==0 else "FINITE")
    return (float(result.statistic),float(ci.low),float(ci.high),s)
or_rows=[conditional_or(int(a),int(n-a),int(c),int(m-c)) for a,n,c,m in zip(S["k"],S["n"],T["k"],T["n"])]
print("OR seconds",time.time()-t)
# All representative combinations for the two previously documented related pairs.
pairids=[["ACH-000997","ACH-001061"],["ACH-000253","ACH-001039"]]
import itertools
pair_sens=[]
for retain in itertools.product(*pairids):
    mask=C.copy()
    for pair,rep in zip(pairids,retain):
        mask &= ~(np.isin(I,[s for s in pair if s!=rep]))
    s=summarize(X,mask);tst=compare(s["n"],s["k"],T["n"],T["k"])
    pair_sens.append({"retained":retain,"mask":mask,**s,**{"p":tst["p"],"q":tst["q"],"delta":tst["delta"]}})
pair_n=np.stack([s["n"] for s in pair_sens]);pair_k=np.stack([s["k"] for s in pair_sens])
pair_freq=np.stack([s["freq"] for s in pair_sens]);pair_q=np.stack([s["q"] for s in pair_sens])
robust=(pair_k>=2).all(0)
print("Related-pair robust",sum(robust),"FDR worst higher",sum((np.max(pair_q,axis=0)<.05)&(test["delta"]>0)))
balance_freq=np.stack([linstats[s]["freq"] for s in eligible_lineages],1)
balance_elig=(S["n"]>=5)&np.all(np.stack([linstats[s]["n"] for s in eligible_lineages],1)>=5,axis=1)
balance_other=np.where(balance_elig,balance_freq.mean(1),np.nan)
balance_delta=S["freq"]-balance_other
loo=S["freq"][:,None]-(balance_freq.sum(1)[:,None]-balance_freq)/(len(eligible_lineages)-1)
loo_min=np.where(balance_elig,np.min(loo,1),np.nan)
loo_max=np.where(balance_elig,np.max(loo,1),np.nan)
broad=np.where(balance_elig,((S["freq"]>=.8)&np.all(balance_freq>=.8,axis=1)).astype(float),np.nan)
print("Balance eligible",sum(balance_elig),"broad",np.nansum(broad))
# Bootstrap entire model columns within each fixed lineage: same algorithm as prior study.
# Pre-create weights, use float BLAS, retain gene-by-model missingness.
boot_rng=np.random.default_rng(20261003)
eligible_i=np.flatnonzero(balance_elig)
B=10000;G=len(eligible_i)
bootstrap_diff=np.zeros((B,G),dtype=np.float64)
boot_masks=[C]+[O&(L==s) for s in eligible_lineages]
print("Bootstrap ready:",B,"x",G, "in",len(boot_masks),"strata")
np.savez_compressed(OUT/"stats_stage1.npz", p=test["p"],q=test["q"],
                   pair_q=pair_q, balance_elig=balance_elig)

# --- Next execution step ---

t=time.time()
for j,mask in enumerate(boot_masks[:13]):
    n_total=int(mask.sum())
    weights=boot_rng.multinomial(n_total,np.full(n_total,1/n_total),size=B).astype(float)
    v=X[eligible_i][:,mask]
    den=weights @ np.isfinite(v).astype(float).T
    num=weights @ (v>.5).astype(float).T
    freq=np.divide(num,den,out=np.full_like(num,np.nan),where=den>0)
    bootstrap_diff += freq if j==0 else -freq/len(eligible_lineages)
print("Bootstrap first13 strata seconds:",round(time.time()-t,2))
del weights,den,num,freq,v

# --- Next execution step ---

t=time.time()
for j in range(13,len(boot_masks)):
    mask=boot_masks[j]; n_total=int(mask.sum())
    weights=boot_rng.multinomial(n_total,np.full(n_total,1/n_total),size=B).astype(float)
    v=X[eligible_i][:,mask]
    den=weights @ np.isfinite(v).astype(float).T
    num=weights @ (v>.5).astype(float).T
    freq=np.divide(num,den,out=np.full_like(num,np.nan),where=den>0)
    bootstrap_diff -= freq/len(eligible_lineages)
boot_ci=np.full((2,len(order)),np.nan)
boot_ci[:,eligible_i]=np.nanquantile(bootstrap_diff,[.025,.975],axis=0)
boot_usable=np.zeros(len(order),int);boot_usable[eligible_i]=np.isfinite(bootstrap_diff).sum(0)
print("Bootstrap end seconds:",round(time.time()-t,2),"usable",Counter(boot_usable[eligible_i]),
      "higher balanced",sum((test["q"]<.05)&(test["delta"]>0)&(boot_ci[0]>0)))
np.savez_compressed(OUT/"bootstrap_summary.npz",low=boot_ci[0],high=boot_ci[1],usable=boot_usable)
del weights,den,num,freq,v,bootstrap_diff
# Check original full universe.
original_gene_lines=(ROOT/"Markdown pegado(20261001-223143).md").read_text().splitlines()
original_genes=[l.strip().strip("|").strip() for l in original_gene_lines if "|" in l and not re.fullmatch(r"[|:\-\s]+",l)]
assert original_genes==genes
# Build per-input aligned data and matrix pairs; updated organoid calculations.
idx_order={g:i for i,g in enumerate(order)}
orgmap={r["Gene_symbol"]:r for r in orgrows}
gh=[h for h in org_h if h.startswith("Chronos_3DO|")]
eh=[h for h in org_h if h.startswith("EB_REEST_3DO|")]
gh2=[h for h in org_h if h.startswith("Chronos_2DO|")]
eh2=[h for h in org_h if h.startswith("EB_REEST_2DO|")]
orgX=np.array([[r[h] if r[h] is not None else np.nan for h in gh] for r in orgrows],float)
orgP=np.array([[r[h] if r[h] is not None else np.nan for h in eh] for r in orgrows],float)
orgN=np.isfinite(orgX).sum(1);orgK=(orgX<-.5).sum(1)
orgcand=orgK>=2;orgAvail=orgN>0;orgCon=((orgX<-.5)&(orgP>.5)).sum(1)
lineN=np.zeros(len(genes),int);lineK=np.zeros(len(genes),int);linefreq=np.full(len(genes),np.nan)
for j,g in enumerate(genes):
    if g in idx_order:
        k=idx_order[g];lineN[j]=S["n"][k];lineK[j]=S["k"][k];linefreq[j]=S["freq"][k]
linecand=lineK>=2;lineAvail=lineN>0
bothAvail=lineAvail&orgAvail
bothcand=linecand&orgcand
union=linecand|orgcand
print("Integration","line",sum(linecand),"org",sum(orgcand),"both avail",sum(bothAvail),
      "both cand",sum(bothcand),"union",sum(union))
patterns=[]
for i,g in enumerate(genes):
    if bothAvail[i]:
        patterns.append("BOTH_CANDIDATES" if bothcand[i] else "CL_ONLY_BOTH_EVALUABLE" if linecand[i] else "ORG_ONLY_BOTH_EVALUABLE" if orgcand[i] else "NEITHER_BOTH_EVALUABLE")
    elif lineAvail[i]:
        patterns.append("CL_CANDIDATE_ORG_UNAVAILABLE" if linecand[i] else "CL_NONCANDIDATE_ORG_UNAVAILABLE")
    elif orgAvail[i]:
        patterns.append("ORG_CANDIDATE_CL_UNAVAILABLE" if orgcand[i] else "ORG_NONCANDIDATE_CL_UNAVAILABLE")
    else:patterns.append("NEITHER_EVALUABLE")
print(Counter(patterns))
print("joint proportion>=80",np.sum((linefreq>=.8)&(orgK/np.maximum(orgN,1)>=.8)))
print("New gene",orgmap["CYP2E1"]["Data_status"],orgmap["CYP2E1"]["Primary_candidate_3DO"])
# Source organoid ScreenID to ModelID metadata preserved in old file.
orgmodels=physical_rows(orgcells,1773,1791,11)
print("Org models first",orgmodels[:3])
print("source mdhash",hashlib.sha256((ROOT/"Texto pegado(1).txt").read_bytes()).hexdigest())

# --- Next execution step ---

t=time.time()
joint_orig_i=np.flatnonzero(bothAvail)
joint_source_i=[idx_order[genes[i]] for i in joint_orig_i]
xcl=X[joint_source_i][:,C]
xorg=orgX[joint_orig_i]
fc=(xcl>.5).sum(1)/np.isfinite(xcl).sum(1)
fo=(xorg<-.5).sum(1)/np.isfinite(xorg).sum(1)
fixed_nonbroad=~((fc>=.8)&(fo>=.8))
def rankcorr_rows(a,b):
    ra=ss.rankdata(a,axis=1);rb=ss.rankdata(b,axis=1)
    ra-=ra.mean(1,keepdims=True);rb-=rb.mean(1,keepdims=True)
    return (ra*rb).sum(1)/np.sqrt((ra*ra).sum(1)*(rb*rb).sum(1))
rho=ss.spearmanr(fc,fo).statistic
rho_non=ss.spearmanr(fc[fixed_nonbroad],fo[fixed_nonbroad]).statistic
rng=np.random.default_rng(20261004); rb=2000
wc=rng.multinomial(63,np.full(63,1/63),size=rb).astype(float)
wo=rng.multinomial(16,np.full(16,1/16),size=rb).astype(float)
cboot=(wc@(xcl>.5).T)/(wc@np.isfinite(xcl).T)
oboot=(wo@(xorg<-.5).T)/(wo@np.isfinite(xorg).T)
rhos=rankcorr_rows(cboot,oboot)
rhos_non=rankcorr_rows(cboot[:,fixed_nonbroad],oboot[:,fixed_nonbroad])
joint_counts={"both":int(bothcand.sum()),"line_only":sum(1 for p in patterns if p=="CL_ONLY_BOTH_EVALUABLE"),
              "org_only":sum(1 for p in patterns if p=="ORG_ONLY_BOTH_EVALUABLE"),
              "neither":sum(1 for p in patterns if p=="NEITHER_BOTH_EVALUABLE")}
integ_summary={
    "both_evaluable":int(bothAvail.sum()),"union":int(union.sum()),"candidate_line":int(linecand.sum()),
    "candidate_org":int(orgcand.sum()),"candidate_both":int(bothcand.sum()),
    "patterns":dict(Counter(patterns)),
    "jaccard_joint":joint_counts["both"]/(joint_counts["both"]+joint_counts["line_only"]+joint_counts["org_only"]),
    "dice_joint":2*joint_counts["both"]/(2*joint_counts["both"]+joint_counts["line_only"]+joint_counts["org_only"]),
    "rho":float(rho),"rho_ci":np.quantile(rhos,[.025,.975]).tolist(),
    "nonbroad_n":int(fixed_nonbroad.sum()),"rho_nonbroad":float(rho_non),
    "rho_nonbroad_ci":np.quantile(rhos_non,[.025,.975]).tolist(),
    "aligned_both":int(np.sum(bothcand&(linefreq>=.125))),
    "both_knownpairs":sum(1 for i,g in enumerate(genes) if bothcand[i] and robust[idx_order[g]]),
    "both_ge25":int(np.sum((linefreq>=.25)&(orgK/np.maximum(orgN,1)>=.25))),
    "both_ge50":int(np.sum((linefreq>=.5)&(orgK/np.maximum(orgN,1)>=.5))),
    "both_ge80":int(np.sum((linefreq>=.8)&(orgK/np.maximum(orgN,1)>=.8))),
}
print(json.dumps(integ_summary,indent=2),"seconds",time.time()-t)
raw_tot=np.array([float(r[2].replace(",",".")) for r in new["records"]])
new_mean=np.nanmean(new["values"],axis=1)
print("total averages mismatch",np.sum(np.abs(new_mean-raw_tot)>.000500001),"max",max(abs(new_mean-raw_tot)))
print("Cell models overlapping org",[r[1] for r in orgmodels[1:] if r[1] in I[C]])
print("allmatrix overlap organoids",[r[1] for r in orgmodels[1:] if r[1] in I])
print("balanced discordant 3",[(g,float(boot_ci[0][j]),float(balance_delta[j])) for j,g in enumerate(order) if test["q"][j]<.05 and test["delta"][j]>0 and not boot_ci[0][j]>0])
json.dump(integ_summary,open(OUT/"integration_summary.json","w"),indent=2)
# Save organoid and main source values for reproducibility.
with (OUT/"organoid_values.csv").open("w",newline="") as h:
    writer=csv.writer(h); writer.writerow(["Gene_symbol","Data_status",*gh,*eh,*gh2,*eh2])
    writer.writerows([[r["Gene_symbol"],r["Data_status"],*[r[h] for h in gh+eh+gh2+eh2]] for r in orgrows])
print("all1 genes test equality p>",np.sum(X[:,C]==.5),np.sum(X[:,C]==.9))
del wc,wo,cboot,oboot

# --- Next execution step ---

# Fresh results in row dictionaries; simple arithmetic will also be visible as Excel formulas.
def scalar(v):
    if isinstance(v,np.generic): v=v.item()
    if isinstance(v,float):
        if np.isnan(v): return None
        if np.isinf(v): return "Inf" if v>0 else "-Inf"
    return v
def rec(k,n):
    if not n:return "NOT_EVALUABLE"
    if k==0:return "NO_POSITIVE_CALL"
    if k==1:return "SINGLE_MODEL"
    f=k/n
    return "GE80PCT" if f>=.8 else "GE50_LT80PCT" if f>=.5 else "GE25_LT50PCT" if f>=.25 else "LT25PCT"
def yes(b):return "YES" if b else "NO"
def pref(q,d):
    if not np.isfinite(q):return "NOT_EVALUABLE"
    return "HIGHER_IN_COREAD" if q<.05 and d>0 else "LOWER_IN_COREAD" if q<.05 and d<0 else "NO_FDR_DIFFERENCE"
lin_slug={s:re.sub(r"[^A-Za-z0-9]+","_",str(s)).strip("_") for s in other_lineages}
meta_fields=["Input_order","Gene_symbol","Source_symbol","Match_status","Data_status","CL_candidate","CL_recurrence_band","Known_pairs_candidate_all4"]
core_fields=["N_COREAD_evaluable","N_COREAD_positive","COREAD_dep_frequency","COREAD_freq_CI95_low","COREAD_freq_CI95_high","COREAD_P_mean","COREAD_P_median","COREAD_P_Q25","COREAD_P_Q75","COREAD_P_SD","COREAD_P_min","COREAD_P_max","N_COREAD_missing","COREAD_P_gt0p9_n","COREAD_P_gt0p9_frequency","COREAD_P_gt0p9_candidate","COREAD_freq_ge12p5_candidate"]
compare_fields=["N_other_evaluable","N_other_positive","Other_dep_frequency","Other_freq_CI95_low","Other_freq_CI95_high","Other_P_mean","Other_P_median","Frequency_difference","Difference_CI95_low","Difference_CI95_high","Conditional_odds_ratio","OR_CI95_low","OR_CI95_high","OR_status","Fisher_two_sided_P","Fisher_BH_FDR","Frequency_preference","Tests_in_primary_BH_family"]
sens_fields=["HighP_Fisher_P","HighP_Fisher_BH_FDR","Balanced_other_frequency","Balanced_frequency_difference","Balanced_diff_boot_CI95_low","Balanced_diff_boot_CI95_high","Balanced_boot_usable_n","Balanced_eligible_lineages","Leave_one_lineage_delta_min","Leave_one_lineage_delta_max","Higher_BH_and_balanced_CI","Broad_in_measured_lineages","Known_pairs_n_min","Known_pairs_n_max","Known_pairs_positive_min","Known_pairs_positive_max","Known_pairs_frequency_min","Known_pairs_frequency_max","Known_pairs_FDR_max","Known_pairs_higher_all4"]
full_fields=["All_nonCORE_N","All_nonCORE_positive","All_nonCORE_frequency","All_nonCORE_Fisher_P","All_nonCORE_BH_FDR","All_nonCORE_preference","Excluded_noncancer_N","Excluded_noncancer_positive"]
raw_fields=[f"COREAD_P|{id_}" for id_ in I[C]]
prov_fields=["Source_text_row","Supplied_Total_average","Recomputed_mean_all1165","Source_filename","Source_release_status"]
clineages_fields=[]
for s in other_lineages:
    slug=lin_slug[s]
    clineages_fields+= [f"{slug}_{x}" for x in ("N","positive","frequency","P_mean")]
    if s in eligible_lineages: clineages_fields+=[f"{slug}_{x}" for x in ("delta","Fisher_P","joint_BH_FDR")]
clcols=meta_fields+core_fields+compare_fields+sens_fields+full_fields+clineages_fields+raw_fields+prov_fields
clrows=[]
for i,g in enumerate(genes):
    r={h:None for h in clcols}
    r.update(Input_order=i+1,Gene_symbol=g,Source_filename="Texto pegado(1).txt",
             Source_release_status="NOT_PROVIDED")
    if g not in idx_order:
        r.update(Match_status="NOT_IN_SUPPLIED_MATRIX",Data_status="NO_MATCH",CL_candidate="NOT_EVALUABLE",
                 CL_recurrence_band="NOT_EVALUABLE",Known_pairs_candidate_all4="NOT_EVALUABLE",N_COREAD_evaluable=0)
        clrows.append(r);continue
    j=idx_order[g];si=sourceix[g]
    r.update(Source_symbol=g,Match_status="EXACT",Data_status="EVALUABLE",CL_candidate=yes(S["k"][j]>=2),
             CL_recurrence_band=rec(S["k"][j],S["n"][j]),Known_pairs_candidate_all4=yes(robust[j]))
    for f,key in {"N_COREAD_evaluable":"n","N_COREAD_positive":"k","COREAD_dep_frequency":"freq",
        "COREAD_freq_CI95_low":"low","COREAD_freq_CI95_high":"high","COREAD_P_mean":"mean",
        "COREAD_P_median":"median","COREAD_P_Q25":"q25","COREAD_P_Q75":"q75","COREAD_P_SD":"sd",
        "COREAD_P_min":"min","COREAD_P_max":"max","N_COREAD_missing":"missing"}.items():r[f]=S[key][j]
    r.update(COREAD_P_gt0p9_n=S9["k"][j],COREAD_P_gt0p9_frequency=S9["freq"][j],
             COREAD_P_gt0p9_candidate=yes(S9["k"][j]>=2),
             COREAD_freq_ge12p5_candidate=yes(S["k"][j]>=2 and S["freq"][j]>=.125))
    for f,key in {"N_other_evaluable":"n","N_other_positive":"k","Other_dep_frequency":"freq",
         "Other_freq_CI95_low":"low","Other_freq_CI95_high":"high","Other_P_mean":"mean","Other_P_median":"median"}.items():r[f]=T[key][j]
    for f,key in {"Frequency_difference":"delta","Difference_CI95_low":"low","Difference_CI95_high":"high",
                  "Fisher_two_sided_P":"p","Fisher_BH_FDR":"q"}.items():r[f]=test[key][j]
    r.update(dict(zip(["Conditional_odds_ratio","OR_CI95_low","OR_CI95_high","OR_status"],or_rows[j])))
    r.update(Frequency_preference=pref(test["q"][j],test["delta"][j]),
             Tests_in_primary_BH_family=int(sum(test["eligible"])),
             HighP_Fisher_P=test9["p"][j],HighP_Fisher_BH_FDR=test9["q"][j],
             Balanced_other_frequency=balance_other[j],Balanced_frequency_difference=balance_delta[j],
             Balanced_diff_boot_CI95_low=boot_ci[0,j],Balanced_diff_boot_CI95_high=boot_ci[1,j],
             Balanced_boot_usable_n=boot_usable[j] if balance_elig[j] else None,
             Balanced_eligible_lineages=24 if balance_elig[j] else None,
             Leave_one_lineage_delta_min=loo_min[j],Leave_one_lineage_delta_max=loo_max[j],
             Higher_BH_and_balanced_CI=yes(test["q"][j]<.05 and test["delta"][j]>0 and boot_ci[0,j]>0) if balance_elig[j] else "NOT_EVALUABLE",
             Broad_in_measured_lineages=yes(broad[j]==1) if balance_elig[j] else "NOT_EVALUABLE",
             Known_pairs_n_min=pair_n[:,j].min(),Known_pairs_n_max=pair_n[:,j].max(),
             Known_pairs_positive_min=pair_k[:,j].min(),Known_pairs_positive_max=pair_k[:,j].max(),
             Known_pairs_frequency_min=pair_freq[:,j].min(),Known_pairs_frequency_max=pair_freq[:,j].max(),
             Known_pairs_FDR_max=pair_q[:,j].max(),
             Known_pairs_higher_all4=yes(np.all((pair_q[:,j]<.05) & np.array([s["delta"][j]>0 for s in pair_sens]))),
             All_nonCORE_N=Tall["n"][j],All_nonCORE_positive=Tall["k"][j],All_nonCORE_frequency=Tall["freq"][j],
             All_nonCORE_Fisher_P=testall["p"][j],All_nonCORE_BH_FDR=testall["q"][j],
             All_nonCORE_preference=pref(testall["q"][j],testall["delta"][j]),
             Excluded_noncancer_N=np.isfinite(X[j,NC]).sum(),Excluded_noncancer_positive=(X[j,NC]>.5).sum())
    for s in other_lineages:
        slug=lin_slug[s];st=linstats[s]
        for f,key in zip(["N","positive","frequency","P_mean"],["n","k","freq","mean"]):
            r[f"{slug}_{f}"]=st[key][j] if st["n"][j] else (0 if f=="N" else None)
        if s in eligible_lineages:
            for f,key in zip(["delta","Fisher_P","joint_BH_FDR"],["delta","p","jointq"]):r[f"{slug}_{f}"]=lintests[s][key][j]
    r.update({f:val for f,val in zip(raw_fields,X[j,C])})
    r.update(Source_text_row=si+5,Supplied_Total_average=float(new["records"][si][2].replace(",",".")),
             Recomputed_mean_all1165=np.nanmean(X[j]))
    clrows.append({h:scalar(v) for h,v in r.items()})
print("CL columns",len(clcols),"cells",len(clcols)*1720,"samples",clrows[0]["CL_candidate"],clrows[-1]["Gene_symbol"])
# Integration reuse existing field vocabulary and refresh each available annotation.
introws=[]
for i,g in enumerate(genes):
    oldr=irows[i];cr=clrows[i];o=orgmap[g]
    r=dict(oldr)  # values subsequently overwritten; no old CL statistical result may survive.
    r.update(Input_order=i+1,Gene_symbol=g,Evidence_pattern=patterns[i],
             Union_candidate=yes(union[i]) if lineAvail[i] or orgAvail[i] else "NOT_EVALUABLE",
             Both_systems_candidate=yes(bothcand[i]) if bothAvail[i] else "NOT_EVALUABLE",
             CL_candidate=cr["CL_candidate"],Organoid_3DO_candidate=yes(orgcand[i]) if orgAvail[i] else "NOT_EVALUABLE",
             Both_evaluable=yes(bothAvail[i]),CL_data_status=cr["Data_status"],Organoid_data_status=o["Data_status"])
    for dst,src in {"CL_N":"N_COREAD_evaluable","CL_positive":"N_COREAD_positive","CL_frequency":"COREAD_dep_frequency",
        "CL_freq_CI95_low":"COREAD_freq_CI95_low","CL_freq_CI95_high":"COREAD_freq_CI95_high",
        "CL_P_mean":"COREAD_P_mean","CL_P_median":"COREAD_P_median","CL_recurrence_band":"CL_recurrence_band",
        "CL_P_gt0p9_n":"COREAD_P_gt0p9_n","CL_P_gt0p9_frequency":"COREAD_P_gt0p9_frequency",
        "CL_P_gt0p9_candidate":"COREAD_P_gt0p9_candidate","CL_broad_measured_lineages":"Broad_in_measured_lineages",
        "CL_knownpairs_robust":"Known_pairs_candidate_all4","CL_other_frequency":"Other_dep_frequency",
        "CL_other_delta":"Frequency_difference","CL_Fisher_P":"Fisher_two_sided_P","CL_Fisher_BH_FDR":"Fisher_BH_FDR",
        "CL_preference":"Frequency_preference","CL_balanced_delta":"Balanced_frequency_difference",
        "CL_balanced_CI95_low":"Balanced_diff_boot_CI95_low","CL_balanced_CI95_high":"Balanced_diff_boot_CI95_high",
        "CL_source_symbol":"Source_symbol","CL_match_status":"Match_status","CL_source_row":"Source_text_row",
        "CL_freq_ge12p5":"COREAD_freq_ge12p5_candidate"}.items():r[dst]=cr[src]
    r["CL_freq_ge12p5"]= cr["COREAD_freq_ge12p5_candidate"] if lineAvail[i] else "NOT_EVALUABLE"
    if bothAvail[i]:
        j=idx_order[g]
        r["Both_CLhighP_Orgconcordance"]=yes(S9["k"][j]>=2 and orgCon[i]>=2)
        r["Both_ge12p5_recurrence"]=yes(linefreq[i]>=.125 and orgK[i]/orgN[i]>=.125)
        r["Both_freq_ge25pct"]=yes(linefreq[i]>=.25 and orgK[i]/orgN[i]>=.25)
        r["Both_freq_ge50pct"]=yes(linefreq[i]>=.5 and orgK[i]/orgN[i]>=.5)
        r["Both_freq_ge80pct"]=yes(linefreq[i]>=.8 and orgK[i]/orgN[i]>=.8)
        r["Both_with_CLknownpairs_robust"]=yes(robust[j] and orgcand[i])
    else:
        for f in ["Both_CLhighP_Orgconcordance","Both_ge12p5_recurrence","Both_freq_ge25pct","Both_freq_ge50pct","Both_freq_ge80pct","Both_with_CLknownpairs_robust"]:r[f]="NOT_EVALUABLE"
    # Recalculate organoid raw summaries to ensure not relying on old flags.
    r.update(Org_N_3DO=int(orgN[i]),Org_positive_3DO=int(orgK[i]) if orgN[i] else None,
             Org_frequency_3DO=float(orgK[i]/orgN[i]) if orgN[i] else None)
    r["CL_other_N"]=cr["N_other_evaluable"]
    r["CL_other_positive"]=cr["N_other_positive"]
    r["CL_All_nonCORE_BH_FDR"]=cr["All_nonCORE_BH_FDR"]
    r["CL_source_filename"]="Texto pegado(1).txt"
    r["CL_release_status"]="NOT_PROVIDED"
    introws.append({h:scalar(v) for h,v in r.items()})
intcols=list(introws[0])
assert all(r["CL_source_filename"]=="Texto pegado(1).txt" for r in introws)
print("Integration cols",len(intcols))
def write_dict_csv(path,cols,rows):
    with path.open("w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=cols);w.writeheader();w.writerows(rows)
write_dict_csv(OUT/"cell_line_results_1165.csv",clcols,clrows)
write_dict_csv(OUT/"integration_results_1165.csv",intcols,introws)
print("csv results written")

# --- Next execution step ---

# Assemble self-contained workbook notes and computational audit, with no invented source version.
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [
    ROOT/"Texto pegado(1).txt",ROOT/"Texto pegado.txt",
    ROOT/"Supplementary_Table_COREAD_1720_Dependency_Candidates.xlsx",
    ROOT/"Markdown pegado(20261001-223143).md"]}
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint, confint_proportions_2indep
import statsmodels
assert np.allclose(test["q"],multipletests(test["p"],method="fdr_bh")[1],atol=1e-14)
assert np.allclose(test9["q"],multipletests(test9["p"],method="fdr_bh")[1],atol=1e-14)
assert np.allclose(testall["q"],multipletests(testall["p"],method="fdr_bh")[1],atol=1e-14)
assert np.allclose(jointq[np.isfinite(joint)],multipletests(joint[np.isfinite(joint)],method="fdr_bh")[1],atol=1e-14)
for j in [0,1,15,30,80,300,900,1500]:
    lo,hi=proportion_confint(S["k"][j],S["n"][j],method="wilson")
    assert np.allclose([lo,hi],[S["low"][j],S["high"][j]],atol=1e-12)
    lo,hi=confint_proportions_2indep(S["k"][j],S["n"][j],T["k"][j],T["n"][j],method="newcomb",compare="diff")
    assert np.allclose([lo,hi],[test["low"][j],test["high"][j]],atol=1e-12)
# Every original COREAD value matches, after the explicitly corrected input symbol names.
assert same.all()
for i,g in enumerate(genes):
    if orgN[i]:
        assert int(orgK[i])==int(orgmap[g]["GE_lt_minus0p5_n_3DO"])
        assert bool(orgcand[i])==(orgmap[g]["Primary_candidate_3DO"]=="YES")
        assert np.allclose(np.nanmean(orgX[i]),orgmap[g]["Chronos_mean_3DO"])
oldset={r["Gene_symbol"] for r in prevrows if r["Frequency_preference"]=="HIGHER_IN_COREAD"}
newset={r["Gene_symbol"] for r in clrows if r["Frequency_preference"]=="HIGHER_IN_COREAD"}
audit={
    "timestamp_analysis":"2026-10-03","input_models":1165,"distinct_ids":1165,"coread":63,
    "nonCORE_total":1102,"comparator_source_not_noncancer":1099,
    "explicit_noncancer_excluded":[{"id":i,"name":n} for i,n,diag in zip(new["ids"],new["names"],new["diagnoses"]) if diag=="Non-Cancerous"],
    "auxiliary_nonmodel_columns_1based":[4,5],"all_input_rows_have_1170_fields":True,
    "input_gene_rows":1620,"exact_query_matches":1620,"query_total":1720,"query_missing":100,
    "numeric_entries":int(new["values"].size),"numeric_missing":int(np.isnan(new["values"]).sum()),
    "shared_with_old_genes":1619,"shared_with_old_models":250,"shared_values_checked":404750,
    "shared_value_discrepancies":0,"source_total_average_all1165_within_rounding":True,
    "candidate":335,"higher_FDR05":23,"lower_FDR05":10,"main_test_family":1620,
    "lineage_test_family":int(np.isfinite(joint).sum()),"tested_comparator_lineages":24,
    "balanced_genes":1569,"balanced_boot_draws":10000,"known_pairs_robust_candidates":323,
    "higher_and_balanced_CI_positive":20,"old_significant_lost":sorted(oldset-newset),"new_significant_added":sorted(newset-oldset),
    "integration":integ_summary,"hashes":hashes,
    "software":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"statsmodels":statsmodels.__version__},
    "statistical_validation":"BH checked against statsmodels; selected Wilson and Newcombe CIs checked; shared raw values and all primary organoid counts checked.",
}
json.dump(audit,open(OUT/"analysis_audit.json","w"),indent=2)
print("Change higher:",audit["old_significant_lost"],audit["new_significant_added"])

sources=[
("Cell-line data","User-supplied Texto pegado(1).txt. Probabilities; DepMap release, screen QC and paired cell-line Chronos matrix not supplied. No release was substituted."),
("Query universe","Markdown pegado(20261001-223143).md. 1,720 unique symbols; APC is the first entry."),
("Organoid source","Supplementary_Table_COREAD_1720_Dependency_Candidates.xlsx. Individual Chronos values and archived re-estimated EB probabilities."),
("DepMap threshold","https://forum.depmap.org/t/dependency-threshold/498/2"),
("Fisher exact test","https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html"),
("Conditional odds ratio","https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.contingency.odds_ratio.html"),
("Wilson / Newcombe intervals","https://www.statsmodels.org/stable/generated/statsmodels.stats.proportion.proportion_confint.html ; https://www.statsmodels.org/stable/generated/statsmodels.stats.proportion.confint_proportions_2indep.html"),
("Benjamini-Hochberg","https://doi.org/10.1111/j.2517-6161.1995.tb02031.x"),
("Organoid study","https://doi.org/10.1038/s41586-026-10843-7"),
("Organoid metadata","https://doi.org/10.6084/m9.figshare.29472362.v1"),
("DLD-1 / HCT-15 relationship","https://www.atcc.org/products/ccl-221 ; https://pubmed.ncbi.nlm.nih.gov/7621404/"),
("COLO 201 / COLO 205 relationship","https://www.cellosaurus.org/CVCL_0218"),
]+[(f"{name} SHA-256",sha) for name,sha in hashes.items()]

clnotes=[
("Scope","Full reanalysis of the expanded supplied probability matrix. No global mean, other-cancer preference, organoid or cancer-driver prefilter is used for the primary candidate catalogue."),
("Matrix audit","1,620 source genes × 1,165 distinct ACH identifiers; 1,887,300 entries with 8,748 missing. Read tab delimiters and decimal commas. All source data rows contain 1,170 fields; all finite probabilities lie in [0,1]."),
("Auxiliary columns","Only columns with ACH identifiers in the DepMap_ID row are measurements. The two AAA columns (source columns 4 and 5, 1-based) are auxiliary values and not models or an alternative gene list; they are not used."),
("Gene mapping","All 1,620 source gene symbols now match the 1,720-gene universe exactly. The new source itself supplies CYP2E1 (previously unresolvable CYP 20), ARL2, NPR2, SOS2, SOS1, MCF2 and VNN1. No new alias mapping was inferred."),
("COREAD cohort","Use exact supplied Cancer_detailed = Colorectal Adenocarcinoma: 63 ModelIDs. Model names, diagnoses and linages are preserved as supplied. No fresh screen QC flags or complete donor grouping were available."),
("Primary comparator","1,099 non-COREAD models after excluding three explicitly labelled Non-Cancerous by the source: ACH-001310 (HA1E), ACH-002210 (ARH-77), ACH-002247 (IM-9). This is a source-label-based analysis choice, not an external reassignment of their diagnoses."),
("All-model sensitivity","Repeat the comparison using every non-COREAD model (n=1,102), including the three source-labelled Non-Cancerous models. These are clearly separated All_nonCORE columns; all 1,165 models are preserved in the metadata."),
("Primary candidate","Supplied probability strictly >0.5 in ≥2 distinct COREAD ModelIDs. Exploratory descriptive rule, not an FDR-controlled discovery set. Missing observations are excluded per gene/group, never converted to 0."),
("Descriptive summaries","Counts, finite-data denominator, mean, median, Q25/Q75, sample SD, range and frequency. Mean probability is a summary over cell lines, not a disease-level probability of essentiality."),
("Confidence intervals","Wilson pointwise 95% frequency intervals; Newcombe pointwise 95% difference intervals based on Wilson bounds. Working assumption: independent model observations. No multiplicity adjustment; donor independence has not been established."),
("Main statistical test","Two-sided Fisher exact test for COREAD versus primary-comparator positive/nonpositive calls. Require ≥5 finite observations/group. Benjamini–Hochberg across all 1,620 eligible genes BEFORE candidate filtering. Tables with no call variation retain P=1."),
("Odds ratio","Conditional maximum likelihood odds ratio, with exact conditional 95% confidence limits. Rows: COREAD, comparator; columns: positive, nonpositive. Inf means unbounded. With no call variation, OR is unidentifiable (blank) and CI=[0,Inf]."),
("FDR interpretation","Fisher BH FDR describes differences in frequency between measured groups; it is not a dependency probability, a cell-line Chronos effect, or the FDR of the ≥2-model catalogue. Both higher and lower frequency directions are retained."),
("High-probability sensitivity","Recalculate calls at probability strictly >0.9; repeat Fisher and BH over its own complete family of 1,620 eligible genes. Supporting sensitivity, not independent validation."),
("Recurrence sensitivity","Report ≥12.5%, ≥25%, ≥50%, ≥80% frequencies with gene-specific denominators. The ≥2/63 line rule has a different recurrence floor from ≥2/16 organoids."),
("Per-lineage context","Descriptive counts/frequencies/means for all 27 comparator lineages. Two-sided Fisher for 24 lineages with ≥5 total models, further requiring ≥5 finite observations for each tested gene/group. Joint BH over 38,755 eligible gene-by-lineage tests."),
("Small-lineage limits","Adrenal Gland (1), Ampulla of Vater (4) and Vulva/Vagina (2) are descriptive only. They remain in the primary pooled comparator but not in per-lineage tests or the equal-lineage sensitivity."),
("Balanced comparator","Equal weight to each of the 24 sufficiently represented comparator lineages. Require ≥5 finite values in every included lineage and COREAD: 1,569 eligible genes. No partial-lineage imputation for the 51 ineligible measured genes."),
("Balanced bootstrap","10,000 draws resampling complete model columns independently within COREAD and each fixed lineage; seed 20261003. Missingness retained; recompute finite denominators per draw. Percentile pointwise 95% intervals; all draws valid for all 1,569 eligible genes. No additional P/FDR is manufactured."),
("Balanced endpoints","An endpoint equal to zero may arise in an empirical bootstrap of a no-variation stratum; it is not proof that the population dependency rate is exactly 0 or 1. Inferential use is conditional on sampled models."),
("Leave-one-lineage-out","Recompute the equal-lineage comparator with each of 24 sufficiently represented lineages omitted in turn. Minimum/maximum differences are sensitivity summaries, not confidence limits."),
("Broad-dependency annotation","Frequency ≥80% in COREAD and EACH of the 24 sufficiently represented comparator lineages, only when balanced eligibility is met. Not an official DepMap pan-dependency label and not evidence about normal-tissue toxicity."),
("Related-model sensitivity","Four 61-model COREAD subsets retain one representative of HCT-15/DLD-1 and one of COLO 201/COLO 205 (all 2×2 choices). Recalculate counts/frequencies and Fisher/BH in each; 323 candidates persist in all four. Worst FDR is the maximum of four separately corrected results."),
("Relatedness limitation","Only the two documented COREAD pairs are tested. Other related lines or derivatives may remain in COREAD or the expanded comparator. No comprehensive donor-independence claim is made."),
("Previous-file comparison","All 404,750 gene×model entries shared with the previous 250-model matrix agree exactly, including missingness, for 1,619 mapped genes. COREAD membership is unchanged. The within-COREAD candidate set remains 335; comparator/FDR and coverage are recalculated."),
("Total_average audit","Supplied Total_average is retained but never used for selection. In this expanded matrix, all 1,620 row means over finite values agree within the rounding tolerance 0.0005."),
("Unperformed analyses","Cell-line Chronos, screen QC refitting, MSI/mutation-stratified analyses and complete donor correction are not available from the supplied matrix. Probabilities are not transformed into gene effects; no unmatched release was substituted."),
("Reproducibility","Computed 3 October 2026; Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0. BH and selected CIs independently checked with statsmodels 0.14.6. Simple arithmetic uses workbook formulas; inferential outputs are frozen. Changing inputs or cutoffs requires rerunning all statistics."),
("Source data","All 63 raw COREAD probability columns are included at the right of the gene table. Other-model values remain in the source text, identified by SHA-256; all 1,165 model metadata are included below."),
]
intnotes=[
("Parallel analysis","Both systems are assessed against the original 1,720 symbols, without serial filtering. Union retains candidate evidence in either; intersection indicates cross-system support, not independent validation."),
("Criteria","Lines: supplied probability >0.5 in ≥2 of up to 63 COREAD ModelIDs. Organoids: recovered individual Chronos <-0.5 in ≥2 of 16 3DO models. The scales and recurrence floors are different."),
("Coverage","1,620 line-evaluable genes; 1,636 organoid-evaluable genes; 1,613 measured in both. Seven measured only in lines, 23 only in organoids, 77 in neither. Unavailable is not equivalent to nondependent."),
("Updated integration","335 line candidates, 324 organoid candidates, 193 in both and 466 in their observed union. Among jointly measured genes: 139 line-only, 130 organoid-only, 1,151 neither."),
("Unavailable cross-system candidates","C1QTNF9, YWHAE and UQCRQ are line candidates without organoid measurements. PRKDC is an organoid candidate without a line row. These four genes are not classified as discordant measured results."),
("Organoid provenance","Reused Supplementary_Table_COREAD_1720_Dependency_Candidates.xlsx and recomputed all 16-model counts and flags. No change to organoid matrix, scoring or model selection. No use of the obsolete 37-gene aggregate-bound inference."),
("Supporting probabilities","Organoid empirical-Bayesian probabilities are archived re-estimates, not the published DepMap probability matrix. Same-model concordance requires Chronos <-0.5 AND probability >0.5 in the same ≥2 models; no new probability model was fitted for this update."),
("No pooling","Do not average line probabilities with re-estimated organoid probabilities. No pooled probability, combined P/FDR or joint Chronos effect is calculated."),
("Recurrence alignment","Sensitivity requires ≥12.5% recurrence and ≥2 models in each system: 172 line candidates and 148 jointly supported genes. This aligns the recurrence floor, not the outcome definitions."),
("Higher recurrence","Shared recurrence ≥25%, ≥50%, ≥80% is reported descriptively: 110, 85 and 43 genes respectively. These nested groups are not a new FDR-controlled selection."),
("Known-pair sensitivity","189 genes remain candidates in both systems under all four choices of representatives of the two documented related COREAD line pairs. This is a limited relationship audit."),
("Model identity overlap","No shared ModelIDs between any of the 1,165 supplied lines and the 16 primary / 18 expanded organoid models listed below. Donor overlap is unassessed."),
("Rank concordance","Spearman correlation of line probability-call frequency and organoid Chronos-call frequency over ALL 1,613 jointly measured genes (no candidate prefilter). Cross-criterion rank association, not equality of gene effects or assays."),
("Rank uncertainty","2,000 bootstrap draws of complete model columns: resample 63 lines and 16 organoids independently; seed 20261004. Recompute finite-value denominators and Spearman correlation. Pointwise percentile 95% CI; gene-dependence is preserved within each resampled model."),
("Nonbroad sensitivity","Recalculate rank association on the FIXED set of 1,570 genes excluding 43 with ≥80% recurrence in both systems. Not a new gene-discovery test."),
("Overlap summaries","Jaccard and Dice use only jointly evaluable genes (candidate union n=462 on this shared universe); the full observed union is 466. No random-set overlap P value is assigned."),
("Separate statistical contexts","Cell-line Fisher/BH columns now use 63 COREAD versus 1,099 other source-labelled cancer models. The all-1,102 sensitivity is separate. Organoid Welch/BH annotations retain the original 16-vs-76 comparison. Neither comparison filters this integration."),
("Source uncertainty","Cell-line DepMap release and paired Chronos matrix are absent. No matched Chronos–Chronos comparison, causal culture-format effect, normal-tissue window, microbiome mechanism or independent biological validation is established."),
("Display","All 1,720 symbols are retained. Filter Union_candidate=YES for 466 or Both_systems_candidate=YES for 193. NOT_EVALUABLE distinguishes missing cross-system information."),
]
print("Notes prepared",len(clnotes),len(intnotes))