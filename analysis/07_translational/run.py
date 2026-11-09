#!/usr/bin/env python3
"""Reproduce deterministic Figure 7 tables from the frozen ST14 appraisal.

Evidence records were manually curated. This program validates those records,
derives counts/coverage/links, and redraws Figure 7; it does not infer efficacy.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import tempfile
import textwrap
os.environ.setdefault("MPLCONFIGDIR",str(Path(tempfile.gettempdir())/"coread_matplotlib"))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.path import Path as MplPath
from matplotlib.colors import ListedColormap, BoundaryNorm

CLASSES = ["A", "B", "C", "D", "E"]
CLASS_COLORS = {"A":"#2166ac", "B":"#1b9e77", "C":"#e6ab02", "D":"#d95f02", "E":"#7570b3"}
NICHE_COLORS = {"N3":"#443983", "N4":"#3b528b", "N6":"#287c8e", "N9":"#35b779"}
NICHES = list(NICHE_COLORS)
MODULE_LABELS = ["Approved biomarker\nbackbone", "Fusobacterium\nreduction",
    "CXCR1/2 +\nimmunotherapy", "BCL-xL\ndegradation", "FASN\ninhibition",
    "WNT/β-catenin\ninhibition", "MYC\ninhibition", "Precision\nantimicrobial"]
STRATEGY_LABELS = ["Nivolumab + ipilimumab", "Encorafenib + cetuximab +\nfluorouracil-based chemotherapy",
    "Sotorasib + panitumumab", "Metronidazole-based reduction", "SX-682 + PD-1 immunotherapy",
    "DT2216", "Denifanstat / TVB-2640", "E7386", "OMO-103 / Omomyc", "Precision antimicrobial /\nnanomedicine"]


def plot_figure(tables, out):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7,"svg.fonttype":"none",
                         "axes.spines.top":False,"axes.spines.right":False,"savefig.facecolor":"white"})
    fig = plt.figure(figsize=(11, 12.5))
    axa = fig.add_axes([.035,.625,.36,.29])
    axb = fig.add_axes([.70,.625,.26,.29])
    axc = fig.add_axes([.34,.30,.62,.235])
    axd = fig.add_axes([.055,.035,.90,.18])
    fig.text(.04,.947,"A  Niche–module correspondence",fontweight="bold",fontsize=10)
    coords_n = {n:np.array([np.cos(t),np.sin(t)]) for n,t in zip(NICHES,np.deg2rad([145,175,205,235]))}
    theta = np.deg2rad([70,50,30,10,-10,-30,-50,-70])
    coords_m = {i:np.array([np.cos(t),np.sin(t)]) for i,t in enumerate(theta)}
    coverage = tables["coverage"]
    for i, row in coverage.iterrows():
        for n in NICHES:
            code = int(row[n])
            if not code: continue
            start, end = coords_n[n], coords_m[i]
            path = MplPath([start, start*.2, end*.2, end], [1,4,4,4])
            axa.add_patch(patches.PathPatch(path, fill=False, lw=1.0,
                color=NICHE_COLORS[n] if code==2 else "#aaaaaa", linestyle="-" if code==2 else (0,(3,2)), alpha=.8))
    for n, point in coords_n.items():
        axa.scatter(*point, s=70, color=NICHE_COLORS[n], zorder=4)
        axa.text(point[0]-.11,point[1],n,ha="right",va="center",fontweight="bold")
    for i, point in coords_m.items():
        cls = "A" if i==0 else "D" if i==7 else "B"
        axa.scatter(*point,s=55,color=CLASS_COLORS[cls],zorder=4)
        label_y = np.linspace(1.13,-1.13,8)[i]
        axa.plot([point[0]+.05,1.15],[point[1],label_y],color="#cccccc",lw=.5)
        axa.text(1.20,label_y,MODULE_LABELS[i],ha="left",va="center",fontsize=7.5)
    axa.set(xlim=(-1.5,2.7),ylim=(-1.4,1.35),aspect="equal");axa.axis("off")
    axa.text(-1.45,-1.18,"8 direct links  •  5 conditional links\nEqual widths; no efficacy encoding",fontsize=6.2,va="top")

    evidence = tables["strategy_evidence"]
    counts = evidence["Evidence class"].value_counts()
    display_classes = CLASSES[::-1]
    positions = {c:i for i,c in enumerate(display_classes)}
    fig.text(.52,.947,"B  Evidence for ten strategy records",fontweight="bold",fontsize=10)
    for j,row in evidence.iterrows():
        cls = row["Evidence class"]
        axb.plot(positions[cls],j,marker={"A":"s","B":"o","C":"^","D":"D","E":"v"}[cls],
                 markersize=6,ls="",color=CLASS_COLORS[cls])
    axb.set_yticks(range(10),STRATEGY_LABELS,fontsize=8)
    axb.set_xticks(range(5),[f"{c}\nn = {int(counts.get(c,0))}" for c in display_classes],fontsize=7)
    axb.xaxis.tick_top();axb.tick_params(axis="both",length=0,pad=5)
    axb.set_xlim(-.5,4.5);axb.set_ylim(9.65,-.7);axb.grid(axis="x",alpha=.15)
    axb.spines[:].set_visible(False)

    fig.text(.055,.572,"C  Deterministic spatial coverage",fontweight="bold",fontsize=10)
    matrix = coverage[NICHES].to_numpy(int)
    cmap=ListedColormap(["#f2f2f2","#b8d3df","#226985"])
    axc.imshow(matrix,cmap=cmap,norm=BoundaryNorm([-.5,.5,1.5,2.5],3),aspect="auto")
    labels=[x.replace("\n"," ") for x in MODULE_LABELS]
    axc.set_yticks(range(8),labels,fontsize=8)
    axc.set_xticks(range(4),NICHES)
    axc.xaxis.tick_top();axc.tick_params(length=0)
    for i in range(8):
        for j in range(4):axc.text(j,i,str(matrix[i,j]),ha="center",va="center",color="white" if matrix[i,j]==2 else "#333333",fontsize=8)
    axc.spines[:].set_visible(False)
    axc.legend(handles=[patches.Patch(color=cmap(k),label=l) for k,l in enumerate(["0: no direct rationale","1: conditional","2: direct match"])],
               loc="lower center",bbox_to_anchor=(.5,-.2),ncol=3,frameon=False,fontsize=6.5)

    axd.set_title("D  Biomarker-gated evaluation framework",loc="left",fontweight="bold",fontsize=9,pad=10)
    axd.set(xlim=(-.02,1.02),ylim=(0,1));axd.axis("off")
    labels=[("Molecular eligibility", "MSI-H/dMMR\nBRAF V600E\nKRAS G12C"),
            ("Approved backbone examples", "Nivolumab + ipilimumab\nEncorafenib + cetuximab + chemotherapy\nSotorasib + panitumumab"),
            ("Investigational niche modules", "N3: microbial reduction\nN4/N9: myeloid remodelling\nN6: BCL-xL / FASN / WNT / MYC")]
    widths=[.22,.40,.30];starts=[0,.25,.70]
    for (title, body),left,width in zip(labels,starts,widths):
        axd.add_patch(patches.FancyBboxPatch((left,.40),width,.49,boxstyle="round,pad=0.01,rounding_size=.02",facecolor="#f4f7f8",edgecolor="#91a4ad",lw=.8))
        axd.text(left+width/2,.79,title,ha="center",fontweight="bold",fontsize=6.6)
        axd.text(left+width/2,.59,body,ha="center",va="center",fontsize=6.2,linespacing=1.65)
    for x in [.235,.675]:axd.annotate("",xy=(x+.013,.64),xytext=(x-.013,.64),arrowprops={"arrowstyle":"->","color":"#555555","lw":.8})
    axd.text(.5,.24,"Prospective evaluation: target engagement • safety • antitumour activity",ha="center",fontsize=7,fontweight="bold")
    axd.text(.5,.07,"Indication-specific eligibility is retained. Spatial correspondence generates hypotheses;\nit does not establish drug sensitivity, benefit or synergy.",ha="center",fontsize=6.6,linespacing=1.5)
    fig.suptitle("Spatially informed therapeutic prioritisation",fontsize=14,fontweight="bold",y=.992)
    fig.text(.5,.968,"Frozen appraisal: evidence through 14 October 2026; targeted checks 28 October 2026",ha="center",fontsize=8,color="#555555")
    fig.savefig(out / "Figure7.svg",metadata={"Title":"Figure 7: deterministic therapeutic evidence mapping"})
    fig.savefig(out / "Figure7.png",dpi=220)
    plt.close(fig)


def run(output_dir):
    base=Path(__file__).resolve().parent
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    tables={p.stem:pd.read_csv(p) for p in sorted((base/"input").glob("*.csv"))}
    evidence=tables["strategy_evidence"];coverage=tables["coverage"]
    if len(evidence)!=10 or evidence.Strategy.duplicated().any():raise ValueError("Expected ten unique strategies")
    if set(evidence["Evidence class"])-set(CLASSES):raise ValueError("Unknown evidence class")
    counts=evidence["Evidence class"].value_counts().reindex(CLASSES,fill_value=0)
    if counts.to_dict()!={"A":3,"B":6,"C":0,"D":1,"E":0}:raise ValueError("Evidence class counts disagree with ST14")
    matrix=coverage[NICHES].to_numpy(float)
    if matrix.shape!=(8,4) or not np.isin(matrix,[0,1,2]).all():raise ValueError("Invalid coverage grid")
    if not np.array_equal((matrix==2).sum(1),coverage["Direct-niche count"].to_numpy(int)):raise ValueError("Direct counts disagree with coverage")
    value_counts={int(v):int((matrix==v).sum()) for v in [0,1,2]}
    if value_counts!={0:19,1:5,2:8}:raise ValueError("Frozen coverage counts disagree")
    links=coverage.melt(id_vars=["Therapeutic module"],value_vars=NICHES,var_name="Niche",value_name="Coverage_code")
    links["Coverage_code"]=links.Coverage_code.astype(int)
    edges=links[links.Coverage_code>0].copy()
    edges["Match"]=edges.Coverage_code.map({1:"Conditional/secondary",2:"Direct mechanistic"})
    edges["Line_width"]=1.0
    edges["Line_style"]=edges.Coverage_code.map({1:"dashed",2:"solid"})
    sources=[]
    for _,r in evidence.iterrows():
        for col in ["Source 1","Source 2"]:
            url=r[col]
            if pd.notna(url):
                if not str(url).startswith("https://"):raise ValueError("Invalid primary reference URL")
                sources.append({"Record_type":"Strategy","Record":r.Strategy,"Source_field":col,"URL":url})
    for _,r in tables["limiting_evidence"].iterrows():
        sources.append({"Record_type":"Limiting evidence","Record":r.Approach,"Source_field":"Source","URL":r.Source})
    plot_data=evidence.copy();plot_data["Plot_order"]=np.arange(1,11)
    plot_data["Marker"]=plot_data["Evidence class"].map({"A":"square","B":"circle","D":"diamond"})
    outputs={"Figure7A_niche_module_edges":edges,"Figure7B_strategy_evidence":evidence,
             "Figure7B_strategy_plot_data":plot_data,
             "Figure7B_evidence_class_counts":counts.rename_axis("Evidence_class").rename("N_strategy_records").reset_index(),
             "Figure7C_spatial_coverage":coverage,"Figure7C_coverage_long":links,
             "Figure7C_coverage_counts":pd.DataFrame({"Coverage_code":list(value_counts),"N_cells":list(value_counts.values())}),
             "Figure7D_biomarker_modules":tables["biomarker_modules"],
             "Figure7_evidence_hierarchy":tables["evidence_hierarchy"],
             "Figure7_interpretive_boundaries":tables["interpretive_boundaries"],
             "Figure7_limiting_evidence":tables["limiting_evidence"],"Figure7_niche_states":tables["niche_states"],
             "Figure7_niche_quantitative_anchors":tables["niche_quantitative_anchors"],
             "Figure7_reference_links":pd.DataFrame(sources)}
    for name,df in outputs.items():df.to_csv(out/(name+".csv"),index=False)
    plot_figure(tables,out)
    provenance={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((base/"input").glob("*.csv"))}
    status={"module":"07_translational","mode":"processed","passed":True,"strategy_records":10,
            "class_counts":counts.to_dict(),"coverage_cells":32,"coverage_counts":value_counts,
            "nonzero_links":len(edges),"all_link_widths_equal":bool(edges.Line_width.nunique()==1),
            "inference":"None: deterministic summary of human-curated frozen evidence",
            "evidence_cutoff":"2026-10-14","targeted_checks":"2026-10-28","input_sha256":provenance}
    (out/"status.json").write_text(json.dumps(status,indent=2)+"\n")
    print(json.dumps(status,indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode",choices=["processed"],default="processed")
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args();run(a.output_dir)
