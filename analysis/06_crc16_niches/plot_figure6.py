"""Editable six-panel figure and exact plotting tables for current Figure 6."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from analyze_crc16 import C9, NICHES

GENE_ORDER = ["MYC","FASN","SOX9","BCL2L1","IRS2","CTNNB1","HMGB1","RNF6","YAP1"]
AXES = ["Immune_axis","Metabolic_axis","Endocrine_axis","Neural_axis","Stromal_matrix_axis"]
AXIS_LABELS = ["Immune","Metabolic","Endocrine","Neural","Stromal–matrix"]
SELECTED = ["N3","N4","N6","N9"]
PROGRAMMES = ["N4_epithelial_antimicrobial","N4_N9_chemokine_myeloid","N9_CXCL8_PTGS2_inflammatory",
              "N6_Convergent9","N6_WNT_stemness","N3_effector_axis","N3_cGAS_STING_IFN","MHCI"]
PROGRAMME_LABELS = ["N4 epithelial antimicrobial","N4/N9 chemokine–myeloid","N9 CXCL8–PTGS2 inflammatory",
    "N6 Convergent9","N6 WNT–stemness","N3 effector","N3 cGAS–STING/IFN","MHC-I"]
LINEAGES = ["Epi","TNKILC","Myeloid","B","Plasma","Strom","Mast"]
LINEAGE_LABELS = ["Epithelial","T/NK/ILC","Myeloid","B cells","Plasma cells","Stromal","Mast cells"]
IMMUNE_LABELS = ["cGAS–STING/IFN","CD8 recruitment /\ncytotoxicity","Antigen presentation /\nMHC-I",
                "PD-1/PD-L1 checkpoint","CXCL8–CXCR2/myeloid","Succinate/HIF/EZH2","TBX21 effector"]
IMMUNE_ORDER = ["cGAS_STING_IFN","CD8_recruitment_cytotoxicity","Antigen_presentation_MHCI",
                "PD1_PDL1_checkpoint","CXCL8_CXCR2_myeloid","Succinate_HIF_EZH2","TBX21_effector"]


def draw(base, out):
    base, out = Path(base), Path(out)
    spots = pd.read_csv(out/"CRC16_spot_analysis_metadata.csv.gz")
    effects = pd.read_csv(out/"Spatial43_Convergent9_gene_effects.csv")
    annotations = pd.read_csv(base/"input/Spatial43_axis_annotations.csv")
    burden = pd.read_csv(out/"bacterial_burden_by_niche.csv").set_index("Niche").loc[NICHES].reset_index()
    taxa = pd.read_csv(out/"taxon_niche_models_publication_precision.csv")
    immune = pd.read_csv(out/"immune_programmes_by_niche.csv")
    attribution = pd.read_csv(out/"GSE178341_Top_programme_summaries.csv")
    taxa["Displayed"] = np.isfinite(taxa.OR_adjusted)&(taxa.Positive_niche>=5)&(taxa.Positive_rest>=5)
    taxa["Adjusted_log2_OR"] = np.log2(taxa.OR_adjusted)
    taxa["Displayed_log2_OR"] = taxa.Adjusted_log2_OR.clip(-3,3).where(taxa.Displayed)
    if int(taxa.Displayed.sum()) != 98:
        raise ValueError("Expected 98 displayed and 162 masked taxon contrasts")
    spotplot=spots[["barcode","Niche","imagecol","imagerow","svg_color"]]
    spotplot.to_csv(out/"Figure6A_spatial_niches.csv",index=False)
    effects.to_csv(out/"Figure6B_gene_expression.csv",index=False)
    annotations.to_csv(out/"Figure6B_axis_annotations.csv",index=False)
    burden.to_csv(out/"Figure6C_bacterial_occupancy_burden.csv",index=False)
    taxa.to_csv(out/"Figure6D_taxon_niche_screen.csv",index=False)
    immune[immune.Niche.isin(SELECTED)].to_csv(out/"Figure6E_immune_context.csv",index=False)
    attribution.to_csv(out/"Figure6F_cellular_attribution.csv",index=False)

    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"svg.fonttype":"none",
        "axes.spines.top":False,"axes.spines.right":False,"savefig.facecolor":"white"})
    fig=plt.figure(figsize=(12,15))
    a=fig.add_axes([.07,.65,.30,.28])
    b=fig.add_axes([.52,.60,.32,.35])
    ba=fig.add_axes([.865,.60,.075,.35])
    c=fig.add_axes([.075,.365,.30,.16])
    d=fig.add_axes([.62,.34,.32,.19])
    e=fig.add_axes([.28,.08,.17,.17])
    f=fig.add_axes([.69,.08,.25,.17])
    fig.text(.05,.967,"A  Host-defined spatial niches",fontweight="bold",fontsize=11)
    fig.text(.05,.950,"CRC_16 · 3,106 locations · 13 frozen niches",fontsize=8,color="#555555")
    fig.text(.47,.975,"B  Spatial43 expression across niches",fontweight="bold",fontsize=11)
    fig.text(.47,.960,"* Convergent9 membership",fontsize=8,color="#555555")
    colors=spotplot.drop_duplicates("Niche").set_index("Niche").svg_color.to_dict()
    a.scatter(spots.imagecol,spots.imagerow,c=spots.svg_color,s=4,linewidths=0)
    for ni in NICHES:
        subset=spots[spots.Niche.eq(ni)]
        a.text(subset.imagecol.median(),subset.imagerow.median(),ni,ha="center",va="center",fontsize=7,fontweight="bold")
    a.invert_yaxis();a.set_aspect("equal");a.set_xlabel("Spatial x coordinate (pixels)");a.set_ylabel("Spatial y coordinate (pixels)")
    a.legend(handles=[Patch(color=colors[n],label=n) for n in NICHES],loc="lower center",bbox_to_anchor=(.5,-.25),ncol=7,
             frameon=False,fontsize=7,columnspacing=.7,handlelength=.8)
    gene_order=GENE_ORDER+sorted(set(effects.Gene)-set(GENE_ORDER))
    ge=effects.pivot(index="Gene",columns="Niche",values="Effect").loc[gene_order,NICHES]
    im=b.imshow(ge,cmap="RdBu_r",vmin=-4,vmax=4,aspect="auto",interpolation="none")
    b.set_yticks(range(len(gene_order)),[g+(" *" if g in C9 else "") for g in gene_order],fontsize=6.5)
    b.set_xticks(range(13),NICHES,rotation=90,fontsize=7);b.tick_params(length=0);b.spines[:].set_visible(False)
    binary=annotations.set_index("Gene").loc[gene_order,AXES].astype(int)
    ba.imshow(binary,cmap=ListedColormap(["#f4f4f4","#444444"]),vmin=0,vmax=1,aspect="auto",interpolation="none")
    ba.set_xticks(range(5),AXIS_LABELS,rotation=90,fontsize=6.5)
    ba.set_yticks([]);ba.tick_params(length=0);ba.spines[:].set_visible(False)
    colorbar=fig.colorbar(im,cax=fig.add_axes([.954,.64,.007,.28]),orientation="vertical")
    colorbar.set_label("Expression effect: niche − rest, log₂(CPM + 1)",fontsize=7);colorbar.ax.tick_params(labelsize=7)

    fig.text(.05,.556,"C  Bacterial occupancy and burden",fontweight="bold",fontsize=11)
    c.bar(range(13),burden.Bacteria_positive_fraction*100,color=[colors[n] for n in NICHES],width=.8)
    c.set_xticks(range(13),NICHES,rotation=45,fontsize=7);c.set_ylabel("Bacteria-positive locations (%)");c.set_xlabel("Host-defined niche")
    twin=c.twinx();twin.plot(range(13),burden.Mean_filtered_bacterial_UMI,"--o",color="#222222",lw=1,ms=3)
    twin.set_ylabel("Mean bacterial UMIs / location");twin.tick_params(labelsize=7)
    c.text(.02,1.04,"Bars: occupancy   Dashed line: mean burden",transform=c.transAxes,fontsize=7)
    fig.text(.50,.556,"D  Bacterial taxon–niche associations",fontweight="bold",fontsize=11)
    qc=pd.read_csv(out/"bacteria_qc_taxa.csv")
    taxon_order=qc.loc[qc.eligible,"Taxon"].tolist()
    heat=taxa.pivot(index="Taxon",columns="Niche",values="Displayed_log2_OR").loc[taxon_order,NICHES]
    palette=plt.get_cmap("RdBu_r").copy();palette.set_bad("#d9d9d9")
    imd=d.imshow(heat,cmap=palette,vmin=-3,vmax=3,aspect="auto",interpolation="none")
    d.set_yticks(range(20),[x.replace("unclassified_Erysipelotrichaceae_.miscellaneous.","Unclassified Erysipelotrichaceae") for x in taxon_order],fontsize=7)
    d.set_xticks(range(13),NICHES,rotation=90,fontsize=7);d.tick_params(length=0);d.spines[:].set_visible(False)
    for _,r in taxa[taxa.Robust_primary & taxa.Displayed].iterrows():
        d.text(NICHES.index(r.Niche),taxon_order.index(r.Taxon),"*",ha="center",va="center",fontsize=13)
    cb=fig.colorbar(imd,cax=fig.add_axes([.954,.34,.007,.19]),orientation="vertical",ticks=[-3,0,3])
    cb.set_label("Adjusted log₂ odds ratio; grey = masked",fontsize=7);cb.ax.tick_params(labelsize=7)

    fig.text(.05,.278,"E  Immune transcriptional context",fontweight="bold",fontsize=11)
    ie=immune.pivot(index="Axis",columns="Niche",values="Effect").loc[IMMUNE_ORDER,SELECTED]
    ime=e.imshow(ie,cmap="RdBu_r",vmin=-.45,vmax=.45,aspect="auto",interpolation="none")
    e.set_yticks(range(7),IMMUNE_LABELS,fontsize=7);e.set_xticks(range(4),SELECTED);e.tick_params(length=0);e.spines[:].set_visible(False)
    for i in range(7):
        for j in range(4):e.text(j,i,f"{ie.iloc[i,j]:.2f}",ha="center",va="center",fontsize=6,color="white" if abs(ie.iloc[i,j])>.30 else "#222222")
    cb=fig.colorbar(ime,cax=fig.add_axes([.28,.043,.17,.007]),orientation="horizontal",ticks=[-.4,0,.4])
    cb.set_label("Programme-score effect: niche − rest",fontsize=7);cb.ax.tick_params(labelsize=7)
    fig.text(.53,.278,"F  Single-cell attribution of programmes",fontweight="bold",fontsize=11)
    at=attribution.pivot(index="Signature",columns="Group",values="Median_score").loc[PROGRAMMES,LINEAGES]
    imf=f.imshow(at,cmap="RdBu_r",vmin=-1.5,vmax=1.5,aspect="auto",interpolation="none")
    npat=attribution.groupby("Group").N_patients.first()
    f.set_xticks(range(7),[f"{label}\n(n = {int(npat[n])})" for label,n in zip(LINEAGE_LABELS,LINEAGES)],rotation=45,ha="right",fontsize=7)
    f.set_yticks(range(8),PROGRAMME_LABELS,fontsize=7);f.tick_params(length=0);f.spines[:].set_visible(False)
    cb=fig.colorbar(imf,cax=fig.add_axes([.69,.027,.25,.007]),orientation="horizontal",ticks=[-1.5,0,1.5])
    cb.set_label("Median patient-level programme score",fontsize=7);cb.ax.tick_params(labelsize=7)
    fig.savefig(out/"Figure6.svg",metadata={"Title":"Figure 6: current Spatial43/Convergent9 CRC16 analysis"})
    fig.savefig(out/"Figure6.png",dpi=180)
    plt.close(fig)
    return {"Figure6B_gene_niche_cells":len(effects),"Figure6D_displayed":int(taxa.Displayed.sum()),
            "Figure6D_masked":int((~taxa.Displayed).sum()),"Figure6E_cells":len(immune[immune.Niche.isin(SELECTED)]),
            "Figure6F_cells":len(attribution)}
