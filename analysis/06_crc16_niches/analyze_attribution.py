"""Counts-based Convergent9 attribution plus archived programme checkpoints."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from analyze_crc16 import C9, zscores, check_numeric
from stats_helpers import bh


def run(base, out):
    base, out = Path(base), Path(out)
    data, ref = base/"input", base/"reference"
    cell = np.load(data/"GSE178341_C9_cell_counts.npz")
    counts = cell["counts"]
    library = cell["library_size"]
    genes = cell["genes"].tolist()
    barcode = cell["barcodes"].astype(str)
    if set(genes) != set(C9) or counts.shape != (370115, 9):
        raise ValueError("Unexpected C9 count checkpoint")
    metadata = pd.read_csv(data/"GSE178341_source_annotations.csv.gz").set_index("cellID").loc[barcode]
    is_tumour = metadata.SPECIMEN_TYPE.eq("T").to_numpy()
    annotation_columns = {"Top":"clTopLevel", "Mid":"clMidwayPr", "Fine":"cl295v11SubFull"}
    patients, rankings, all_summaries, pseudobulks = [], [], [], []
    for level, column in annotation_columns.items():
        identities = pd.MultiIndex.from_arrays([metadata.loc[is_tumour,"PID"], metadata.loc[is_tumour,column]], names=["PID","Group"])
        codes, groups = pd.factorize(identities, sort=True)
        groups = groups.set_names(["PID", "Group"])
        c = np.zeros((len(groups),9),dtype=np.int64)
        total = np.zeros(len(groups),dtype=np.int64)
        ncells = np.bincount(codes,minlength=len(groups))
        np.add.at(c,codes,counts[is_tumour])
        np.add.at(total,codes,library[is_tumour])
        keep = ncells >= 20
        c,total,ncells,groups = c[keep],total[keep],ncells[keep],groups[keep]
        ex = np.log2(c * 1e6 / total[:,None] + 1)
        score = zscores(ex).mean(1)
        frame = pd.DataFrame({"PID":groups.get_level_values("PID"),"Group":groups.get_level_values("Group"),
            "N_cells":ncells,"Signature":"N6_Convergent9","Score":score,"N_genes":9,"Level":level})
        patients.append(frame)
        bulk = frame[["Level","PID","Group","N_cells"]].copy()
        bulk["Full_library_size"] = total
        for j,gene in enumerate(genes):
            bulk[gene+"_count"] = c[:,j]
            bulk[gene+"_log2CPM1"] = ex[:,j]
        pseudobulks.append(bulk)
        summary = frame.groupby("Group").agg(N_patients=("PID","nunique"),Median_score=("Score","median"),
                                             Mean_score=("Score","mean"),Median_cells=("N_cells","median")).reset_index()
        summary = summary[summary.N_patients >= 5].sort_values(["Median_score","Group"],ascending=[False,True])
        summary.insert(0,"Level",level);summary.insert(1,"Signature","N6_Convergent9")
        summary.insert(2,"Rank",np.arange(1,len(summary)+1))
        all_summaries.append(summary)
        rankings.append(summary if level=="Top" else summary.head(10))
    c9_patient=pd.concat(patients,ignore_index=True)
    c9_rank=pd.concat(rankings,ignore_index=True)
    rank_ref=pd.read_csv(ref/"lineage_rankings.csv")
    rank_error=check_numeric(c9_rank,rank_ref[rank_ref.Signature.eq("N6_Convergent9")],
                             ["Level","Signature","Group"],["N_patients","Median_score","Mean_score","Median_cells"])
    c9_patient.to_csv(out/"GSE178341_C9_patient_population_scores.csv",index=False)
    pd.concat(pseudobulks,ignore_index=True).to_csv(out/"GSE178341_C9_pseudobulk_counts_expression.csv.gz",index=False)
    pd.concat(all_summaries,ignore_index=True).to_csv(out/"GSE178341_C9_all_eligible_population_summaries.csv",index=False)
    # The seven other programmes are the archived patient-score input. Their
    # original gene lists were not recovered; no replacement list is inferred.
    checkpoint=pd.read_csv(data/"GSE_Patient_TopScores.csv")
    checkpoint=checkpoint[~checkpoint.Signature.eq("N6_Convergent6")].copy()
    checkpoint["Level"]="Top"
    checkpoint["Calculation_origin"]="Archived patient-score checkpoint; gene definitions unavailable"
    fresh=c9_patient[c9_patient.Level.eq("Top")].copy()
    fresh["Calculation_origin"]="Recomputed from nine-gene cell counts and complete cell libraries"
    top=pd.concat([checkpoint,fresh],ignore_index=True)
    if top.Signature.nunique()!=8 or top.PID.nunique()!=62:
        raise ValueError("Expected eight programmes and 62 tumour patients")
    top.to_csv(out/"GSE178341_Top_patient_scores_current.csv",index=False)
    broad=top.groupby(["Group","Signature"]).agg(N_patients=("PID","nunique"),Median_score=("Score","median"),
                                                 Mean_score=("Score","mean"),Median_cells=("N_cells","median")).reset_index()
    broad=broad[broad.N_patients>=5]
    broad.to_csv(out/"GSE178341_Top_programme_summaries.csv",index=False)
    broad_ref=rank_ref[rank_ref.Level.eq("Top")]
    broad_error=check_numeric(broad,broad_ref,["Group","Signature"],["N_patients","Median_score","Mean_score","Median_cells"])
    targets={"N4_epithelial_antimicrobial":"Epi","N6_Convergent9":"Epi","N6_WNT_stemness":"Epi",
             "N4_N9_chemokine_myeloid":"Myeloid","N9_CXCL8_PTGS2_inflammatory":"Myeloid","N3_effector_axis":"TNKILC"}
    tests,paired=[],[]
    for signature,target in targets.items():
        wide=top[top.Signature.eq(signature)].pivot(index="PID",columns="Group",values="Score")
        other=wide.drop(columns=target).mean(axis=1)
        difference=(wide[target]-other).dropna()
        probability=wilcoxon(difference.to_numpy(),alternative="two-sided").pvalue
        tests.append({"Signature":signature,"Target_lineage":target,"N_patients":len(difference),
            "Median_target_minus_rest":difference.median(),"Mean_target_minus_rest":difference.mean(),
            "Wilcoxon_P":probability,"Fraction_target_gt_rest":(difference>0).mean()})
        paired.extend({"Signature":signature,"Target_lineage":target,"PID":pid,"Target_score":wide.loc[pid,target],
                       "Mean_other_lineages":other.loc[pid],"Target_minus_rest":value} for pid,value in difference.items())
    tests=pd.DataFrame(tests);tests["BH_FDR"]=bh(tests.Wilcoxon_P,6)
    lineage_error=check_numeric(tests,pd.read_csv(ref/"lineage_tests.csv"),["Signature","Target_lineage"],
        ["N_patients","Median_target_minus_rest","Mean_target_minus_rest","Wilcoxon_P","Fraction_target_gt_rest","BH_FDR"])
    tests.to_csv(out/"GSE178341_lineage_contrasts.csv",index=False)
    pd.DataFrame(paired).to_csv(out/"GSE178341_lineage_patient_differences.csv",index=False)
    # Mid/Fine older programme summaries also remain explicitly archived.
    old=pd.read_csv(base/"provenance/top_cellular_attribution.csv")
    old=old[~old.Signature.eq("N6_Convergent6")]
    current_rank=pd.concat([old,c9_rank],ignore_index=True)
    current_rank.to_csv(out/"GSE178341_population_rankings_current.csv",index=False)
    report={"source_cells":len(barcode),"full_feature_library_UMI":int(library.sum()),"tumour_patients":int(top.PID.nunique()),
            "gene_count":9,"patient_profile_min_cells":20,"population_min_patients":5,"programme_count":8,
            "paired_tests":6,"C9_rank_max_errors":rank_error,"Top_summary_max_errors":broad_error,"lineage_test_max_errors":lineage_error,
            "recomputed_from_counts":["N6_Convergent9"],"checkpoint_programmes":sorted(checkpoint.Signature.unique().tolist()),
            "checkpoint_limitation":"Original gene definitions for seven historical attribution programmes were not recovered; patient scores, summaries and paired tests remain reproducible from archived scores."}
    (out/"GSE178341_attribution_validation.json").write_text(json.dumps(report,indent=2)+"\n")
    return report
