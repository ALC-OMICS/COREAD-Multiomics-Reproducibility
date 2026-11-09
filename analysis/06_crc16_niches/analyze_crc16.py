"""Recompute host effects, bacterial QC, regressions and revised CMH blocks."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from threadpoolctl import threadpool_limits
from stats_helpers import bh, logistic, ols, cmh

C9 = ["BCL2L1", "CTNNB1", "FASN", "HMGB1", "IRS2", "MYC", "RNF6", "SOX9", "YAP1"]
NICHES = [f"N{i}" for i in range(1, 14)]


def expression(host, genes, library, wanted):
    aliases = {"EPRS1": "EPRS"}
    index = {g: i for i, g in enumerate(genes)}
    resolved = [aliases.get(g, g) for g in wanted]
    if any(g not in index for g in resolved):
        raise ValueError("Unavailable requested CRC16 gene")
    counts = host[[index[g] for g in resolved]].toarray().T
    return counts, np.log2(counts * 1e6 / library[:, None] + 1)


def zscores(x):
    sd = x.std(axis=0, ddof=0)
    return np.divide(x - x.mean(axis=0), sd, out=np.zeros_like(x), where=sd > 0)


def check_numeric(actual, reference, keys, columns, tolerance=1e-8):
    a = actual.set_index(keys).sort_index()
    r = reference.set_index(keys).sort_index()
    if not a.index.equals(r.index):
        raise ValueError(f"Reference keys differ: {keys}")
    report = {}
    for col in columns:
        av, rv = a[col].to_numpy(float), r[col].to_numpy(float)
        if not np.array_equal(np.isfinite(av), np.isfinite(rv)):
            raise ValueError(f"Finite-mask mismatch for {col}")
        finite = np.isfinite(av)
        error = float(np.max(np.abs(av[finite] - rv[finite]))) if finite.any() else 0.
        if not np.allclose(av[finite], rv[finite], atol=tolerance, rtol=tolerance):
            raise ValueError(f"Reference disagreement {col}: max absolute {error}")
        report[col] = error
    return report


def run(base, out):
    base, out = Path(base), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    data, ref, archived = base / "input", base / "reference", base / "provenance"
    host = sparse.load_npz(data / "CRC16_host_counts.npz")
    genes = pd.read_csv(data / "CRC16_host_features.csv").gene.tolist()
    metadata = pd.read_csv(data / "CRC16_metadata_raw.csv.gz", index_col="barcode")
    spots = pd.read_csv(data / "CRC16_frozen_niches.csv", index_col="barcode").loc[metadata.index]
    library = np.asarray(host.sum(axis=0)).ravel()
    detected = host.getnnz(axis=0)
    if host.shape != (36601, 3106) or not np.array_equal(library, metadata.nCount_Spatial):
        raise ValueError("Source counts/library integrity failed")
    niche = spots.Niche.to_numpy()
    xy = spots[["imagecol", "imagerow"]].to_numpy(float)
    scaled = (xy - xy.mean(0)) / xy.std(0, ddof=0)
    x, y = scaled.T
    common = np.column_stack([np.log1p(library), x, y, x*x, y*y, x*y])
    with threadpool_limits(limits=1):
        km = KMeans(n_clusters=12, init="k-means++", random_state=42, n_init=20,
                    max_iter=300, tol=0.0001, algorithm="lloyd").fit(xy)
    blocks = km.labels_ + 1
    blockmap = spots.reset_index()[["barcode", "imagecol", "imagerow", "Niche"]]
    blockmap["Block"] = blocks
    blockmap.to_csv(out / "CRC16_spatial_blocks.csv", index=False)
    centers = pd.DataFrame(km.cluster_centers_, columns=["imagecol", "imagerow"])
    centers.insert(0, "Block", np.arange(1, 13))
    centers["N_spots"] = np.bincount(km.labels_, minlength=12)
    centers.to_csv(out / "CRC16_block_centres.csv", index=False)
    (out / "CMH_partition_parameters.json").write_text(json.dumps({
        "n_clusters":12,"random_state":42,"n_init":20,"max_iter":300,"tol":0.0001,
        "algorithm":"lloyd","initialization":"k-means++","threads":1,
        "coordinates":["imagecol","imagerow"],"coordinate_scaling":False,
        "provenance":"Documented revised sensitivity partition; not a recovered historical CMH assignment",
        "continuity_correction":False,"zero_cell_pseudocount":False}, indent=2) + "\n")

    qc_reference = pd.read_csv(archived / "bacteria_qc_taxa.csv").rename(columns={"Unnamed: 0":"Taxon"})
    raw_taxa = metadata.columns[5:-1].tolist()
    if len(raw_taxa) != 316 or metadata.columns[-1] != "TotalMicrobe":
        raise ValueError("Unexpected raw taxon metadata schema")
    raw = metadata[raw_taxa].fillna(0).to_numpy(float)
    if not np.array_equal(raw, np.round(raw)) or (raw < 0).any():
        raise ValueError("Taxon counts are not nonnegative integers")
    filtered = metadata[qc_reference.Taxon].fillna(0).to_numpy(float)
    filtered[filtered == 1] = 0
    filtered = filtered.astype(np.int32)
    totals, prevalence = filtered.sum(0), (filtered > 0).sum(0)
    eligible = (totals >= 20) & (prevalence >= 10)
    qc = pd.DataFrame({"Taxon":qc_reference.Taxon,"total_umi":totals,
                       "prevalence_n":prevalence,"prevalence_fraction":prevalence/len(spots),"eligible":eligible})
    check_numeric(qc, qc_reference, ["Taxon"], ["total_umi","prevalence_n","prevalence_fraction"])
    qc.to_csv(out / "bacteria_qc_taxa.csv", index=False)
    pd.DataFrame({"Taxon":raw_taxa,"Retained":pd.Series(raw_taxa).isin(qc.Taxon)}).to_csv(out / "bacterial_qc_membership.csv", index=False)
    burden = filtered.sum(1)
    if burden.sum() != 80272 or (burden > 0).sum() != 720 or eligible.sum() != 20:
        raise ValueError("Bacterial QC totals differ from the source checkpoint")
    pd.DataFrame(filtered, index=metadata.index, columns=qc.Taxon).to_csv(out / "CRC16_bacterial_counts_postQC.csv.gz")
    spotout = spots.copy()
    spotout["host_library_size"] = library
    spotout["host_detected_features"] = detected
    spotout["bacterial_UMI"] = burden
    spotout["bacterial_occupancy"] = burden > 0
    spotout["CMH_block"] = blocks
    spotout.to_csv(out / "CRC16_spot_analysis_metadata.csv.gz")
    neighbors = cKDTree(xy).query(xy, k=7)[1][:, 1:]
    same = (niche[neighbors] == niche[:, None]).mean(1)

    gene_ref = pd.read_csv(ref / "gene_effects.csv")
    spatial43 = sorted(gene_ref.Gene.unique())
    counts, expr = expression(host, genes, library, spatial43)
    gene_rows, programme_rows = [], []
    for ni in NICHES:
        m = niche == ni
        for j, gene in enumerate(spatial43):
            gene_rows.append({"Gene":gene,"Niche":ni,"Set":"Convergent9" if gene in C9 else "Spatial43_only",
                "MeanLog2CPM_niche":expr[m,j].mean(),"MeanLog2CPM_rest":expr[~m,j].mean(),
                "Effect":expr[m,j].mean()-expr[~m,j].mean(),
                "Prevalence_niche":(counts[m,j]>0).mean(),"Prevalence_rest":(counts[~m,j]>0).mean()})
    gene_effects = pd.DataFrame(gene_rows)
    gene_error = check_numeric(gene_effects, gene_ref, ["Gene","Niche"],
                              ["MeanLog2CPM_niche","MeanLog2CPM_rest","Effect","Prevalence_niche","Prevalence_rest"])
    old_gene_reference = pd.read_csv(archived / "Spatial35_Convergent6_by_niche_final.csv")
    frozen_assignment_error = check_numeric(gene_effects[gene_effects.Gene.isin(old_gene_reference.Gene)], old_gene_reference,
        ["Gene","Niche"], ["MeanLog2CPM_niche","MeanLog2CPM_rest","Effect","Prevalence_niche","Prevalence_rest"])
    gene_effects.to_csv(out / "Spatial43_Convergent9_gene_effects.csv", index=False)
    immune_ref = pd.read_csv(ref / "immune_programmes.csv")
    definitions = immune_ref[["Axis","genes_available"]].drop_duplicates()
    definitions.to_csv(out / "immune_programme_gene_definitions.csv", index=False)
    for _, row in definitions.iterrows():
        wanted = row.genes_available.split(";")
        _, ex = expression(host, genes, library, wanted)
        score = zscores(ex).mean(1)
        for ni in NICHES:
            m = niche == ni
            programme_rows.append({"Axis":row.Axis,"Niche":ni,"n_genes":len(wanted),"genes_available":row.genes_available,
                "Mean_niche":score[m].mean(),"Mean_rest":score[~m].mean(),"Effect":score[m].mean()-score[~m].mean()})
    immune = pd.DataFrame(programme_rows)
    immune_error = check_numeric(immune, immune_ref, ["Axis","Niche"], ["Mean_niche","Mean_rest","Effect"])
    immune.to_csv(out / "immune_programmes_by_niche.csv", index=False)
    _, ex9 = expression(host, genes, library, C9)
    c9score = zscores(ex9).mean(1)
    c9rows, burdenrows, globalrows = [], [], []
    for ni in NICHES:
        m = niche == ni
        c9rows.append({"Niche":ni,"Convergent9_mean_score":c9score[m].mean(),
                      "Convergent9_score_effect":c9score[m].mean()-c9score[~m].mean()})
        burdenrows.append({"Niche":ni,"n_spots":int(m.sum()),"Bacteria_positive_fraction":(burden[m]>0).mean(),
            "Mean_filtered_bacterial_UMI":burden[m].mean(),"Median_filtered_bacterial_UMI":np.median(burden[m]),
            "Mean_log1p_filtered_bacterial_UMI":np.log1p(burden[m]).mean(),
            "Mean_host_UMI":library[m].mean(),"Median_host_UMI":np.median(library[m]),
            "Spatial_neighbor_enrichment":same[m].mean()/m.mean()})
        design = np.column_stack([np.ones(len(spots)),m.astype(float),common])
        beta, odds, probability, status = logistic(burden>0, design)
        b, p = ols(np.log1p(burden), design)
        globalrows.append({"Niche":ni,"n_spots":int(m.sum()),"AnyBacteria_niche":(burden[m]>0).mean(),
            "AnyBacteria_rest":(burden[~m]>0).mean(),"AnyBacteria_OR_adjusted":odds,"P_any":probability,
            "GeomeanBacterialBurden_niche":np.expm1(np.log1p(burden[m]).mean()),
            "GeomeanBacterialBurden_rest":np.expm1(np.log1p(burden[~m]).mean()),
            "LogBacterialBurden_beta_adjusted":b,"P_burden":p,"Occupancy_model_status":status})
    c9table, burden_table, global_table = pd.DataFrame(c9rows), pd.DataFrame(burdenrows), pd.DataFrame(globalrows)
    global_table["FDR_any"] = bh(global_table.P_any, 13)
    global_table["FDR_burden"] = bh(global_table.P_burden, 13)
    c9_error = check_numeric(c9table, pd.read_csv(ref / "niche_summary.csv"), ["Niche"],
                             ["Convergent9_mean_score","Convergent9_score_effect"])
    global_error = check_numeric(global_table, pd.read_csv(ref / "global_bacteria.csv"), ["Niche"],
        ["AnyBacteria_OR_adjusted","P_any","LogBacterialBurden_beta_adjusted","P_burden","FDR_any","FDR_burden"], tolerance=1e-6)
    c9table.to_csv(out / "Convergent9_programme_by_niche.csv", index=False)
    burden_table.to_csv(out / "bacterial_burden_by_niche.csv", index=False)
    global_table.to_csv(out / "global_bacterial_models.csv", index=False)
    global_publication = global_table.copy()
    gref = pd.read_csv(ref / "global_bacteria.csv").set_index("Niche").loc[global_publication.Niche]
    for col in ["AnyBacteria_OR_adjusted","P_any","FDR_any"]:
        global_publication[col] = gref[col].to_numpy()
    global_publication.to_csv(out / "global_bacterial_models_publication_precision.csv", index=False)

    rows, strata, hqrows = [], [], []
    high = (library >= 500) & (detected >= 200)
    hqeligible = (filtered[high].sum(0) >= 20) & ((filtered[high] > 0).sum(0) >= 10)
    for j, taxon in enumerate(qc.Taxon):
        if not eligible[j] and not hqeligible[j]:
            continue
        focal = filtered[:,j]
        positive = focal > 0
        nonfocal = np.log1p(burden - focal)
        for ni in NICHES:
            m = niche == ni
            design = np.column_stack([np.ones(len(spots)),m.astype(float),nonfocal,common])
            if eligible[j]:
                beta, odds, probability, status = logistic(positive, design)
                mo, mp, stat, orstatus, tables = cmh(m, positive, blocks)
                pin, pout = int(positive[m].sum()), int(positive[~m].sum())
                conditional, conditionalp = np.nan, np.nan
                if pin >= 3 and pout >= 3:
                    vin, vout = np.log2(focal[m & positive]+1), np.log2(focal[~m & positive]+1)
                    conditional = vin.mean()-vout.mean()
                    conditionalp = mannwhitneyu(vin, vout, alternative="two-sided").pvalue
                rows.append({"Taxon":taxon,"Niche":ni,"Niche_spots":int(m.sum()),"Positive_niche":pin,"Positive_rest":pout,
                    "Prevalence_niche":positive[m].mean(),"Prevalence_rest":positive[~m].mean(),
                    "Prevalence_fold":(positive[m].mean()+1e-9)/(positive[~m].mean()+1e-9),
                    "Mean_UMI_per_spot_niche":focal[m].mean(),"Mean_UMI_per_spot_rest":focal[~m].mean(),
                    "Logit_beta_adjusted":beta,"OR_adjusted":odds,"P_adjusted":probability,
                    "CMH_OR_spatial_blocks":mo,"P_spatial_CMH":mp,"CMH_statistic":stat,"CMH_OR_status":orstatus,
                    "Conditional_log2UMI_effect":conditional,"P_conditional_abundance":conditionalp,"Occupancy_model_status":status})
                strata.extend({"Taxon":taxon,"Niche":ni,**entry} for entry in tables)
            if hqeligible[j]:
                beta, odds, probability, status = logistic(positive[high], design[high])
                hqrows.append({"Taxon":taxon,"Niche":ni,"HQ_niche_spots":int((m&high).sum()),
                    "HQ_positive_niche":int((m&high&positive).sum()),"HQ_positive_rest":int((~m&high&positive).sum()),
                    "HQ_OR_adjusted":odds,"HQ_P":probability,"Occupancy_model_status":status})
    primary, hq = pd.DataFrame(rows), pd.DataFrame(hqrows)
    primary["FDR_adjusted"] = bh(primary.P_adjusted, 260)
    primary["FDR_spatial_CMH"] = bh(primary.P_spatial_CMH)
    primary["FDR_conditional_abundance"] = bh(primary.P_conditional_abundance)
    defined = np.isfinite(primary.Logit_beta_adjusted)&np.isfinite(primary.CMH_OR_spatial_blocks)&(primary.CMH_OR_spatial_blocks>0)
    concordant = np.sign(primary.Logit_beta_adjusted)==np.sign(np.log(primary.CMH_OR_spatial_blocks.clip(lower=1e-300)))
    primary["Direction_concordant"] = pd.Series(np.where(defined, concordant, None), dtype=object)
    primary["Robust_primary"] = defined&concordant&(primary.FDR_adjusted<0.05)&(primary.FDR_spatial_CMH<0.10)&(primary.Positive_niche>=5)&(primary.Positive_rest>=5)
    hq["HQ_FDR"] = bh(hq.HQ_P, 104)
    primary_reference = pd.read_csv(ref / "taxon_models.csv")
    primary_error = check_numeric(primary, primary_reference, ["Taxon","Niche"],
        ["Logit_beta_adjusted","OR_adjusted","P_adjusted","FDR_adjusted"], tolerance=5e-6)
    cmh_error = check_numeric(primary, primary_reference, ["Taxon","Niche"],
        ["CMH_OR_spatial_blocks","P_spatial_CMH","FDR_spatial_CMH","Conditional_log2UMI_effect","P_conditional_abundance","FDR_conditional_abundance"], tolerance=1e-8)
    hq_error = check_numeric(hq, pd.read_csv(ref / "high_complexity.csv"), ["Taxon","Niche"], ["HQ_OR_adjusted","HQ_P","HQ_FDR"], tolerance=5e-6)
    primary.to_csv(out / "taxon_niche_models_refitted.csv", index=False)
    hq.to_csv(out / "high_complexity_models_refitted.csv", index=False)
    hq_publication = hq.copy()
    href = pd.read_csv(ref / "high_complexity.csv").set_index(["Taxon","Niche"]).loc[pd.MultiIndex.from_frame(hq[["Taxon","Niche"]])]
    for col in ["HQ_OR_adjusted","HQ_P","HQ_FDR"]:
        hq_publication[col] = href[col].to_numpy()
    hq_publication.to_csv(out / "high_complexity_models_publication_precision.csv", index=False)
    pd.DataFrame(strata).to_csv(out / "CMH_stratum_contingency_tables.csv", index=False)
    # Preserve the already reviewed source precision only after independent refit checks.
    published = primary.copy()
    matched = primary_reference.set_index(["Taxon","Niche"]).loc[pd.MultiIndex.from_frame(primary[["Taxon","Niche"]])]
    for col in ["Logit_beta_adjusted","OR_adjusted","P_adjusted","FDR_adjusted"]:
        published[col] = matched[col].to_numpy()
    if not np.array_equal(published.Robust_primary.to_numpy(), matched.Robust_primary.to_numpy()):
        raise ValueError("Robust association identities disagree with the current reference")
    published.to_csv(out / "taxon_niche_models_publication_precision.csv", index=False)
    validation = {"host_genes":host.shape[0],"spots":len(spots),"niches":len(NICHES),"Spatial43":len(spatial43),"Convergent9":len(C9),
        "raw_taxa":len(raw_taxa),"raw_taxonomic_UMI":int(raw.sum()),"retained_bacterial_taxa":len(qc),"bacterial_UMI":int(burden.sum()),"occupied_spots":int((burden>0).sum()),
        "primary_planned":len(primary),"primary_estimable":int(primary.P_adjusted.notna().sum()),
        "CMH_finite":int(primary.P_spatial_CMH.notna().sum()),"CMH_unavailable":int(primary.P_spatial_CMH.isna().sum()),
        "CMH_stratum_tables":len(strata),"conditional_abundance_evaluable":int(primary.P_conditional_abundance.notna().sum()),
        "high_complexity_spots":int(high.sum()),"high_complexity_planned":len(hq),"high_complexity_estimable":int(hq.HQ_P.notna().sum()),
        "robust_primary":int(primary.Robust_primary.sum()),"CMH_OR_status_counts":primary.CMH_OR_status.value_counts().to_dict(),
        "gene_effect_max_errors":gene_error,"archived_35_gene_assignment_verification_max_errors":frozen_assignment_error,"immune_programme_max_errors":immune_error,
        "Convergent9_max_errors":c9_error,"global_model_max_errors":global_error,"primary_refit_max_errors":primary_error,
        "CMH_and_conditional_max_errors":cmh_error,"high_complexity_max_errors":hq_error}
    required={"raw_taxonomic_UMI":84044,"primary_planned":260,"primary_estimable":166,"CMH_finite":240,"CMH_unavailable":20,"CMH_stratum_tables":3120,
              "conditional_abundance_evaluable":119,"high_complexity_planned":104,"high_complexity_estimable":55}
    if any(validation[k] != v for k,v in required.items()):
        raise ValueError(f"Analysis denominators differ: {validation}")
    (out / "CRC16_validation.json").write_text(json.dumps(validation, indent=2)+"\n")
    return validation
