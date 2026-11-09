#!/usr/bin/env python3
"""Process one 8-µm Visium HD tumour from official 10x H5 and metadata.

The script extracts the frozen 466-gene catalogue, calculates log2(CPM + 1)
using library sizes from all 18,085 source features, aligns source spatial
annotations by barcode, calculates tumour/rest-of-tissue and 50-µm
periphery/interior effects, and computes global Moran's I on a binary
four-neighbour square-grid graph.

Raw H5 and metadata files are intentionally not distributed in the repository.
See data/README.md for official accessions and expected checksums.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

ALIASES = {
    "EPRS1": "EPRS",
    "PHB1": "PHB",
}


def decode_strings(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [value.decode("utf-8") if isinstance(value, (bytes, np.bytes_)) else str(value) for value in values],
        dtype=object,
    )


def sha256_file(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def parse_grid_coordinates(barcodes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    pattern = re.compile(r"s_008um_(\d+)_(\d+)-\d+$")
    rows = np.empty(len(barcodes), dtype=np.int32)
    cols = np.empty(len(barcodes), dtype=np.int32)
    for index, barcode in enumerate(barcodes):
        match = pattern.match(str(barcode))
        if match is None:
            raise ValueError(f"Unexpected Visium HD barcode format: {barcode}")
        rows[index] = int(match.group(1))
        cols[index] = int(match.group(2))
    return rows, cols


def build_rook_edges(rows: np.ndarray, cols: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return each undirected rook edge once, using right/down neighbours."""
    coordinate_to_index = {(int(row), int(col)): idx for idx, (row, col) in enumerate(zip(rows, cols))}
    left: list[int] = []
    right: list[int] = []
    for idx, (row, col) in enumerate(zip(rows, cols)):
        for neighbour in ((int(row) + 1, int(col)), (int(row), int(col) + 1)):
            other = coordinate_to_index.get(neighbour)
            if other is not None:
                left.append(idx)
                right.append(other)
    return np.asarray(left, dtype=np.int32), np.asarray(right, dtype=np.int32)


def extract_candidate_counts(
    h5_path: Path,
    genes: list[str],
    chunk_barcodes: int = 10_000,
) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray, pd.DataFrame, int]:
    """Return bins x genes raw counts and full-library sizes."""
    gene_to_output = {gene: idx for idx, gene in enumerate(genes)}

    with h5py.File(h5_path, "r") as handle:
        group = handle["matrix"]
        n_features, n_barcodes = (int(value) for value in group["shape"][:])
        feature_names = decode_strings(group["features"]["name"][:])
        feature_ids = decode_strings(group["features"]["id"][:])
        barcodes = decode_strings(group["barcodes"][:])

        source_to_output = np.full(n_features, -1, dtype=np.int32)
        mapping_rows: list[dict[str, object]] = []

        feature_lookup = {name: idx for idx, name in enumerate(feature_names)}
        for gene, output_idx in gene_to_output.items():
            source_name = gene if gene in feature_lookup else ALIASES.get(gene)
            if source_name is None or source_name not in feature_lookup:
                continue
            source_idx = feature_lookup[source_name]
            source_to_output[source_idx] = output_idx
            mapping_rows.append(
                {
                    "Gene": gene,
                    "source_feature_name": source_name,
                    "source_feature_id": feature_ids[source_idx],
                    "source_feature_index": source_idx,
                    "match_method": "exact_symbol" if source_name == gene else "historical_symbol",
                }
            )

        data_ds = group["data"]
        indices_ds = group["indices"]
        indptr_ds = group["indptr"]

        library_sizes = np.zeros(n_barcodes, dtype=np.float64)
        out_rows: list[np.ndarray] = []
        out_cols: list[np.ndarray] = []
        out_values: list[np.ndarray] = []

        for column_start in range(0, n_barcodes, chunk_barcodes):
            column_end = min(column_start + chunk_barcodes, n_barcodes)
            pointers = np.asarray(indptr_ds[column_start : column_end + 1], dtype=np.int64)
            data_start, data_end = int(pointers[0]), int(pointers[-1])
            source_rows = np.asarray(indices_ds[data_start:data_end], dtype=np.int32)
            values = np.asarray(data_ds[data_start:data_end])
            counts_per_column = np.diff(pointers)
            global_columns = np.repeat(
                np.arange(column_start, column_end, dtype=np.int32),
                counts_per_column,
            )
            local_columns = global_columns - column_start
            library_sizes[column_start:column_end] = np.bincount(
                local_columns,
                weights=values.astype(np.float64, copy=False),
                minlength=column_end - column_start,
            )

            mapped_rows = source_to_output[source_rows]
            keep = mapped_rows >= 0
            if np.any(keep):
                out_rows.append(mapped_rows[keep])
                out_cols.append(global_columns[keep])
                out_values.append(values[keep])

    if out_values:
        gene_by_bin = sparse.coo_matrix(
            (
                np.concatenate(out_values),
                (np.concatenate(out_rows), np.concatenate(out_cols)),
            ),
            shape=(len(genes), n_barcodes),
        ).tocsr()
        gene_by_bin.sum_duplicates()
    else:
        gene_by_bin = sparse.csr_matrix((len(genes), n_barcodes), dtype=np.int32)

    return (
        gene_by_bin.transpose().tocsr(),
        barcodes,
        library_sizes,
        pd.DataFrame(mapping_rows),
        n_features,
    )


def align_metadata(metadata_path: Path, barcodes: np.ndarray) -> pd.DataFrame:
    import gzip, io
    if metadata_path.name.endswith('.gz'):
        with gzip.open(metadata_path, 'rb') as f:
            metadata = pd.read_parquet(io.BytesIO(f.read()))
    else:
        metadata = pd.read_parquet(metadata_path)
    required = {"barcode", "tissue", "X", "Y", "Periphery"}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"Metadata lacks required columns: {sorted(missing)}")
    metadata = metadata.copy()
    metadata["barcode"] = metadata["barcode"].astype(str)
    if metadata['barcode'].duplicated().any():
        raise ValueError('Duplicate source metadata barcodes; ambiguous join rejected')
    metadata = metadata.set_index('barcode')
    aligned = metadata.reindex(pd.Index(barcodes.astype(str), name="barcode")).reset_index()
    if aligned["tissue"].isna().any():
        missing_count = int(aligned["tissue"].isna().sum())
        raise ValueError(f"{missing_count:,} matrix barcodes were absent from metadata")
    return aligned


def transform_counts(counts: sparse.csr_matrix, library_sizes: np.ndarray) -> sparse.csc_matrix:
    transformed = counts.tocoo(copy=True)
    denominators = library_sizes[transformed.row]
    positive = denominators > 0
    transformed.data = np.where(
        positive,
        np.log2(transformed.data.astype(np.float64) / denominators * 1e6 + 1.0),
        0.0,
    )
    transformed.eliminate_zeros()
    return transformed.tocsc()


def summarise_gene(
    vector: np.ndarray,
    region: np.ndarray,
    edge_left: np.ndarray,
    edge_right: np.ndarray,
) -> dict[str, float]:
    masks = {
        "all": np.ones(len(region), dtype=bool),
        "tissue": region == "Tissue",
        "interior": region == "Tumor",
        "periphery": region == "50 micron",
    }
    masks["tumor_combined"] = masks["interior"] | masks["periphery"]

    values: dict[str, float] = {}
    for name, mask in masks.items():
        n = int(mask.sum())
        if n == 0:
            values[f"mean_{name}"] = np.nan
            values[f"detect_{name}"] = np.nan
            continue
        subset = vector[mask]
        values[f"mean_{name}"] = float(subset.mean())
        values[f"detect_{name}"] = float(np.count_nonzero(subset) / n)

    values["tumor_minus_tissue"] = values["mean_tumor_combined"] - values["mean_tissue"]
    values["periphery_minus_interior"] = values["mean_periphery"] - values["mean_interior"]

    centered = vector - vector.mean()
    denominator = float(np.dot(centered, centered))
    if denominator == 0.0 or len(edge_left) == 0:
        moran = np.nan
    else:
        edge_sum = float(np.sum(centered[edge_left] * centered[edge_right]))
        moran = float(len(vector) / len(edge_left) * edge_sum / denominator)
    values["Moran_I"] = moran
    return values


def prepare_checkpoint(sample, h5_path, metadata_path, genes, output_dir):
    """Extract only candidate counts while preserving COMPLETE library sizes."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    counts, barcodes, libraries, mapping, n_features = extract_candidate_counts(Path(h5_path), genes)
    metadata = align_metadata(Path(metadata_path), barcodes)
    region = metadata['Periphery'].astype('string').fillna('NA').to_numpy(dtype=object)
    valid = np.isin(region, ['Tissue','Tumor','50 micron'])
    n_excluded = int((~valid).sum())
    if not valid.all():
        counts = counts[valid]
        barcodes = barcodes[valid]
        libraries = libraries[valid]
        region = region[valid]
    rows, cols = parse_grid_coordinates(barcodes)
    if len(set(zip(rows, cols))) != len(rows):
        raise ValueError('Duplicate grid coordinates')
    left, right = build_rook_edges(rows, cols)
    # Atomic compressed checkpoints include source barcode alignment and all-feature denominators.
    temp = output_dir/f'{sample}_candidate_counts.tmp.npz'
    final = output_dir/f'{sample}_candidate_counts.npz'
    sparse.save_npz(temp, counts, compressed=True)
    temp.replace(final)
    meta_path = output_dir/f'{sample}_bin_metadata.npz'
    np.savez_compressed(meta_path, library_sizes=libraries, grid_rows=rows, grid_cols=cols,
                        regions=np.array([{'Tissue':0,'Tumor':1,'50 micron':2}[x] for x in region], dtype=np.uint8),
                        barcodes=barcodes.astype(str), genes=np.asarray(genes, dtype=str))
    mapping.to_csv(output_dir/f'{sample}_feature_mapping.csv',index=False)
    qc = {'sample':sample,'source_h5_sha256':sha256_file(Path(h5_path)),
          'source_metadata_gzip_sha256':sha256_file(Path(metadata_path)),
          'n_features':int(n_features),'n_bins':len(barcodes),
          'region_counts':{x:int((region==x).sum()) for x in ['Tissue','Tumor','50 micron']},
          'zero_library_bins':int((libraries==0).sum()),'excluded_unannotated_bins':n_excluded,
          'rook_edges':len(left),'candidate_universe':len(genes),'evaluable_candidates':len(mapping),
          'count_checkpoint_sha256':sha256_file(final),'metadata_checkpoint_sha256':sha256_file(meta_path)}
    (output_dir/f'{sample}_extraction_QC.json').write_text(json.dumps(qc,indent=2)+'\n')
    return qc


def metrics_from_checkpoint(sample, checkpoint_dir):
    """Recompute all 466 regional and Moran endpoints from sparse source counts."""
    path = Path(checkpoint_dir)
    counts = sparse.load_npz(path/f'{sample}_candidate_counts.npz')
    metadata = np.load(path/f'{sample}_bin_metadata.npz', allow_pickle=False)
    mapping = pd.read_csv(path/f'{sample}_feature_mapping.csv').set_index('Gene')
    genes = metadata['genes'].tolist()
    libraries = metadata['library_sizes']
    region = np.asarray(['Tissue','Tumor','50 micron'])[metadata['regions']]
    assert counts.shape == (len(libraries), len(genes))
    left, right = build_rook_edges(metadata['grid_rows'], metadata['grid_cols'])
    matrix = transform_counts(counts, libraries)
    results = []
    for i, gene in enumerate(genes):
        record = {'Gene':gene, 'evaluable':gene in mapping.index,
                  'source_feature':mapping.loc[gene,'source_feature_name'] if gene in mapping.index else None}
        if gene in mapping.index:
            column=matrix.getcol(i)
            vector=np.zeros(matrix.shape[0],dtype=float)
            vector[column.indices]=column.data
            summary=summarise_gene(vector,region,left,right)
            record.update(summary)
            record['mean_log2CPM1_all']=record.pop('mean_all')
            record['detection_fraction_all']=record.pop('detect_all')
        results.append(record)
    qc = json.loads((path/f'{sample}_extraction_QC.json').read_text())
    qc['rook_edges_recomputed'] = len(left)
    qc['counts_integrity_passed'] = sha256_file(path/f'{sample}_candidate_counts.npz') == qc['count_checkpoint_sha256']
    qc['metadata_integrity_passed'] = sha256_file(path/f'{sample}_bin_metadata.npz') == qc['metadata_checkpoint_sha256']
    if not (qc['counts_integrity_passed'] and qc['metadata_integrity_passed']):
        raise ValueError('Checkpoint integrity check failed')
    return pd.DataFrame(results), qc

