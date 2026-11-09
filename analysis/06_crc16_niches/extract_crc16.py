#!/usr/bin/env python3
"""Extract the public CRC16 object and align the frozen CSV niche map.

This is a source-object extractor, not a replacement clustering algorithm.
Requires rdata 1.0.0 and its dependencies for the optional RDS extraction step.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
from scipy import sparse

SOURCE_SHA256 = "bb04b5f057f173f789fd9a36c2b30b17f31329b677e8812dfe19014aa43babc3"
SOURCE_URL = "https://raw.githubusercontent.com/ying-ge/FigureYa/main/FigureYa323STpathseq/STpathseq/Results/CRC_16/outs/ST.seu.rds"


def deref(obj):
    while obj is not None and obj.info.type.name == "REF":
        obj = obj.referenced_object
    return obj


def scalar(obj):
    obj = deref(obj)
    if obj is None:
        return None
    if obj.info.type.name == "SYM":
        return scalar(obj.value)
    if isinstance(obj.value, bytes):
        return obj.value.decode("utf-8")
    return obj.value


def attributes(obj):
    result = {}
    pair = deref(obj.attributes)
    while pair is not None and pair.info.type.name not in ("NILVALUE", "NIL"):
        result[scalar(pair.tag)] = deref(pair.value[0])
        pair = deref(pair.value[1])
    return result


def vector(obj):
    obj = deref(obj)
    if obj.info.type.name == "STR":
        return [scalar(item) for item in obj.value]
    values = np.asarray(obj.value)
    levels = attributes(obj).get("levels")
    if levels is not None:
        labels = vector(levels)
        return [labels[int(v) - 1] if v > 0 else None for v in values]
    if obj.info.type.name == "INT" and (values == np.iinfo(np.int32).min).any():
        values = values.astype(float)
        values[values == np.iinfo(np.int32).min] = np.nan
    return values


def dataframe(obj):
    attr = attributes(obj)
    return pd.DataFrame(
        {name: vector(item) for name, item in zip(vector(attr["names"]), obj.value)},
        index=pd.Index(vector(attr["row.names"]), name="barcode"),
    )


def recover_niches(frozen_path, coordinates, reference_sizes):
    """Align frozen host-only assignments by barcode; never recluster or infer labels."""
    frozen = pd.read_csv(frozen_path, index_col="barcode")
    if not coordinates.index.is_unique or not frozen.index.is_unique:
        raise ValueError("Duplicate source or frozen-assignment barcode")
    if set(frozen.index) != set(coordinates.index):
        raise ValueError("Frozen-assignment barcodes do not match source coordinates")
    frozen = frozen.loc[coordinates.index]
    columns = ["imagecol", "imagerow"]
    if not np.allclose(frozen[columns].to_numpy(float), coordinates[columns].to_numpy(float), rtol=0, atol=1e-6):
        raise ValueError("Frozen-assignment coordinates differ from source object")
    expected = reference_sizes.set_index("Niche")["n_spots"].astype(int).to_dict()
    if frozen["Niche"].value_counts().to_dict() != expected:
        raise ValueError("Frozen niche sizes differ from the recorded partition")
    out = coordinates.copy()
    for column in ["Niche", "svg_color", "svg_coordinate_error_pixels"]:
        out[column] = frozen[column]
    return out


def extract(rds_path, frozen_path, reference_sizes_path, out_dir):
    import rdata

    digest = hashlib.sha256(Path(rds_path).read_bytes()).hexdigest()
    if digest != SOURCE_SHA256:
        raise ValueError("Source RDS checksum mismatch; no extraction performed")
    source = rdata.parser.parse_file(rds_path).object
    slots = attributes(source)
    assays = dict(zip(vector(attributes(slots["assays"])["names"]), slots["assays"].value))
    count = attributes(attributes(assays["Spatial"])["counts"])
    matrix = sparse.csc_matrix((count["x"].value, count["i"].value, count["p"].value), shape=tuple(count["Dim"].value)).tocsr()
    matrix.data = matrix.data.astype(np.int32)
    names = count["Dimnames"].value
    genes, barcodes = vector(names[0]), vector(names[1])
    metadata = dataframe(slots["meta.data"]).loc[barcodes]
    coordinates = dataframe(attributes(slots["images"].value[0])["coordinates"]).loc[barcodes]
    if matrix.shape != (36601, 3106) or metadata.shape[1] != 322:
        raise ValueError("Unexpected source dimensions")
    libraries = np.asarray(matrix.sum(axis=0)).ravel()
    if not np.array_equal(libraries, metadata.nCount_Spatial):
        raise ValueError("All-feature library sizes disagree with source metadata")
    niche_map = recover_niches(frozen_path, coordinates, pd.read_csv(reference_sizes_path))
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(out_dir / "CRC16_host_counts.npz", matrix)
    pd.DataFrame({"gene": genes}).to_csv(out_dir / "CRC16_host_features.csv", index=False)
    metadata.to_csv(out_dir / "CRC16_metadata_raw.csv.gz")
    niche_map.to_csv(out_dir / "CRC16_frozen_niches.csv", index=True)
    record = {
        "source_url": SOURCE_URL, "source_sha256": digest,
        "source_features": matrix.shape[0], "source_spots": matrix.shape[1],
        "source_metadata_columns": metadata.shape[1],
        "frozen_niches": niche_map.Niche.nunique(),
        "maximum_svg_coordinate_error_pixels": float(niche_map.svg_coordinate_error_pixels.max()),
        "assignment_provenance": "Previously archived host-only CSV assignments aligned by unique source barcode and coordinates; no reclustering",
    }
    (out_dir / "CRC16_extraction_provenance.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rds", type=Path, required=True)
    parser.add_argument("--frozen-niches", type=Path, default=Path(__file__).resolve().parent/"input/CRC16_frozen_niches.csv")
    parser.add_argument("--reference-sizes", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    a = parser.parse_args()
    extract(a.rds, a.frozen_niches, a.reference_sizes, a.output_dir)
