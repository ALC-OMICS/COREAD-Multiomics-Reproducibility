#!/usr/bin/env python3
"""Join sequential HDF5 transfer parts and verify the expected SHA-256."""
from __future__ import annotations
import argparse, hashlib, re
from pathlib import Path

EXPECTED = "f435bb2651ff5297d0c24a99daf58850ed67ae1ed6c5ef05fad48fa3f0186670"
DEFAULT_NAME = "GSE178341_crc10x_full_c295v4_submit.h5"

def sha256(path: Path, block: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def part_number(path: Path) -> int:
    m = re.search(r"\.part(\d+)", path.name)
    if not m:
        raise ValueError(f"Cannot parse part number: {path.name}")
    return int(m.group(1))

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts-dir", type=Path, required=True)
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--expected-sha256", default=EXPECTED)
    args = ap.parse_args()
    out = args.output or (args.parts_dir / DEFAULT_NAME)
    parts = sorted(args.parts_dir.glob(f"{DEFAULT_NAME}.part*"), key=part_number)
    if not parts:
        raise SystemExit(f"No sequential parts found in {args.parts_dir}")
    nums = [part_number(p) for p in parts]
    if nums != list(range(1, len(parts) + 1)):
        raise SystemExit(f"Non-contiguous parts: {nums}")
    with out.open("wb") as dst:
        for p in parts:
            print(f"Appending {p.name} ({p.stat().st_size:,} bytes)")
            with p.open("rb") as src:
                while True:
                    b = src.read(8 * 1024 * 1024)
                    if not b:
                        break
                    dst.write(b)
    digest = sha256(out)
    print(f"SHA-256: {digest}")
    if args.expected_sha256 and digest.lower() != args.expected_sha256.lower():
        raise SystemExit("Checksum mismatch")
    print(f"Verified: {out}")

if __name__ == "__main__":
    main()
