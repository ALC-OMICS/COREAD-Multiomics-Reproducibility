#!/usr/bin/env python3
"""Compatibility entrypoint for the canonical CRC16 source extractor."""
from pathlib import Path
import runpy
if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).resolve().parents[1]/"06_crc16_niches/extract_crc16.py"), run_name="__main__")
