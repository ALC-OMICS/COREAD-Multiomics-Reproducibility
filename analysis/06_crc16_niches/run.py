#!/usr/bin/env python3
"""Run current CRC16 analyses from included source measurement checkpoints."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import tempfile
os.environ.setdefault("MPLCONFIGDIR",str(Path(tempfile.gettempdir())/"coread_matplotlib"))
from threadpoolctl import threadpool_limits
from analyze_crc16 import run as analyze_crc16
from analyze_attribution import run as analyze_attribution
from plot_figure6 import draw


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode",choices=["processed"],default="processed")
    parser.add_argument("--output-dir",required=True,type=Path)
    args=parser.parse_args()
    base=Path(__file__).resolve().parent
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    status={"module":"06_crc16_niches","mode":args.mode,"passed":False}
    try:
        with threadpool_limits(limits=1):
            status["CRC16"]=analyze_crc16(base,out)
            status["GSE178341"]=analyze_attribution(base,out)
            status["figure"]=draw(base,out)
        status["input_sha256"]={str(p.relative_to(base)):hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in sorted((base/"input").glob("*")) if p.is_file()}
        status["passed"]=True
    except Exception as error:
        status["error"]=str(error)
        (out/"status.json").write_text(json.dumps(status,indent=2)+"\n")
        raise
    (out/"status.json").write_text(json.dumps(status,indent=2)+"\n")
    print(json.dumps({"module":status["module"],"passed":status["passed"],"output_dir":str(out)},indent=2))


if __name__=="__main__":
    main()
