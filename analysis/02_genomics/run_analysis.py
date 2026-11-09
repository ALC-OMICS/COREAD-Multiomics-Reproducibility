from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
for name in ["parse_exports.py","analyze_patients.py","pivot_external_table.py","make_patient_figures.py"]:
    subprocess.run([sys.executable,str(root/"code"/name)],cwd=root,check=True)
print("Numerical analysis and figures reproduced; see README for optional Excel export.")
