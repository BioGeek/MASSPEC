"""Write the packaged InstaNovo predictions that the notebook ends with.

The predictions come from a GPU, so they are made outside the notebook:

    # InstaNovo 1.2.2 (the current release) in its own environment; on a
    # Blackwell GPU (sm_120) torch needs the CUDA 12.8 build
    uv venv --python 3.12 .venv-instanovo
    uv pip install --python .venv-instanovo/bin/python \
        "instanovo==1.2.2" "torch==2.8.0" \
        --extra-index-url https://download.pytorch.org/whl/cu128 \
        --index-strategy unsafe-best-match
    .venv-instanovo/bin/instanovo predict -d data/raw/Mn_0s_1.mzML \
        -o Mn_0s_1_instanovo-v1.2.0.csv -i instanovo-v1.2.0

That is the release's default: the instanovo-v1.2.0 transformer, beam
search with 5 beams, and InstaNovo+ (instanovoplus-v1.1.0) refining every
prediction below 0.9 confidence. About 13 minutes on an RTX 5060 Ti.

    uv run python scripts/package_instanovo_pxd071026.py \
        Mn_0s_1_instanovo-v1.2.0.csv

keeps, for every predicted spectrum, its scan, precursor m/z and charge,
the final prediction (ProForma, UNIMOD modifications), its confidence (the
exponentiated log probability) and its precursor mass error in ppm, in
src/masspec/data/pxd071026_Mn_0s_1_instanovo.tsv.
"""
import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "src" / "masspec" / "data" / "pxd071026_Mn_0s_1_instanovo.tsv"

with open(sys.argv[1], newline="") as handle:
    rows = list(csv.DictReader(handle))
with open(OUT, "w") as handle:
    handle.write("scan\tprecursor_mz\tcharge\tprediction\tconfidence\t"
                 "ppm_error\n")
    for row in rows:
        confidence = math.exp(float(row["log_probs"]))
        ppm_error = (f"{float(row['delta_mass_ppm']):.2f}"
                     if row["delta_mass_ppm"] else "")
        handle.write(f"{row['scan_number']}\t{float(row['precursor_mz']):.5f}"
                     f"\t{row['precursor_charge']}\t{row['predictions']}\t"
                     f"{confidence:.4f}\t{ppm_error}\n")
print(f"wrote {len(rows)} predictions to {OUT.relative_to(ROOT)}")
