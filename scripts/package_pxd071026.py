"""Write the packaged copy of the real run that the notebook uses.

    uv run --extra real python scripts/package_pxd071026.py

Reads data/raw/Mn_0s_1.mzML (see scripts/real_pxd071026.py for how to get
it) and writes src/masspec/data/pxd071026_Mn_0s_1.tsv: the run's MS1 scans
merged into one list of neutral masses by features.masses_from_mzml(), the
same call real_pxd071026.py makes, at full float precision so the packaged
list gives the same numbers as the mzML. Uncalibrated: calibration is a step
the notebook shows. The two sequences beside it are copied in by hand once.
"""
from pathlib import Path

from masspec.features import masses_from_mzml

ROOT = Path(__file__).resolve().parent.parent
RUN = ROOT / "data" / "raw" / "Mn_0s_1.mzML"
OUT = ROOT / "src" / "masspec" / "data" / "pxd071026_Mn_0s_1.tsv"

features = masses_from_mzml(RUN, ppm=10, min_scans=2)
with open(OUT, "w") as handle:
    handle.write("mass\tintensity\tscans\tcharges\n")
    for mass, intensity, scans, charges in features:
        charge_list = ",".join(str(charge) for charge in charges)
        handle.write(f"{float(mass)!r}\t{float(intensity)!r}\t{scans}\t"
                     f"{charge_list}\n")
print(f"wrote {len(features)} masses to {OUT.relative_to(ROOT)}")
