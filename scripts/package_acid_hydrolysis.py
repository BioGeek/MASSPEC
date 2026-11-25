"""Write the packaged copies of two partial acid hydrolysates of myoglobin.

    uv run python scripts/package_acid_hydrolysis.py

Both PRIDE projects hydrolysed horse myoglobin in hydrochloric acid in a
microwave oven (microwave-assisted acid hydrolysis, MAAH) and deposited
PEAKS's de novo reading of every MS/MS spectrum:

* PXD003804 (Savidor et al. 2017): 6 M HCl, 4 min; Q Exactive, MS1 scans
  from 300 m/z. File Myo_Rep1_de-Novo_Peaks-output.csv.
* PXD030203: MAAH beside proteinase K and papain; Orbitrap Fusion Lumos,
  MS1 scans from 350 m/z. File apo_MB_MAAH_1_peaks.csv.

The CSVs are downloaded into data/raw/ if they are not there yet. For each,
this keeps the de novo peptides read with an average local confidence (ALC)
of 50% or more, one row per distinct peptide (from its most confident
spectrum), and writes them to src/masspec/data/ with the precursor's
MEASURED neutral mass, (m/z - proton) * z. PEAKS's own "Mass" column is the
predicted peptide's mass, not the measured one, so it is not used.
"""
import csv
import urllib.request
from pathlib import Path

from masspec.masses import PROTON

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
OUT = ROOT / "src" / "masspec" / "data"
FTP = "https://ftp.pride.ebi.ac.uk/pride/data/archive"
RUNS = {
    "pxd003804_myo_rep1.tsv":
        f"{FTP}/2017/04/PXD003804/Myo_Rep1_de-Novo_Peaks-output.csv",
    "pxd030203_apo_mb_maah_1.tsv":
        f"{FTP}/2022/02/PXD030203/apo_MB_MAAH_1_peaks.csv",
}
MIN_ALC = 50.0


def fetch(url):
    """The local copy of the URL's file under data/raw/, downloaded once."""
    path = RAW / url.rsplit("/", 1)[1]
    if not path.exists():
        RAW.mkdir(parents=True, exist_ok=True)
        print(f"downloading {url}")
        urllib.request.urlretrieve(url, path)
    return path


def best_reads(path):
    """Each distinct peptide's most confident spectrum, at ALC >= 50."""
    best = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            alc = float(row["ALC (%)"])
            if alc < MIN_ALC:
                continue
            peptide = row["Peptide"]
            if peptide not in best or alc > best[peptide][0]:
                charge = int(row["z"])
                mass = (float(row["m/z"]) - PROTON) * charge
                best[peptide] = (alc, mass, charge)
    return best


for name, url in RUNS.items():
    reads = best_reads(fetch(url))
    with open(OUT / name, "w") as handle:
        handle.write("peptide\tmass\tcharge\talc\n")
        by_mass = sorted((mass, peptide, charge, alc)
                         for peptide, (alc, mass, charge) in reads.items())
        for mass, peptide, charge, alc in by_mass:
            handle.write(f"{peptide}\t{mass:.5f}\t{charge}\t{alc:g}\n")
    print(f"wrote {len(reads)} peptides to {(OUT / name).relative_to(ROOT)}")
