"""MASSPEC on a real present-day run.

The undeuterated pepsin digest of sperm whale myoglobin in PRIDE PXD071026
(an HDX-MS study, LTQ Orbitrap).

    # data: Mn_0s_1.raw and rcsb_pdb_6Z4T.fasta from the PRIDE project,
    # pepsin A from UniProt P00791
    mono ThermoRawFileParser.exe -i=data/raw/Mn_0s_1.raw -o=data/raw/ -f=1
    uv run python scripts/real_pxd071026.py

Prints: the features found, the systematic mass error (estimated from
pepsin's own autolysis peptides, never from myoglobin), how many features
are myoglobin fragments, and the reconstruction, scored against myoglobin
and against pepsin. Then the same reconstruction on ONLY the features that
match myoglobin: an oracle with every noise mass removed, to separate the
effect of noise from that of coverage.
"""
import bisect
from pathlib import Path

from masspec.evaluate import assemble, score
from masspec.features import calibrate, masses_from_mzml
from masspec.masses import H2O, MASSES
from masspec.modern import FastObserved, run_extend

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "raw"


def read_sequence(path):
    """The FASTA file's sequence, upper case, with I written as L."""
    residues = "".join(
        line.strip() for line in open(path) if not line.startswith(">")
    )
    return residues.upper().replace("I", "L")


def rounded(scores):
    """The score dictionary with every float rounded to two decimals."""
    return {
        key: round(value, 2) if isinstance(value, float) else value
        for key, value in scores.items()
    }


def fragments(protein, max_len=60):
    out = set()
    for start in range(len(protein)):
        fragment_mass = H2O
        for end in range(start, min(len(protein), start + max_len)):
            fragment_mass += MASSES[protein[end]]
            out.add(fragment_mass)
    return sorted(out)


def matching(masses, protein, ppm=2.0):
    fragment_masses = fragments(protein)
    keep = []
    for observed in masses:
        index = bisect.bisect_left(fragment_masses, observed)
        neighbours = fragment_masses[max(0, index - 1):index + 1]
        if any(
            abs(observed - candidate) / candidate * 1e6 <= ppm
            for candidate in neighbours
        ):
            keep.append(observed)
    return keep


myoglobin = read_sequence(DATA_DIR / "rcsb_pdb_6Z4T.fasta")
pepsin = read_sequence(DATA_DIR / "pepsin_A_mature.fasta")
features = masses_from_mzml(DATA_DIR / "Mn_0s_1.mzML", ppm=10, min_scans=2)
masses = [feature[0] for feature in features]
offset, peak, background = calibrate(masses, pepsin)
calibrated = [mass / (1 + offset * 1e-6) for mass in masses]
decoy = len(matching([mass + 0.3 for mass in calibrated], myoglobin))
print(
    f"{len(masses)} MS1 features; offset from pepsin autolysis peptides "
    f"{offset:+.2f} ppm ({peak} matches in the peak bin, "
    f"{background:.0f} per background bin)"
)
print(
    "features matching a myoglobin fragment within 2 ppm: "
    f"{len(matching(calibrated, myoglobin))} "
    f"(by chance, at +0.3 Da: {decoy}); "
    f"a pepsin fragment: {len(matching(calibrated, pepsin))}"
)
for tolerance in (3.0, 5.0):
    answer = assemble(
        run_extend(FastObserved(calibrated, tolerance), MASSES).contigs
    )
    print(
        f"reconstruction at {tolerance} ppm vs myoglobin: "
        f"{rounded(score(answer, myoglobin))}"
    )
    print(
        f"reconstruction at {tolerance} ppm vs pepsin:    "
        f"{rounded(score(answer, pepsin))}"
    )
only = matching(calibrated, myoglobin)
answer = assemble(run_extend(FastObserved(only, 3.0), MASSES).contigs)
print(
    f"oracle, the {len(only)} myoglobin-matching masses alone: "
    f"{rounded(score(answer, myoglobin))}"
)
