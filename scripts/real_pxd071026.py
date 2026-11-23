"""MASSPEC on a real present-day run: the undeuterated pepsin digest of sperm whale myoglobin in
PRIDE PXD071026 (an HDX-MS study, LTQ Orbitrap).

    # data: Mn_0s_1.raw and rcsb_pdb_6Z4T.fasta from the PRIDE project, pepsin A from UniProt P00791
    mono ThermoRawFileParser.exe -i=data/raw/Mn_0s_1.raw -o=data/raw/ -f=1
    uv run python scripts/real_pxd071026.py

Prints: the features found, the systematic mass error (estimated from pepsin's own autolysis
peptides, never from myoglobin), how many features are myoglobin fragments, and the reconstruction,
scored against myoglobin and against pepsin. Then the same reconstruction on ONLY the features that
match myoglobin: an oracle with every noise mass removed, to separate the effect of noise from that
of coverage.
"""
import bisect
from pathlib import Path

from masspec.evaluate import assemble, score
from masspec.features import calibrate, masses_from_mzml
from masspec.masses import H2O, MASSES
from masspec.modern import FastObserved, run_extend

D = Path(__file__).resolve().parent.parent / "data" / "raw"
seq = lambda p: "".join(x.strip() for x in open(p) if not x.startswith(">")).upper().replace("I", "L")
myo, pepsin = seq(D / "rcsb_pdb_6Z4T.fasta"), seq(D / "pepsin_A_mature.fasta")
fmt = lambda s: {k: round(v, 2) if isinstance(v, float) else v for k, v in s.items()}


def fragments(protein, max_len=60):
    out = set()
    for i in range(len(protein)):
        m = H2O
        for j in range(i, min(len(protein), i + max_len)):
            m += MASSES[protein[j]]
            out.add(m)
    return sorted(out)


def matching(masses, protein, ppm=2.0):
    fr = fragments(protein)
    keep = []
    for x in masses:
        k = bisect.bisect_left(fr, x)
        if any(abs(x - c) / c * 1e6 <= ppm for c in fr[max(0, k - 1):k + 1]):
            keep.append(x)
    return keep


feats = masses_from_mzml(D / "Mn_0s_1.mzML", ppm=10, min_scans=2)
masses = [f[0] for f in feats]
off, peak, bg = calibrate(masses, pepsin)
cal = [m / (1 + off * 1e-6) for m in masses]
decoy = len(matching([m + 0.3 for m in cal], myo))
print(f"{len(masses)} MS1 features; offset from pepsin autolysis peptides {off:+.2f} ppm "
      f"({peak} matches in the peak bin, {bg:.0f} per background bin)")
print(f"features matching a myoglobin fragment within 2 ppm: {len(matching(cal, myo))} "
      f"(by chance, at +0.3 Da: {decoy}); a pepsin fragment: {len(matching(cal, pepsin))}")
for tol in (3.0, 5.0):
    ans = assemble(run_extend(FastObserved(cal, tol), MASSES).contigs)
    print(f"reconstruction at {tol} ppm vs myoglobin: {fmt(score(ans, myo))}")
    print(f"reconstruction at {tol} ppm vs pepsin:    {fmt(score(ans, pepsin))}")
only = matching(cal, myo)
ans = assemble(run_extend(FastObserved(only, 3.0), MASSES).contigs)
print(f"oracle, the {len(only)} myoglobin-matching masses alone: {fmt(score(ans, myo))}")
