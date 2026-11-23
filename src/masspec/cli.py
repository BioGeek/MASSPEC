"""Command line.

    masspec reproduce                      the paper's insulin experiment, at 0.02 and 0.05 Da
    masspec simulate SEQ --out masses.txt  present-day mass list for a protein (truth to masses.truth.json)
    masspec run masses.txt --ppm 5         reconstruct from a mass list (or a centroided .mzML)
    masspec run masses.txt --truth SEQ     ... and score the answer against a known sequence
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .evaluate import assemble, score
from .fragments import ideal_masses
from .masses import H2O_1970, MASSES, MASSES_1970, MODIFIED, normalise
from .modern import FastObserved, run_extend
from .mspec import Observed, run


def _fasta_or_seq(arg):
    p = Path(arg)
    if p.exists():
        lines = p.read_text().splitlines()
        return normalise("".join(x.strip() for x in lines if x and not x.startswith(">")))
    return normalise(arg)


def cmd_reproduce(a):
    A = "GJVEQCCASVCSJOQJENOCN"
    B = "FVNQHJCGSHJVEAJOJVCGERGFFOTPKA"
    print("Bovine insulin, every fragment up to 9 residues below 1000 a.m.u., the paper's mass table.")
    print("Paper alphabet: J = Leu/Ile, O = Tyr. Chains come back in either direction.\n")
    for tol in (0.02, 0.05):
        r = run(Observed(ideal_masses([A, B], 9, 1000.0, MASSES_1970, H2O_1970), tol_da=tol),
                MASSES_1970, H2O_1970)
        counts = {(s, n): c for s, n, c in r.history}
        print(f"+-{tol} a.m.u.: compositions 1/2/3 long {counts[('MASFRG', 1)]}/{counts[('MASFRG', 2)]}/"
              f"{counts[('MASFRG', 3)]} (paper 16/66/197 at 0.02), after PRUNE {counts[('PRUNE', 3)]} "
              f"(93), after CQ {counts[('CQ', 3)]} (172), 10 long {counts[('OVERLAP', 10)]} (33)")
        for c in r.chains(15):
            print("   ", c, len(c))
        print()


def cmd_simulate(a):
    from .simulate import coverage, simulate
    sim = simulate(_fasta_or_seq(a.sequence), enzyme=a.enzyme, ppm=a.ppm_error, dropout=a.dropout,
                   n_noise=a.noise, mass_lo=a.lo, mass_hi=a.hi, alkylated=not a.no_alkylation,
                   met_ox=a.met_ox, seed=a.seed)
    out = Path(a.out)
    out.write_text("\n".join(f"{m:.6f}" for m in sim.masses) + "\n")
    out.with_suffix(".truth.json").write_text(json.dumps(
        {"protein": sim.protein, "params": sim.params,
         "peptides": [(i, j, p, round(ab, 6)) for i, j, p, ab in sim.peptides]}, indent=1))
    print(f"{len(sim.masses)} masses ({len(sim.peptides)} true fragments, coverage of fragments up to "
          f"12 residues {coverage(sim, 12):.2f}) -> {out}")


def cmd_run(a):
    path = Path(a.masses)
    if path.suffix.lower() == ".mzml":
        from .features import masses_from_mzml
        feats = masses_from_mzml(path, ppm=a.feature_ppm, min_scans=a.min_scans)
        masses = [m for m, _, _, _ in feats]
        print(f"{len(masses)} features from {path.name}", file=sys.stderr)
    else:
        from .features import read_mass_list
        masses = read_mass_list(path, a.column)
    table = dict(MASSES)
    if not a.no_alkylation:
        table["C"] = MODIFIED["c"]
    for mod in a.mods or []:
        table[mod] = MODIFIED[mod]
    r = run_extend(FastObserved(masses, a.ppm), table, lo=a.lo, hi=a.hi, beam=a.beam,
                   support=a.support, support_len=a.support_len)
    answer = assemble(r.contigs)
    for seq, sc in answer:
        print(f"{sc:6d}  {seq}")
    if a.truth:
        s = score(answer, _fasta_or_seq(a.truth), table)
        print(json.dumps({k: round(v, 3) if isinstance(v, float) else v for k, v in s.items()}, indent=1))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="masspec", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("reproduce", help="the paper's insulin experiment").set_defaults(fn=cmd_reproduce)

    s = sub.add_parser("simulate", help="present-day mass list for a protein")
    s.add_argument("sequence", help="sequence or FASTA file")
    s.add_argument("--out", default="masses.txt")
    s.add_argument("--enzyme", default="pepsin", choices=["random", "pepsin", "proteinase_k", "acid"])
    s.add_argument("--ppm-error", type=float, default=3.0)
    s.add_argument("--dropout", type=float, default=0.2)
    s.add_argument("--noise", type=int, default=200)
    s.add_argument("--met-ox", type=float, default=0.3)
    s.add_argument("--lo", type=float, default=350.0)
    s.add_argument("--hi", type=float, default=6000.0)
    s.add_argument("--no-alkylation", action="store_true")
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(fn=cmd_simulate)

    r = sub.add_parser("run", help="reconstruct from a mass list or centroided mzML")
    r.add_argument("masses")
    r.add_argument("--column", help="mass column of a CSV/TSV")
    r.add_argument("--ppm", type=float, default=10.0, help="matching tolerance")
    r.add_argument("--lo", type=float, default=350.0, help="lowest mass the run can see")
    r.add_argument("--hi", type=float, default=6000.0)
    r.add_argument("--support", type=float, default=0.5)
    r.add_argument("--support-len", type=int, default=20)
    r.add_argument("--beam", type=int, default=2000)
    r.add_argument("--no-alkylation", action="store_true", help="Cys was not carbamidomethylated")
    r.add_argument("--mods", nargs="*", choices=["m", "n", "q"], help="variable mods as extra letters")
    r.add_argument("--feature-ppm", type=float, default=10.0)
    r.add_argument("--min-scans", type=int, default=2)
    r.add_argument("--truth", help="known sequence or FASTA, to score the answer")
    r.set_defaults(fn=cmd_run)

    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
