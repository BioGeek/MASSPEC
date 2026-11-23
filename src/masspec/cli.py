"""Command line.

    masspec reproduce
        the paper's insulin experiment, at 0.02 and 0.05 Da
    masspec simulate SEQ --out masses.txt
        present-day mass list for a protein (truth to masses.truth.json)
    masspec run masses.txt --ppm 5
        reconstruct from a mass list (or a centroided .mzML)
    masspec run masses.txt --truth SEQ
        ... and score the answer against a known sequence
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
    path = Path(arg)
    if path.exists():
        lines = path.read_text().splitlines()
        return normalise("".join(line.strip() for line in lines
                                 if line and not line.startswith(">")))
    return normalise(arg)


def cmd_reproduce(args):
    chain_a = "GJVEQCCASVCSJOQJENOCN"
    chain_b = "FVNQHJCGSHJVEAJOJVCGERGFFOTPKA"
    print("Bovine insulin, every fragment up to 9 residues below 1000 "
          "a.m.u., the paper's mass table.")
    print("Paper alphabet: J = Leu/Ile, O = Tyr. "
          "Chains come back in either direction.\n")
    for tol in (0.02, 0.05):
        observed = Observed(
            ideal_masses([chain_a, chain_b], 9, 1000.0, MASSES_1970,
                         H2O_1970),
            tol_da=tol)
        result = run(observed, MASSES_1970, H2O_1970)
        counts = {(stage, length): count
                  for stage, length, count in result.history}
        print(f"+-{tol} a.m.u.: compositions 1/2/3 long "
              f"{counts[('MASFRG', 1)]}/{counts[('MASFRG', 2)]}/"
              f"{counts[('MASFRG', 3)]} (paper 16/66/197 at 0.02), "
              f"after PRUNE {counts[('PRUNE', 3)]} "
              f"(93), after CQ {counts[('CQ', 3)]} (172), "
              f"10 long {counts[('OVERLAP', 10)]} (33)")
        for chain in result.chains(15):
            print("   ", chain, len(chain))
        print()


def cmd_simulate(args):
    from .simulate import coverage, simulate
    sim = simulate(_fasta_or_seq(args.sequence), enzyme=args.enzyme,
                   ppm=args.ppm_error, dropout=args.dropout,
                   n_noise=args.noise, mass_lo=args.lo, mass_hi=args.hi,
                   alkylated=not args.no_alkylation, met_ox=args.met_ox,
                   seed=args.seed)
    out = Path(args.out)
    out.write_text("\n".join(f"{mass:.6f}" for mass in sim.masses) + "\n")
    out.with_suffix(".truth.json").write_text(json.dumps(
        {"protein": sim.protein, "params": sim.params,
         "peptides": [(start, end, peptide, round(abundance, 6))
                      for start, end, peptide, abundance in sim.peptides]},
        indent=1))
    print(f"{len(sim.masses)} masses ({len(sim.peptides)} true fragments, "
          f"coverage of fragments up to "
          f"12 residues {coverage(sim, 12):.2f}) -> {out}")


def cmd_run(args):
    path = Path(args.masses)
    if path.suffix.lower() == ".mzml":
        from .features import masses_from_mzml
        feats = masses_from_mzml(path, ppm=args.feature_ppm,
                                 min_scans=args.min_scans)
        masses = [mass for mass, _, _, _ in feats]
        print(f"{len(masses)} features from {path.name}", file=sys.stderr)
    else:
        from .features import read_mass_list
        masses = read_mass_list(path, args.column)
    if args.calibrant:
        from .features import calibrate
        offset, peak, background = calibrate(
            masses, _fasta_or_seq(args.calibrant))
        print(f"calibrant: offset {offset:+.2f} ppm ({peak} matches in the "
              f"peak bin, {background:.0f} per background bin)",
              file=sys.stderr)
        masses = [mass / (1 + offset * 1e-6) for mass in masses]
    table = dict(MASSES)
    if not args.no_alkylation:
        table["C"] = MODIFIED["c"]
    for mod in args.mods or []:
        table[mod] = MODIFIED[mod]
    result = run_extend(FastObserved(masses, args.ppm), table, lo=args.lo,
                        hi=args.hi, beam=args.beam, support=args.support,
                        support_len=args.support_len)
    answer = assemble(result.contigs)
    for seq, support in answer:
        print(f"{support:6d}  {seq}")
    if args.truth:
        scores = score(answer, _fasta_or_seq(args.truth), table)
        print(json.dumps(
            {key: round(value, 3) if isinstance(value, float) else value
             for key, value in scores.items()},
            indent=1))


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="masspec", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser(
        "reproduce", help="the paper's insulin experiment",
    ).set_defaults(fn=cmd_reproduce)

    sim_parser = sub.add_parser(
        "simulate", help="present-day mass list for a protein")
    sim_parser.add_argument("sequence", help="sequence or FASTA file")
    sim_parser.add_argument("--out", default="masses.txt")
    sim_parser.add_argument(
        "--enzyme", default="pepsin",
        choices=["random", "pepsin", "proteinase_k", "acid"])
    sim_parser.add_argument("--ppm-error", type=float, default=3.0)
    sim_parser.add_argument("--dropout", type=float, default=0.2)
    sim_parser.add_argument("--noise", type=int, default=200)
    sim_parser.add_argument("--met-ox", type=float, default=0.3)
    sim_parser.add_argument("--lo", type=float, default=350.0)
    sim_parser.add_argument("--hi", type=float, default=6000.0)
    sim_parser.add_argument("--no-alkylation", action="store_true")
    sim_parser.add_argument("--seed", type=int, default=0)
    sim_parser.set_defaults(fn=cmd_simulate)

    run_parser = sub.add_parser(
        "run", help="reconstruct from a mass list or centroided mzML")
    run_parser.add_argument("masses")
    run_parser.add_argument("--column", help="mass column of a CSV/TSV")
    run_parser.add_argument("--ppm", type=float, default=10.0,
                            help="matching tolerance")
    run_parser.add_argument("--lo", type=float, default=350.0,
                            help="lowest mass the run can see")
    run_parser.add_argument("--hi", type=float, default=6000.0)
    run_parser.add_argument("--support", type=float, default=0.5)
    run_parser.add_argument("--support-len", type=int, default=20)
    run_parser.add_argument("--beam", type=int, default=2000)
    run_parser.add_argument("--no-alkylation", action="store_true",
                            help="Cys was not carbamidomethylated")
    run_parser.add_argument("--mods", nargs="*", choices=["m", "n", "q"],
                            help="variable mods as extra letters")
    run_parser.add_argument("--feature-ppm", type=float, default=10.0)
    run_parser.add_argument("--min-scans", type=int, default=2)
    run_parser.add_argument(
        "--calibrant",
        help="sequence or FASTA of a protein known to be in the run (pepsin "
        "in a pepsin digest), to remove the systematic mass error before "
        "matching")
    run_parser.add_argument(
        "--truth", help="known sequence or FASTA, to score the answer")
    run_parser.set_defaults(fn=cmd_run)

    args = parser.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
