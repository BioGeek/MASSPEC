"""Where does MASSPEC work on present-day data?

The paper's 'Practical considerations', redone.

    uv run python scripts/study.py > docs/study.md

Each row simulates one condition (simulate.py) for two proteins and three
seeds, runs the modern reconstruction (modern.run_extend), picks the answer
(evaluate.assemble) and scores it against the truth (evaluate.score). One
factor varies per block; the rest stay at the baseline: pepsin, 3 ppm mass
error matched at 9 ppm, 20% dropout, 350 Da lower limit, 200 noise masses.
"""
import statistics
import time

from masspec.evaluate import assemble, score
from masspec.masses import MASSES
from masspec.modern import FastObserved, run_extend
from masspec.simulate import coverage, simulate

PROTEINS = {
    "ubiquitin (76)": (
        "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKE"
        "STLHLVLRLRGG"
    ),
    "sperm whale myoglobin (154)": (
        "VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKV"
        "GVTALTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQ"
        "GAMNKALELFRKDIAAKYKELGYQGG"
    ),
}
BASE = dict(
    enzyme="pepsin", ppm=3.0, dropout=0.2, mass_lo=350.0, n_noise=200
)
BLOCKS = [
    ("enzyme", [
        ("random", {"enzyme": "random"}),
        ("pepsin", {}),
        ("proteinase K", {"enzyme": "proteinase_k"}),
        ("dilute acid", {"enzyme": "acid"}),
    ]),
    ("mass error", [
        ("1 ppm", {"ppm": 1.0}),
        ("3 ppm", {}),
        ("10 ppm", {"ppm": 10.0}),
        ("30 ppm", {"ppm": 30.0}),
    ]),
    ("fragments missed", [
        ("0%", {"dropout": 0.0}),
        ("20%", {}),
        ("40%", {"dropout": 0.4}),
        ("60%", {"dropout": 0.6}),
    ]),
    ("lowest mass seen", [
        ("100 Da", {"mass_lo": 100.0}),
        ("350 Da", {}),
        ("600 Da", {"mass_lo": 600.0}),
    ]),
    ("noise masses", [
        ("0", {"n_noise": 0}),
        ("200", {}),
        ("2000", {"n_noise": 2000}),
    ]),
]
SEEDS = (0, 1, 2)


def column_mean(rows, column):
    """The mean of one column over the rows."""
    return statistics.mean(row[column] for row in rows)


print(
    "| factor | value | protein | fragment coverage | recall "
    "| recall incl. mass-equivalent | precision | contigs "
    "| longest contig | seconds |"
)
print("|---|---|---|---|---|---|---|---|---|---|")
for factor, levels in BLOCKS:
    for label, overrides in levels:
        params = {**BASE, **overrides}
        for name, sequence in PROTEINS.items():
            rows = []
            for seed in SEEDS:
                sim = simulate(
                    sequence,
                    enzyme=params["enzyme"],
                    ppm=params["ppm"],
                    dropout=params["dropout"],
                    mass_lo=params["mass_lo"],
                    n_noise=params["n_noise"],
                    alkylated=False,
                    seed=seed,
                )
                started = time.time()
                result = run_extend(
                    FastObserved(sim.masses, 3 * params["ppm"]),
                    MASSES,
                    lo=params["mass_lo"],
                )
                answer = assemble(result.contigs)
                scores = score(answer, sim.protein)
                longest = max(
                    (len(contig) for contig, _ in answer), default=0
                )
                rows.append((
                    coverage(sim, 12),
                    scores["recall_exact"],
                    scores["recall_with_equivalents"],
                    scores["precision_exact"],
                    len(answer),
                    longest,
                    time.time() - started,
                ))
            means = [column_mean(rows, column) for column in range(7)]
            print(
                f"| {factor} | {label} | {name} | {means[0]:.2f} "
                f"| {means[1]:.2f} | {means[2]:.2f} | {means[3]:.2f} | "
                f"{means[4]:.0f} | {means[5]:.0f} | {means[6]:.1f} |",
                flush=True,
            )
