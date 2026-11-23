"""Where does MASSPEC work on present-day data? The paper's 'Practical considerations', redone.

    uv run python scripts/study.py > docs/study.md

Each row simulates one condition (simulate.py) for two proteins and three seeds, runs the modern
reconstruction (modern.run_extend), picks the answer (evaluate.assemble) and scores it against the
truth (evaluate.score). One factor varies per block; the rest stay at the baseline: pepsin, 3 ppm
mass error matched at 9 ppm, 20% dropout, 350 Da lower limit, 200 noise masses.
"""
import statistics
import time

from masspec.evaluate import assemble, score
from masspec.masses import MASSES
from masspec.modern import FastObserved, run_extend
from masspec.simulate import coverage, simulate

PROTEINS = {
    "ubiquitin (76)": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
    "sperm whale myoglobin (154)": "VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKVGVTALTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQGG",
}
BASE = dict(enzyme="pepsin", ppm=3.0, dropout=0.2, mass_lo=350.0, n_noise=200)
BLOCKS = [
    ("enzyme", [("random", {"enzyme": "random"}), ("pepsin", {}), ("proteinase K", {"enzyme": "proteinase_k"}), ("dilute acid", {"enzyme": "acid"})]),
    ("mass error", [("1 ppm", {"ppm": 1.0}), ("3 ppm", {}), ("10 ppm", {"ppm": 10.0}), ("30 ppm", {"ppm": 30.0})]),
    ("fragments missed", [("0%", {"dropout": 0.0}), ("20%", {}), ("40%", {"dropout": 0.4}), ("60%", {"dropout": 0.6})]),
    ("lowest mass seen", [("100 Da", {"mass_lo": 100.0}), ("350 Da", {}), ("600 Da", {"mass_lo": 600.0})]),
    ("noise masses", [("0", {"n_noise": 0}), ("200", {}), ("2000", {"n_noise": 2000})]),
]
SEEDS = (0, 1, 2)

print("| factor | value | protein | fragment coverage | recall | recall incl. mass-equivalent | precision | contigs | longest contig | seconds |")
print("|---|---|---|---|---|---|---|---|---|---|")
for factor, levels in BLOCKS:
    for label, over in levels:
        p = {**BASE, **over}
        for name, seq in PROTEINS.items():
            rows = []
            for seed in SEEDS:
                sim = simulate(seq, enzyme=p["enzyme"], ppm=p["ppm"], dropout=p["dropout"],
                               mass_lo=p["mass_lo"], n_noise=p["n_noise"], alkylated=False, seed=seed)
                t = time.time()
                r = run_extend(FastObserved(sim.masses, 3 * p["ppm"]), MASSES, lo=p["mass_lo"])
                ans = assemble(r.contigs)
                s = score(ans, sim.protein)
                rows.append((coverage(sim, 12), s["recall_exact"], s["recall_with_equivalents"],
                             s["precision_exact"], len(ans), max((len(a) for a, _ in ans), default=0),
                             time.time() - t))
            mean = lambda k: statistics.mean(r[k] for r in rows)
            print(f"| {factor} | {label} | {name} | {mean(0):.2f} | {mean(1):.2f} | {mean(2):.2f} | {mean(3):.2f} | "
                  f"{mean(4):.0f} | {mean(5):.0f} | {mean(6):.1f} |", flush=True)
