"""Write and execute notebooks/masspec_1970.ipynb.

    uv run --group notebook python notebooks/build_notebook.py

The notebook's cells live here, as plain strings, so a change to the
walkthrough is a reviewable diff of this file rather than of notebook JSON.
Executing it stores every output in the .ipynb, so it reads on GitHub
without being run.
"""
from pathlib import Path

import nbformat
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent
cells = []


def add_markdown(text):
    """Append a markdown cell holding the text."""
    cells.append(nbformat.v4.new_markdown_cell(text.strip()))


def add_code(source):
    """Append a code cell holding the source."""
    cells.append(nbformat.v4.new_code_cell(source.strip()))


add_markdown(r"""
# MASSPEC (1970), step by step

Margaret O. Dayhoff and Richard V. Eck, "MASSPEC: A computer program for
complete sequence analysis of large proteins from mass spectrometry data of a
single sample", *Computers in Biology and Medicine* 1, 5-28 (1970),
[doi:10.1016/0010-4825(70)90013-2][doi].

[doi]: https://doi.org/10.1016/0010-4825(70)90013-2

Their proposal: break a protein near-randomly (partial acid hydrolysis), make
the fragments volatile, and measure the whole mixture in **one**
high-resolution mass spectrum. An exact mass fixes a fragment's elemental
formula but not its sequence, nor even its amino acid composition (Ala-Val and
Gly-Leu weigh the same). The sequence is recovered as the one set of
overlapping fragments that is mutually consistent among the thousands of
spurious candidates the same masses also allow.

They demonstrated it on simulated, ideal data for bovine insulin. This notebook
redoes that experiment with the reimplementation in this repository
(`src/masspec`), one program at a time, and regenerates each figure of the
paper next to the values the paper prints. **The code is a reimplementation
from the paper's description, not the original program, which appears to be
lost.**

The paper's one-letter code is used throughout: **J** is leucine or isoleucine
(they weigh the same, so the method cannot tell them apart) and **O** is
tyrosine.
""")

add_code(r"""
import random
from itertools import combinations_with_replacement

import matplotlib.pyplot as plt

from masspec.fragments import add_noise, fragen, ideal_masses
from masspec.masses import _FREE_1970, H2O, H2O_1970, MASSES, MASSES_1970
from masspec.mspec import (
    Observed,
    canon,
    cq,
    masfrg,
    overlap_step,
    prune,
    run,
)


def peptide_mass(residues):
    # peptide mass: residues plus one water
    return sum(MASSES_1970[residue] for residue in residues) + H2O_1970


def composition(residues):
    # composition, letters in alphabetical order
    return "".join(sorted(residues))


def difference(first, second):
    # first - second, rounded to the five decimals the paper prints. Adding
    # 0.0 turns a negative zero (a leftover such as -1e-12) into a plain 0.
    return round(first - second, 5) + 0.0


# insulin A chain, bovine (21 residues)
INSULIN_A = "GJVEQCCASVCSJOQJENOCN"
# insulin B chain, bovine (30 residues)
INSULIN_B = "FVNQHJCGSHJVEAJOJVCGERGFFOTPKA"
""")

add_markdown(r"""
## Fig. 2: the masses

The paper lists the monoisotopic masses of the free amino acids; a peptide
weighs their sum minus one water per peptide bond, which is the same as the sum
of RESIDUE masses plus one water. Below, the paper's values beside present-day
ones (free amino acid = residue + water).

One entry is wrong in the paper: tryptophan is printed as 188.09489, but free
tryptophan (C11H12N2O2) weighs 204.0899. The program evidently ran on the
printed value, because the counts in Figs. 6-8 are only reproduced with it (see
below). Insulin has no tryptophan, so the slip never showed.
""")

add_code(r"""
# Fig. 2's rows, in its order. It also lists Ile (I) and Leu (L), with the
# same mass as J; the program itself only uses J.
ACIDS = {
    "A": "Ala", "C": "Cys", "E": "Glu", "F": "Phe", "G": "Gly", "D": "Asp",
    "H": "His", "J": "Lux (Ile or Leu)", "K": "Lys", "M": "Met",
    "N": "Asn", "O": "Tyr", "P": "Pro", "Q": "Gln", "R": "Arg", "S": "Ser",
    "T": "Thr", "V": "Val", "W": "Trp",
}
modern = {"J": MASSES["L"], "O": MASSES["Y"], **MASSES}
print(
    f"{'Acid':<18}{'Code':<6}{'Mass':>11}{'Present-Day':>13}"
    f"{'Difference':>12}"
)
for code, acid in ACIDS.items():
    now = modern[code] + H2O  # free amino acid = residue + water
    printed = _FREE_1970[code]
    gap = difference(printed, now)
    print(f"{acid:<18}{code:<6}{printed:11.5f}{now:13.5f}{gap:12.5f}")
water_gap = difference(H2O_1970, H2O)
print(f"{'H2O':<24}{H2O_1970:11.5f}{H2O:13.5f}{water_gap:12.5f}")
""")

add_markdown(r"""
## Fig. 1: FRAGEN, fragmenting the protein

FRAGEN simulates a random hydrolysis by listing every fragment of the chain.
The masses will only say what a fragment is MADE of, so the paper prints the
amino acid composition of each, in alphabetical order. Fig. 1 shows the 24
fragments seven residues long from the B chain (the paper writes L for Leu
here, J elsewhere).
""")

add_code(r"""
def as_printed(residues):
    # a composition as Figs. 1, 3 and 4 print it: L for J, letters spaced
    return " ".join(residues.replace("J", "L"))


print("Structure of Initial Polymer, Bovine Insulin\n")
print("A chain - " + " ".join(INSULIN_A))
print("B chain - " + " ".join(INSULIN_B))
print("\nFragments from B Chain\n")
print(f"{'Fragments of Length 7':>21}   Amino Acid Composition")
seven = [INSULIN_B[start:start + 7] for start in range(len(INSULIN_B) - 6)]
for number, fragment in enumerate(seven, 1):
    print(f"{number:21}   {as_printed(composition(fragment))}")
""")

add_markdown(r"""
## Figs. 3 and 4: MASFRG, the masses of the fragments

MASFRG computes each fragment's mass. Fig. 3 lists the B chain's dipeptides,
Fig. 4 its octapeptides; the paper's values are copied below.
""")

add_code(r"""
paper_fig3 = {
    "AE": 218.09016, "AK": 217.14251, "AJ": 202.13162, "CG": 178.04113,
    "CJ": 234.10369, "CV": 220.08805, "EG": 204.07452, "ER": 303.15411,
    "EV": 246.12144, "FF": 312.14725, "FG": 222.10034, "FO": 328.14217,
    "FV": 264.14726, "GR": 231.13301, "GS": 162.06397, "HJ": 268.15340,
    "HQ": 283.12791, "HS": 242.10139, "KP": 243.15815, "JO": 294.15781,
    "JV": 230.16290, "NQ": 260.11193, "NV": 231.12178, "OT": 282.12144,
    "PT": 216.11089,
}


def paper_order(residues):
    # Fig. 3 sorts its rows as printed, with L for J: AE, AK, AL, CG...
    return residues.replace("J", "L")


dipeptides = sorted(
    {
        composition(INSULIN_B[start:start + 2])
        for start in range(len(INSULIN_B) - 1)
    },
    key=paper_order,
)
print("Peptides Two Amino Acids Long\n")
print(
    f"{'':3}{'Mass':>12}{'Paper':>12}{'Difference':>12}"
    "   Amino Acid Composition of Peptide"
)
for number, dipeptide in enumerate(dipeptides, 1):
    here = peptide_mass(dipeptide)
    printed = paper_fig3[dipeptide]
    gap = difference(here, printed)
    print(
        f"{number:3}{here:12.5f}{printed:12.5f}{gap:12.5f}"
        f"   {as_printed(dipeptide)}"
    )
""")

add_code(r"""
paper_fig4 = [
    916.42204, 856.38567, 894.37617, 893.41729, 864.42713, 856.41081,
    814.36389, 824.43874, 930.48057, 956.53258, 918.54208, 908.46723,
    866.42031, 866.42031, 951.48426, 895.42170, 879.42679, 913.41115,
    977.40606, 975.44454, 1015.47582, 1014.52817, 929.46422,
]
eight = [INSULIN_B[start:start + 8] for start in range(len(INSULIN_B) - 7)]
print("Peptides Eight Amino Acids Long\n")
print(
    f"{'':3}{'Mass':>12}{'Paper':>12}{'Difference':>12}"
    "   Amino Acid Composition of Peptide"
)
for number, (fragment, printed) in enumerate(zip(eight, paper_fig4), 1):
    here = peptide_mass(fragment)
    gap = difference(here, printed)
    print(
        f"{number:3}{here:12.5f}{printed:12.5f}{gap:12.5f}"
        f"   {as_printed(composition(fragment))}"
    )
""")

add_markdown(r"""
Every mass agrees with the paper to one unit in the fifth decimal: 22 of the
25 dipeptides exactly (FF, HQ and NQ are one unit lower in the paper), and
every octapeptide one unit lower in the paper.
Rounding rules alone cannot explain that: every mass in Fig. 2 has five
decimals, so these sums have no sixth digit to round or truncate. The original
program evidently carried more digits than it printed. Fig. 2 itself is
computed from elemental formulas with C = 12, H = 1.00782, N = 14.00307, O =
15.99491 and S = 31.97207 (all 18 values come out exactly), and the dipeptides
of Fig. 3 and the tetrapeptides of Fig. 5 all match if the water it subtracts
is 18.010565 rather than the printed 18.01056; the octapeptides fit no single
value. The rest is presumably the floating-point arithmetic of a 1970 machine,
which cannot be recovered from the printed numbers. The table also shows why
the method is hard: rows 13 and 14 are the same composition, and many different
peptides share a mass.

## Fig. 5: PEPMS, the table of all peptide masses

To go back from a mass to what could have made it, the program needs every
composition sorted by mass. PEPMS tabulates them (only up to four residues:
longer ones are generated later, from overlaps). Fig. 5 shows the stretch from
506.24 to 507.27. Each line is a group of compositions with identical formulas;
the step column is the gap to the previous mass. Unlike Figs. 1 and 4, which
are alphabetical, Fig. 5 writes each composition lightest residue first: VNQF
is V (99) < N (114) < Q (128) < F (147). So does this table.
""")

add_code(r"""
def lightest_first(residue):
    # Fig. 5's order within a composition: by residue mass, then letter.
    return (MASSES_1970[residue], residue)


table = {}
for length in range(1, 5):
    for residues in combinations_with_replacement(sorted(MASSES_1970), length):
        table.setdefault(round(peptide_mass(residues), 5), []).append(
            "".join(sorted(residues, key=lightest_first))
        )
window = sorted(value for value in table if 506.24 <= value <= 507.27)
print(
    f"{'Mass':>11}   {'Peptide Composition':<22}"
    "Increase in Mass from Preceding Entry"
)
previous = None
for group_mass in window:
    first, *others = table[group_mass]
    step = f"{group_mass - previous:.5f}" if previous else ""
    print(f"{group_mass:11.5f}   {first:<22}{step:>9}")
    for other in others:
        print(f"{'':11}   {other}")
    previous = group_mass
anfr_tjhh = difference(peptide_mass("ANFR"), peptide_mass("TJHH"))
jnnf_vnqf = difference(peptide_mass("JNNF"), peptide_mass("VNQF"))
print(f"\nANFR - TJHH = {anfr_tjhh:.5f}    JNNF - VNQF = {jnnf_vnqf:.5f}")
""")

add_markdown(r"""
The groups are the paper's, entry for entry: VNQF with JNNF at 506.24863, TQKM
with CKKE, then PPFF, PJMF and JJMM, and ANFR, GQFR and TJHH together at
506.25986, the paper's own example of a degeneracy (identical elemental
formulas). Within a group the paper's order follows no rule we could find
(probably the order its program generated them); here it is alphabetical. The
masses here are one unit higher in the fifth decimal, for the
reason given under Fig. 4. At a resolution of 0.02 a.m.u. many of these
neighbouring groups also blur together.

## The input: one spectrum of ideal masses

The insulin problem: every fragment of both chains up to nine residues long,
each mass observed once, exactly, and nothing above 1000 a.m.u.
""")

add_code(r"""
masses = ideal_masses(
    [INSULIN_A, INSULIN_B],
    max_len=9,
    mass_limit=1000.0,
    table=MASSES_1970,
    water=H2O_1970,
)
print(
    len(masses), "distinct masses, from", round(masses[0], 3),
    "to", round(masses[-1], 3),
)
observed = Observed(masses, tol_da=0.02)
""")

add_markdown(r"""
## Fig. 6: MASFRG, from masses back to compositions

Every composition up to four residues whose mass lies within ±0.02 of an
observed mass. "There are actually only 45 true fragments 3 long. However,
their masses correspond to any of 197 amino acid compositions."
""")

add_code(r"""
compositions = masfrg(observed, MASSES_1970, H2O_1970, max_len=4)
true_compositions = {
    length: {
        tuple(sorted(fragment))
        for fragment in fragen([INSULIN_A, INSULIN_B], 4)
        if len(fragment) == length
    }
    for length in range(1, 5)
}
print("Insulin of Beef 'A' and 'B'   Limits -.02 to +.02 a.m.u.\n")
for length in (1, 2, 3):
    listed = ["".join(residues) for residues in sorted(compositions[length])]
    print(f"There Are {len(listed)} Fragments of Length {length}")
    for row in range(0, len(listed), 10):
        chunk = listed[row:row + 10]
        print("   " + " ".join(f"{entry:<5}" for entry in chunk).rstrip())
""")

add_code(r"""
print(f"\n{'Length':>6}{'Fragments':>11}{'True Ones':>11}{'Paper':>7}")
for length, paper in zip(range(1, 5), (16, 66, 197, None)):
    print(
        f"{length:6}{len(compositions[length]):11}"
        f"{len(true_compositions[length]):11}{paper if paper else '':>7}"
    )
""")

add_markdown(r"""
## Fig. 7: PRUNE

"If the peptide (F,G,H,I) is true in any of its possible permutations, there
must also be at least two of the following: (F,G,H), (F,G,I), (F,H,I),
(G,H,I)." A true fragment's sub-fragments are all in the data; a spurious
composition's often are not.
""")

add_code(r"""
pruned = prune(compositions, min_subs=2)
for length in range(1, 5):
    print(
        f"length {length}: {len(compositions[length]):4} -> "
        f"{len(pruned[length]):4}"
    )
print(f"\ntripeptides after PRUNE: {len(pruned[3])} (paper: 93)")
""")

add_markdown(r"""
WEED, which removes compositions contained in no longer one, is not used: the
paper notes that the dipeptide list after the next step is "identical with the
previous list", i.e. that step ran on PRUNE's output. With WEED in between, the
next count comes out at 170 instead of 172.

## Fig. 8: CQ, from compositions to sequences

Each composition is expanded to every order of its residues, and a sequence of
length N is kept only if both sequences of length N-1 inside it (its first N-1
and its last N-1 residues) are kept. A mass carries no direction, so a sequence
and its reverse are one entry: "AEJ means both AEJ and JEA".
""")

add_code(r"""
sequences = cq(pruned)
print("Insulin of Beef 'A' and 'B'   Limits -.02 to +.02\n=CQ=")
print("Completely Sequenced Fragments Which May Be Present\n")
for length in range(1, 5):
    print(f"length {length}: {len(sequences[length])} sequences")
print(f"\nsequenced tripeptides: {len(sequences[3])} (paper: 172)")
true_tripeptides = {
    canon(fragment)
    for fragment in fragen([INSULIN_A, INSULIN_B], 3)
    if len(fragment) == 3
}
print(
    f"of which true: {len(sequences[3] & true_tripeptides)} "
    f"of {len(true_tripeptides)}"
)
""")

add_markdown(r"""
## Figs. 9 to 11: OVERLAP, growing the chain

From here on the candidates are built rather than looked up. Two sequences of
length N-1 that share N-2 residues overlap into one of length N, which is kept
only if its mass was observed. "Masses above a present limit are all assumed to
have been observed": above 1000 a.m.u. only consistency counts. A sequence that
overlaps no other cannot be extended, and is reported: "PIECES N LONG WHICH DO
NOT OVERLAP ANY OTHER PIECES".
""")

add_code(r"""
def accept(sequence):
    sequence_mass = peptide_mass(sequence)
    return sequence_mass > 1000.0 or observed.contains(sequence_mass)


sizes, finished = {}, {}
current, length = sequences[4], 5
while current:
    current, done = overlap_step(current, length, accept)
    sizes[length] = len(current)
    finished[length - 1] = sorted(done)
    length += 1

# The list sizes Figs. 9 to 11 print, by fragment length.
paper = {
    6: 241, 7: 92, 8: 49, 9: 37, 10: 33, 21: 11, 22: 9, 23: 8, 29: 2, 30: 1,
}
print("Insulin 'A' and 'B', Bovine")
print("Maximum Mass Measured = 1000.0\n")
for length in sorted(sizes):
    if sizes[length] == 0:
        continue
    noun = "Fragment" if sizes[length] == 1 else "Fragments"
    line = f"This Results in a List of {sizes[length]} {noun} of Size {length}"
    if length in paper:
        line += f"   (paper: {paper[length]})"
    print(line)
    pieces = finished.get(length, [])
    # The last list is the answer, printed below as the Final List.
    if pieces and 15 <= length < max(sizes) - 1:
        print(f"   Pieces {length} Long Which Do Not Overlap Any Other Pieces")
        for piece in pieces:
            print(f"      {piece}")
last = max(length for length in sizes if sizes[length])
print("\nFinal List")
for piece in finished[last]:
    print(f"   {piece}")
""")

add_markdown(r"""
The list grows, then shrinks to exactly the true fragments: 33 of length 10
("which is the correct number"), 11 of length 21. At length 21 one piece can no
longer be extended: the A chain, "written backwards, as it happens" in the
paper; here it comes out forwards, since only the direction of writing differs.
The program continues until one sequence 30 long remains: the B chain, here
written backwards.

From length 10 on, every list size the paper prints matches. The short lengths
do not: the paper keeps **241** hexapeptides (Fig. 9), 92 heptapeptides, 49
octapeptides and 37 nonapeptides (Fig. 10), against 136, 62, 40 and 35 here.
The gap closes as N grows, and the lists before (Figs. 6-8) and the answer
agree exactly, so the original evidently kept more candidates at those lengths,
by a rule the paper does not describe closely enough to reproduce. (Skipping
the mass check at length 6 alone gives 245 hexapeptides.)
""")

add_code(r"""
lengths = sorted(sizes)
figure, axes = plt.subplots(figsize=(8, 3.6))
axes.plot(
    lengths,
    [sizes[length] for length in lengths],
    color="#2563eb",
    lw=2,
    label="this reimplementation",
)
axes.scatter(
    list(paper),
    list(paper.values()),
    color="#dc2626",
    zorder=3,
    s=36,
    label="printed in the paper",
)
axes.axhline(0, color="#94a3b8", lw=0.8)
axes.set_xlabel("fragment length N")
axes.set_ylabel("candidate sequences kept")
axes.set_title("OVERLAP: candidates of each length, insulin at ±0.02 a.m.u.")
axes.legend(frameon=False)
plt.show()
""")

add_markdown(r"""
## The tryptophan misprint, demonstrated

Run the same lookup with the true tryptophan mass, and W (residue 186.079)
falls within 0.02 of the dipeptide G+E (186.064). Every "EG" in insulin then
also reads as a W.
""")

add_code(r"""
# Free tryptophan at its true mass, as a residue of the paper's table.
free_tryptophan = MASSES["W"] + H2O
masses_true_w = dict(MASSES_1970, W=free_tryptophan - H2O_1970)
compositions_true_w = masfrg(observed, masses_true_w, H2O_1970, 4)
print(
    "with the printed W:",
    [len(compositions[length]) for length in (1, 2, 3)],
    " (paper: 16, 66, 197)",
)
print(
    "with the true W:   ",
    [len(compositions_true_w[length]) for length in (1, 2, 3)],
)
print(
    "\nanswer with the true W:",
    run(observed, masses_true_w, H2O_1970).chains(25),
)
""")

add_markdown(r"""
## Resolution: ±0.02 against ±0.05

"If the maximum observable mass was limited to 1000 a.m.u., and the maximum
length of peptides to 9 links, the insulin problem could be solved with a
resolution of ±0.02 a.m.u., but could not be solved at ±0.05 a.m.u."
""")

add_code(r"""
for tolerance in (0.02, 0.05):
    result = run(Observed(masses, tol_da=tolerance), MASSES_1970, H2O_1970)
    chains = result.chains(20)
    print(f"±{tolerance}: {len(chains)} chains of 20 residues or more")
    for chain in chains:
        print("   ", chain)
""")

add_markdown(r"""
At ±0.05 glutamine and lysine (0.036 apart) can no longer be told apart, and
every Q/K position of the answer splits in two: four A chains, four B chains,
and some debris. The answer is not unique, so the problem is not solved.

## Practical considerations: noise


"Random noise masses can be tolerated in surprisingly large quantity because
the resulting peptides so seldom seem to overlap by chance ... If there were no
more than about 4000 actually observed masses in this range [100 to 1100],
including all noise, the spectrum would not be too crowded."
""")

add_code(r"""
for n_noise in (0, 500, 1000, 2000, 4000, 8000):
    noisy = add_noise(masses, n_random=n_noise, lo=100, hi=1100, seed=0)
    result = run(
        Observed(noisy, tol_da=0.02),
        MASSES_1970,
        H2O_1970,
        max_list=300_000,
    )
    chains = result.chains(15)
    both_found = INSULIN_A in chains and (
        INSULIN_B in chains or INSULIN_B[::-1] in chains
    )
    print(
        f"{n_noise:5} noise masses ({len(noisy):5} in all): "
        f"{len(chains)} chains of 15+, both insulin chains found: "
        f"{both_found}"
        + ("" if result.converged else f"  [{result.reason}]")
    )
""")

add_markdown(r"""
The reimplementation tolerates even more than the paper expected: both chains
still come back with 8000 random masses, twice the crowding the paper thought
the limit.

## Practical considerations: longer proteins

"We conclude that as an ideal problem, with true masses only, and no noise, a
protein of 100 links could be reconstructed, using the masses of fragments up
to only about 7 or 8 links long, having a maximum molecular mass in the order
of 800 a.m.u." Three random 100-residue proteins, exact masses, every fragment
below the mass limit:
""")

add_code(r"""
letters = sorted(MASSES_1970)
for seed in (1, 2, 3):
    rng = random.Random(seed)
    protein = "".join(rng.choice(letters) for _ in range(100))
    row = []
    for limit in (700, 800, 1000):
        spectrum = Observed(
            ideal_masses([protein], 60, limit, MASSES_1970, H2O_1970),
            tol_da=0.001,
        )
        result = run(
            spectrum,
            MASSES_1970,
            H2O_1970,
            mass_limit=limit,
            max_list=300_000,
        )
        found = any(
            chain in (protein, protein[::-1])
            for chain in result.chains(100)
        )
        row.append(f"{limit}: {'solved' if found else 'not solved'}")
    print(f"protein {seed}:  " + ",  ".join(row))
""")

add_markdown(r"""
The paper's estimate holds: with masses up to about 800-1000 a.m.u. a
100-residue protein comes back from its fragment masses alone, and below that
it may not.

## And today?

The paper ends: "The actual characteristics of a working system will determine
the practicability of the whole scheme." The rest of this repository asks that
question of present-day instruments: `simulate.py` makes the mass list of an
LC-MS run of a non-specific digest, `modern.py` adapts the reconstruction to
it, and `scripts/real_pxd071026.py` runs it on a real pepsin digest of
myoglobin. The README summarises what came out: the accuracy is there (a few
ppm), but a standard proteomics or HDX run does not contain the many short,
overlapping fragments the method lives on.
""")

notebook = nbformat.v4.new_notebook(cells=cells, metadata={
    "kernelspec": {
        "name": "python3",
        "display_name": "Python 3",
        "language": "python",
    },
})
NotebookClient(
    notebook,
    timeout=900,
    kernel_name="python3",
    resources={"metadata": {"path": str(HERE)}},
).execute()
out = HERE / "masspec_1970.ipynb"
nbformat.write(notebook, out)
print("wrote", out)
