# MASSPEC, reverse-engineered

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/BioGeek/MASSPEC/blob/main/notebooks/masspec_1970.ipynb)

A reimplementation of **MASSPEC**, [Margaret O. Dayhoff](https://en.wikipedia.org/wiki/Margaret_Oakley_Dayhoff)
and [Richard V. Eck](https://pmc.ncbi.nlm.nih.gov/articles/PMC2206694/)'s 1970 program for
sequencing a whole protein from one high-resolution mass spectrum of its partial hydrolysate
(*Computers in Biology and Medicine* 1, 5-28, 1970,
[doi:10.1016/0010-4825(70)90013-2](https://doi.org/10.1016/0010-4825(70)90013-2)), and an extension
that runs it on present-day data.

The [notebook](notebooks/masspec_1970.ipynb) walks through the paper program by program and
rebuilds every figure, then runs the method on present-day data.

**This is not Dayhoff and Eck's code**, which appears to be lost. It was written from the paper's
description and its appendix, which describes each program step by step. Where the paper is
ambiguous, the choice made here is stated below.

```bash
uv sync                        # core: stdlib only
uv sync --extra real           # + pyteomics, lxml, numpy, for reading mzML
uv run masspec reproduce       # the paper's insulin experiment
uv run masspec simulate SEQ_OR_FASTA --out masses.txt
uv run masspec run masses.txt --ppm 9 --truth SEQ_OR_FASTA
uv run masspec run run.mzML --ppm 9 --truth SEQ_OR_FASTA
uv run pytest                  # includes the style check below
uv run ruff check .            # PEP 8, 79 columns
```

## The 1970 method (`mspec.py`)

The idea: digest a protein near-randomly, make the fragments volatile, and measure the whole mixture
in one high-resolution spectrum. Each exact mass gives a few possible COMPOSITIONS, the sequence is
the one set of overlapping fragments that stays mutually consistent among thousands of spurious
candidates. "It is somewhat like a 1000-piece jig-saw puzzle, of which all the pieces are known to be
present, but mixed with several thousand other pieces taken at random from other puzzles."

| Program | What it does |
|---|---|
| `masfrg` | every composition up to 4 residues whose mass matches an observed mass |
| `prune` | keep a composition only if 2 of its one-shorter sub-compositions are kept |
| `weed` | keep a composition only if a one-longer kept composition contains it |
| `cq` | compositions to sequences whose prefix and suffix are both kept |
| `overlap_step` | two length N-1 sequences sharing N-2 residues make a length N candidate, kept if its mass was observed (masses above the limit are assumed observed) |

A sequence and its reverse are one entry: masses carry no direction.

### Reproducing the paper

`masspec reproduce` runs the insulin experiment: every fragment of bovine insulin's two chains up to
9 residues long, masses below 1000, ideal (exact) data, the paper's own mass table.

| | paper | here |
|---|---|---|
| compositions 1 / 2 / 3 residues long, at +-0.02 | 16 / 66 / 197 (Fig. 6) | 16 / 66 / 197 |
| tripeptide compositions after PRUNE | 93 (Fig. 7) | 93 |
| tripeptide sequences after CQ | 172 (Fig. 8) | 172 |
| sequences of length 6 / 7 / 8 / 9 | 241 / 92 / 49 / 37 (Figs. 9, 10) | 136 / 62 / 40 / 35; with the hexapeptides unchecked, 245 / 96 / 50 / 37 |
| 10-residue sequences | 33, "the correct number" | 33 |
| answer at +-0.02 | both chains | both chains, A forwards, B backwards |
| answer at +-0.05 | "could not be solved" | 4 A chains and 4 B chains, differing in Q/K |

Three things had to be decided:

* **The tryptophan mass.** Fig. 2 prints free tryptophan as 188.09489; it is 204.0899. With the true
  value W sits within 0.02 of G+E, the counts become 17 / 71 / 208, and a second B chain with "W" for
  "EG" survives. With the printed value everything matches, so the program evidently ran on it.
  Insulin has no tryptophan, which is why nobody noticed.
* **WEED is off by default.** The paper says the dipeptide list after CQ is "identical with the
  previous list", i.e. CQ ran on PRUNE's output; with WEED the counts are 170 and 55 instead of 172
  and 56. `run(use_weed=True)` turns it on.
* **The short OVERLAP lists.** Followed to the letter, the text (each candidate
  kept only if its mass was observed) keeps fewer candidates of length 6 to 9
  than the paper prints. The printed figures show what the original did
  instead: Fig. 9, headed INPUT FRAGMENTS, lists hexapeptides whose masses were
  never observed (its first row holds ACGEVA, 548.226), and the original built
  its heptapeptides from all of them. With the hexapeptide mass check skipped,
  this code lists 245 input hexapeptides whose first forty are the paper's
  first forty, and keeps 96, 50 and 37 of lengths 7, 8 and 9 against the
  paper's 92, 49 and 37. The default follows the text; the notebook shows both.

`notebooks/masspec_1970.ipynb` walks through all of this step by step, regenerating each figure of
the paper (Figs. 1 to 11) next to the values the paper prints, and redoes its "practical
considerations": noise, resolution, and 100-residue proteins. It is built and executed by
`notebooks/build_notebook.py` (`uv sync --group notebook` first). On Colab it installs the
package itself: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/BioGeek/MASSPEC/blob/main/notebooks/masspec_1970.ipynb)

## On present-day data (`modern.py`, `simulate.py`, `features.py`)

Dayhoff and Eck's data never existed: their simulation assumed every fragment observed, exactly. An
LC-MS run of a non-specific digest (pepsin, as in every HDX-MS experiment, or proteinase K) breaks
that in three ways, and the method was changed in three ways to match:

1. **Short fragments are not seen.** Below roughly 350 Da peptides hardly retain on a C18 column. The
   paper already assumed masses ABOVE its instrument's limit observed; here masses outside an
   observable window `[lo, hi]` are, at both ends.
2. **Many fragments are missed.** The 1970 OVERLAP builds a length-N candidate only from two
   surviving length N-1 windows, so one missing fragment breaks every chain through it. On simulated
   pepsin data that recovered 17 of ubiquitin's 76 residues. The modern mode keeps MASSPEC's seeds
   (MASFRG, PRUNE and CQ with the window) but grows each chain one residue at a time, judging an
   extension by the fragments that END at the new residue: at least one, and half of those in the
   window, must be observed. A mass may support a chain only once: it is one fragment of the
   protein. (Without that rule, low-complexity chains such as ALASALSAL... grew forever on the same
   few masses.)
3. **Masses are measured in ppm**, and Cys is carbamidomethylated (on by default), with optional
   variable modifications as extra residue letters (`--mods m n q`: oxidised Met, deamidated N/Q).

The answer is a set of contigs. Masses alone cannot order the contigs, cannot tell I from L, and
cannot tell Q from GA or N from GG (identical formulas); `evaluate.score` counts such mass-identical
swaps apart from real errors.

### Where it works (`scripts/study.py`, `docs/study_table.md`)

Simulated pepsin digests of ubiquitin (76 residues) and sperm whale myoglobin (154), three seeds
each. Baseline: 3 ppm mass error, 20% of detectable fragments missed, nothing seen below 350 Da,
200 noise masses. One factor varied at a time. Recall: the fraction of the protein's residues
covered correctly. Precision: the fraction of reported residues that are correct.

| condition | ubiquitin recall / precision | myoglobin recall / precision |
|---|---|---|
| baseline (pepsin, 3 ppm, 20% missed) | 0.71 / 0.46 | 0.86 / 0.85 |
| uniform cleavage (the paper's assumption) | 0.91 / 0.74 | 0.98 / 0.69 |
| dilute acid | 0.88 / 0.95 | 0.97 / 0.97 |
| 1 ppm | 0.71 / 0.58 | 0.81 / 0.87 |
| 10 ppm | 0.62 / 0.09 | 0.42 / 0.05 |
| 30 ppm | 0 / 0 | 0 / 0 |
| no fragments missed | 0.86 / 0.70 | 0.91 / 0.94 |
| 40% missed | 0.56 / 0.33 | 0.51 / 0.20 |
| nothing seen below 600 Da | 0 / 0 | 0 / 0 |
| 2000 noise masses | 0.71 / 0.41 | 0.85 / 0.85 |

What it says, in the paper's own terms:

* **Mass accuracy decides it.** The paper needed +-0.01 a.m.u. for a 150-residue protein; at the
  masses a modern digest yields, that is a few ppm. At 1-3 ppm (any current Orbitrap or FT-ICR) the
  method works; at 10 ppm it barely does; at 30 ppm (an older TOF) not at all.
* **Completeness is the other half.** The method lives on fragments that overlap and nest. With
  everything observed, 86-91% of each protein comes back; with 40% of fragments missed it falls
  apart. Short fragments matter most: raising the lowest observable mass from 350 to 600 Da leaves
  nothing.
* **Noise barely matters**, as the paper predicted ("random noise masses can be tolerated in
  surprisingly large quantity"): ten times more noise masses than the baseline change almost nothing.

## Real data

`masspec run` takes a centroided mzML directly: `features.py` finds the isotope envelopes in every
MS1 scan and merges them into one list of neutral masses, MASSPEC's "single spectrum". It also reads
mass lists from other feature finders (MaxQuant `allPeptides.txt`, Dinosaur `.features.tsv`, any
CSV/TSV with a mass column).

The test case is the undeuterated pepsin digest of sperm whale myoglobin in
[PXD071026](https://www.ebi.ac.uk/pride/archive/projects/PXD071026), an HDX-MS study on an LTQ
Orbitrap: one small protein, a non-specific enzyme, one run. The protein's sequence is in the
deposit (`rcsb_pdb_6Z4T.fasta`). `scripts/real_pxd071026.py` reproduces every number below; its
output is in `docs/real_pxd071026.txt`. The run's merged mass list ships with the package
(`masspec.example.pxd071026()`, CC0 like the deposit; provenance in `src/masspec/data/README.md`),
so the notebook's real-data section runs without the 148 MB raw file, on Colab too.

```bash
mono ThermoRawFileParser.exe -i=data/raw/Mn_0s_1.raw -o=data/raw/ -f=1   # Thermo .raw to centroided mzML
uv run masspec run data/raw/Mn_0s_1.mzML --ppm 3 --no-alkylation \
    --calibrant data/raw/pepsin_A_mature.fasta --truth data/raw/rcsb_pdb_6Z4T.fasta
```

**Calibration without the answer.** The run's masses sit about 5 ppm low. In a pepsin digest pepsin
digests itself, its sequence is public (UniProt P00791), and its autolysis fragments are in every
run, so `--calibrant` estimates the offset from them, never from the protein being sequenced:
-5.06 ppm, with 285 pepsin matches in the peak 2-ppm bin against about 54 in a background bin.

**What the run contains.** 5440 MS1 features. After calibration, 180 match a myoglobin fragment
within 2 ppm, of which about 35 would match by chance, so roughly 145 real myoglobin fragments,
mostly long ones (up to 49 residues). Pepsin's own fragments account for 408.

**What MASSPEC makes of it: nothing.** No residue of myoglobin (or of pepsin) is recovered, at 3 or
5 ppm. Removing every noise mass does not help: on the 180 myoglobin-matching masses alone, recall
is 6%, which is what the same contigs score against shuffled myoglobin (up to 6% over ten
shuffles): chance.

**What a present-day *de novo* model makes of the same run.** The run also recorded 4505 MS/MS
spectra, which MASSPEC never had: each isolates one peptide and fragments it. InstaNovo 1.2.2
(model instanovo-v1.2.0, with InstaNovo+ refinement, on a GPU) reads them; kept at 10 ppm and
confidence 0.9, its 248 distinct peptides include 33 exact stretches of myoglobin, covering 82% of
it, and 124 of pepsin. Against shuffled myoglobin no prediction is an exact stretch. The predictions
ship with the package (`masspec.example.pxd071026_instanovo()`) and close the notebook;
`scripts/package_instanovo_pxd071026.py` gives the command that made them.

The study above explains why MASSPEC fails here. The simulated pepsin digest of the same protein yields about 2000
detectable fragments; this run holds about 145, some 7%. That is the "fragments missed" axis taken
far past the 60% at which the method had already collapsed. Dayhoff and Eck's method needs many
overlapping and nested fragments, including short ones: about six masses per residue (302 for
insulin's 51) in their insulin problem. An HDX peptide map is built for the opposite purpose: a
modest set of well-separated peptides that tile the protein once.

So does a method that rests on accurate mass work better now that instruments measure mass so much
better? Only **if the experiment is designed for it**. The accuracy is there: a few ppm, after
calibration. What is missing is the data: a near-random partial digest, with the short fragments
kept. A standard proteomics or HDX run is not that experiment.

### Partial acid hydrolysis (PXD003804, PXD030203)

The digest the paper had in mind is partial acid hydrolysis. Two PRIDE deposits come close: horse
myoglobin (153 residues) hydrolysed in HCl in a microwave oven, then LC-MS/MS on an Orbitrap, both
made for *de novo* protein sequencing from MS/MS:
[PXD003804](https://www.ebi.ac.uk/pride/archive/projects/PXD003804) (Savidor et al., *MCP* 2017,
[doi:10.1074/mcp.O116.065417](https://doi.org/10.1074/mcp.O116.065417); 6 M HCl, Q Exactive, MS1
from 300 m/z) and [PXD030203](https://www.ebi.ac.uk/pride/archive/projects/PXD030203) (Mai et al.,
*Anal. Chem.* 2022, [doi:10.1021/acs.analchem.1c03718](https://doi.org/10.1021/acs.analchem.1c03718);
Orbitrap Fusion Lumos, MS1 from 350 m/z). Both deposits include PEAKS's *de novo* reading of every MS/MS
spectrum. The peptides read with an average local confidence of 50% or more, with their measured
precursor masses, ship with the package (`masspec.example.acid_hydrolysates()`, written by
`scripts/package_acid_hydrolysis.py`). The notebook's section on them gives every number below.

Keeping only the reads that are exact stretches of myoglobin gives MASSPEC its best case: every mass
a real fragment, no noise.

| | PXD003804 | PXD030203 | simulated acid hydrolysate |
|---|---|---|---|
| distinct myoglobin fragments | 739 | 1092 | 2364 |
| per residue | 4.8 | 7.1 | 15.5 |
| shortest | 3 residues | 3 residues | 3 residues |
| lighter than 500 Da | 20 | 25 | 163 |
| MASSPEC recall / precision, as deposited | 6% / 6% | 7% / 6% | |
| + every 2-4-residue fragment (1.5 ppm error) | 41% / 30% | 51% / 81% | |
| + the same with a fifth of them missed | 23% / 13% | 33% / 46% | |

Acid gives what pepsin did not: five to seven fragments per residue, about what the paper's insulin
problem had, starting at nearly every residue. But the short end is still missing: the MS1 scans
start at 300 or 350 m/z, the cleanup (solid-phase extraction) and the C18 column lose small polar
peptides, and the instrument rarely chooses singly charged ions for MS/MS. As deposited, MASSPEC
recovers 6-7% of the protein, where the same contigs score up to 4% against shuffled myoglobin.
Adding back the short fragments that a run with a low m/z limit would see (made up from the known
sequence, so a "what if") brings back about half. The missing experiment is a partial acid
hydrolysate measured for this: MS1 down to about 100 m/z, no cleanup step that loses small polar
peptides, and singly charged ions kept.

## Layout

```
src/masspec/masses.py     residue masses: the paper's table, and present-day monoisotopic masses
src/masspec/fragments.py  FRAGEN and MASFRG's data half: the paper's ideal data
src/masspec/mspec.py      the 1970 reconstruction
src/masspec/simulate.py   present-day data: enzyme, detection window, dropout, ppm, mods, noise
src/masspec/modern.py     the reconstruction for present-day data
src/masspec/features.py   mzML to one mass list; readers for other tools' mass lists
src/masspec/evaluate.py   picking the answer, and scoring it against a known sequence
scripts/study.py          the table above
```

## License

MIT, see [LICENSE](LICENSE). The license covers this code, not the 1970 paper:
nothing of the paper is reproduced here beyond short quotations and the
published masses and counts used to check the reimplementation.

## Why?

Prompted by
[issue #4](https://github.com/BioGeek/awesome_de_novo_peptide_sequencing/issues/4) of the awesome
*de novo* peptide sequencing catalog.