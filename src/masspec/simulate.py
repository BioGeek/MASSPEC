"""Present-day data for MASSPEC: the peptide masses an LC-MS run of a non-specific digest would give.

Dayhoff and Eck simulated ideal data: every fragment up to nine residues, exact masses, nothing
missing. A modern experiment differs in ways that all matter to the method:

* Cleavage is not uniform. Pepsin prefers to cut after F, L, Y, W and M; proteinase K after
  aliphatic and aromatic residues; dilute-acid hydrolysis after D. So some fragments are abundant
  and many are rare.
* Detection has a window. Short peptides (roughly below 300-400 Da) do not retain on a C18 column
  and are seldom seen, while long ones (up to several kDa) are, since the instrument sees every
  charge state. Low-abundance fragments drop out.
* Masses are measured to a few ppm, not to 0.01-0.02 Da.
* Cysteines are usually carbamidomethylated, and some methionines are oxidised.
* There is noise: other proteins, isotope peaks assigned as monoisotopic (+1.00336 Da), and unrelated
  features.

`simulate()` produces such a mass list from a protein sequence, with the truth kept beside it so a
reconstruction can be scored.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field

from .masses import MASSES, MODIFIED, H2O, ISOTOPE_SPACING, normalise

# Relative cleavage rate after each residue (the C-terminal side of the bond), per enzyme. The
# values are deliberately coarse: what matters to MASSPEC is that cleavage is uneven, not its
# exact preferences.
ENZYMES = {
    "random": {},                                        # every bond equally likely
    "pepsin": {"F": 1.0, "L": 1.0, "Y": 0.8, "W": 0.8, "M": 0.6, "E": 0.3},
    "proteinase_k": {"A": 0.8, "F": 1.0, "Y": 1.0, "W": 1.0, "L": 0.9, "V": 0.8, "M": 0.6},
    "acid": {"D": 1.0},                                  # dilute acid: Asp-X bonds first
}


@dataclass
class Simulated:
    protein: str
    masses: list                      # what the "instrument" reports
    peptides: list                    # the true peptides detected, (start, end, sequence, abundance)
    params: dict = field(default_factory=dict)


def bond_probs(seq, enzyme="pepsin", p_max=0.35, background=0.03):
    """Probability that each internal bond (after residue k) is cut in one protein molecule."""
    pref = ENZYMES[enzyme]
    if not pref:
        return [p_max] * (len(seq) - 1)
    return [background + (p_max - background) * pref.get(seq[k], 0.0) for k in range(len(seq) - 1)]


def peptide_abundances(seq, probs, min_len=2, max_len=40):
    """Expected fraction of molecules yielding each fragment i..j (inclusive), for every fragment.

    A fragment survives when the bonds at its two ends are cut (a protein terminus counts as cut)
    and every bond inside it is intact.
    """
    n = len(seq)
    out = []
    for i in range(n):
        left = 1.0 if i == 0 else probs[i - 1]
        intact = 1.0
        for j in range(i, min(n, i + max_len)):
            if j > i:
                intact *= 1.0 - probs[j - 1]
            right = 1.0 if j == n - 1 else probs[j]
            if j - i + 1 >= min_len:
                out.append((i, j, left * intact * right))
    return out


def simulate(protein, enzyme="pepsin", p_max=0.35, background=0.03, min_len=1, max_len=40,
             mass_lo=350.0, mass_hi=6000.0, rel_threshold=1e-3, dropout=0.2, ppm=3.0,
             alkylated=True, met_ox=0.3, isotope_errors=0.05, n_noise=200, seed=0):
    """A present-day mass list for one protein.

    rel_threshold: fragments below this fraction of the most abundant one are not detected.
    dropout: of the rest, this fraction is missed anyway (ionisation, co-elution).
    met_ox: fraction of Met-containing fragments that are ALSO seen oxidised.
    isotope_errors: fraction of detected fragments ALSO reported 1.00336 Da heavy.
    n_noise: unrelated masses, uniform over the window.
    """
    rng = random.Random(seed)
    seq = normalise(protein)
    probs = bond_probs(seq, enzyme, p_max, background)
    ab = peptide_abundances(seq, probs, min_len, max_len)
    top = max(a for _, _, a in ab)
    table = dict(MASSES)
    if alkylated:
        table["C"] = MODIFIED["c"]
    masses, detected = [], []
    for i, j, a in ab:
        if a < rel_threshold * top or rng.random() < dropout:
            continue
        pep = seq[i:j + 1]
        m = sum(table[x] for x in pep) + H2O
        if not mass_lo <= m <= mass_hi:
            continue
        detected.append((i, j, pep, a / top))
        variants = [m]
        if "M" in pep and rng.random() < met_ox:
            variants.append(m + 15.994915)
        if rng.random() < isotope_errors:
            variants.append(m + ISOTOPE_SPACING)
        for v in variants:
            masses.append(v * (1.0 + rng.gauss(0.0, ppm) * 1e-6))
    masses += [rng.uniform(mass_lo, mass_hi) for _ in range(n_noise)]
    params = dict(enzyme=enzyme, p_max=p_max, background=background, min_len=min_len,
                  max_len=max_len, mass_lo=mass_lo, mass_hi=mass_hi, rel_threshold=rel_threshold,
                  dropout=dropout, ppm=ppm, alkylated=alkylated, met_ox=met_ox,
                  isotope_errors=isotope_errors, n_noise=n_noise, seed=seed)
    return Simulated(seq, sorted(masses), detected, params)


def coverage(sim: Simulated, max_len=None):
    """Of all fragments up to `max_len` inside the window, the fraction that was detected."""
    n = len(sim.protein)
    max_len = max_len or sim.params["max_len"]
    lo, hi = sim.params["mass_lo"], sim.params["mass_hi"]
    table = dict(MASSES)
    if sim.params["alkylated"]:
        table["C"] = MODIFIED["c"]
    total = sum(1 for i in range(n) for j in range(i, min(n, i + max_len))
                if lo <= sum(table[x] for x in sim.protein[i:j + 1]) + H2O <= hi)
    seen = sum(1 for i, j, _, _ in sim.peptides if j - i + 1 <= max_len)
    return seen / total if total else math.nan
