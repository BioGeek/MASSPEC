"""FRAGEN and the simulation half of MASFRG: ideal mass-spectrometer data, as the paper made it.

FRAGEN "generates short peptides from polymers ... the results of a random acid hydrolysis may be
simulated; all possible peptide compositions of specified lengths may be generated". MASFRG then
"produces the masses of peptides and of noise which might be obtained from a mass spectrometer".
The insulin experiment used every fragment of both chains up to 9 residues long, keeping masses
below 1000.
"""
from __future__ import annotations

import random

from .masses import MASSES, H2O, peptide_mass


def fragen(chains, max_len=9):
    """Every contiguous fragment of every chain, 1 to `max_len` residues long (with repeats)."""
    out = []
    for chain in chains:
        for n in range(1, max_len + 1):
            out += [chain[i:i + n] for i in range(len(chain) - n + 1)]
    return out


def ideal_masses(chains, max_len=9, mass_limit=1000.0, table=MASSES, water=H2O):
    """The paper's ideal data: the distinct masses of all fragments below `mass_limit`."""
    return sorted({round(peptide_mass(f, table, water), 6) for f in fragen(chains, max_len)
                   if peptide_mass(f, table, water) < mass_limit})


def add_noise(masses, n_random=0, lo=100.0, hi=1000.0, isotopes=False, seed=0):
    """MASFRG's noise: `n_random` uniform random masses, and optionally each mass's 13C isotope peak."""
    rng = random.Random(seed)
    out = list(masses)
    if isotopes:
        out += [m + 1.0033548378 for m in masses]
    out += [rng.uniform(lo, hi) for _ in range(n_random)]
    return sorted(out)
