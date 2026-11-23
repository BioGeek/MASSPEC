"""Residue masses.

MASSPEC (Dayhoff and Eck, Comput. Biol. Med. 1, 5-28, 1970) tabulates the masses of the FREE
amino acids (its Fig. 2) and computes a peptide's mass as their sum minus one water per peptide
bond. That is the same as the sum of residue masses plus one water, which is how this package
computes it: `peptide_mass(seq) = sum(residue) + H2O`.

Two tables:

* ``MASSES_1970``: the paper's Fig. 2, converted to residue masses (free mass minus 18.01056),
  in the paper's own one-letter code: J is Leu/Ile (the method cannot tell them apart), O is Tyr.
  Used to reproduce the paper's insulin experiment number for number.
* ``MASSES``: present-day monoisotopic residue masses (Unimod / NIST), standard one-letter code,
  with L standing for Leu/Ile. Used for everything else.
"""

H2O_1970 = 18.01056
H2O = 18.0105646863

# Fig. 2 of the paper: free amino acid masses (a.m.u.).
_FREE_1970 = {
    "A": 89.04763, "C": 121.01970, "D": 133.03745, "E": 147.05309, "F": 165.07891,
    "G": 75.03199, "H": 155.06941, "J": 131.09455, "K": 146.10544, "M": 149.05098,
    "N": 132.05343, "O": 181.07382, "P": 115.06327, "Q": 146.06907, "R": 174.11158,
    "S": 105.04254, "T": 119.05818, "V": 117.07891, "W": 188.09489,
}
# W is as PRINTED, 188.09489, although free tryptophan weighs 204.0899 (C11H12N2O2): the paper's
# own counts need it. With the true value, W's residue (186.079) sits within 0.02 of the dipeptide
# G+E (186.064), so every "GE" in insulin also reads as a W: the lists come out at 17 single residues
# and 71 dipeptide compositions instead of the paper's 16 and 66, and a second B chain with "W" for
# "EG" survives. With the printed value both disappear. Insulin has no tryptophan, so the misprint
# was invisible to its authors.
# Not rounded: the paper computes in free-acid masses (sum minus water per bond), so rounding each
# residue would put the last digit off by one.
MASSES_1970 = {aa: m - H2O_1970 for aa, m in _FREE_1970.items()}

MASSES = {
    "G": 57.02146372, "A": 71.03711379, "S": 87.03202841, "P": 97.05276385, "V": 99.06841391,
    "T": 101.0476785, "C": 103.0091845, "L": 113.0840640, "N": 114.0429275, "D": 115.0269431,
    "Q": 128.0585775, "K": 128.0949630, "E": 129.0425931, "M": 131.0404846, "H": 137.0589119,
    "F": 147.0684139, "R": 156.1011110, "Y": 163.0633285, "W": 186.0792630,
}

# Common modifications, as extra residue letters the modern mode can switch on. Lower case so they
# never collide with an amino acid.
MODIFIED = {
    "c": MASSES["C"] + 57.021464,   # carbamidomethyl Cys (fixed, after alkylation)
    "m": MASSES["M"] + 15.994915,   # oxidised Met
    "n": MASSES["N"] + 0.984016,    # deamidated Asn
    "q": MASSES["Q"] + 0.984016,    # deamidated Gln
}

ISOTOPE_SPACING = 1.0033548378  # 13C - 12C


def peptide_mass(seq, table=MASSES, water=H2O):
    """Neutral monoisotopic mass of a linear peptide."""
    return sum(table[a] for a in seq) + water


def normalise(seq):
    """A sequence in the modern alphabet with I folded into L (mass cannot separate them)."""
    return seq.upper().replace("I", "L")
