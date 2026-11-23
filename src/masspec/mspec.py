"""MASSPEC's reconstruction, as Dayhoff and Eck describe it.

Comput. Biol. Med. 1, 5-28, 1970.

The input is a list of observed peptide masses: every fragment of a protein's
near-random partial hydrolysate, measured together in one high-resolution mass
spectrum. The output is the protein's sequence, found as the one set of
overlapping fragments that stays mutually consistent among the many spurious
compositions the same masses also allow.

The programs, in the paper's order (its Appendix 1):

* ``MASFRG`` (reconstruction half): every amino acid COMPOSITION up to four
  residues long whose mass is within the tolerance of an observed mass. Order
  is unknown at this point, and many compositions share a mass (AV = GL), so
  the list is mostly spurious.
* ``PRUNE``: a composition of N residues is kept only if at least two of the
  compositions you get by removing one residue are in the list of length N-1.
  "If the peptide (F,G,H,I) is true in any of its possible permutations, there
  must also be at least two of the following: (F,G,H), (F,G,I), (F,H,I),
  (G,H,I)."
* ``WEED``: a composition is kept only if it is contained in at least one kept
  composition one residue longer, longest first.
* ``CQ``: the compositions become SEQUENCES; a sequence of length N is kept
  only if both of the length N-1 sequences inside it (its prefix and its
  suffix) are kept. Masses carry no direction, so a sequence and its reverse
  are the same entry ("AEJ means both AEJ and JEA").
* ``OVERLAP``: from length 5 on, two kept sequences of length N-1 sharing N-2
  residues overlap into a candidate of length N, which is kept only if its
  mass was observed. "Masses above a present limit are all assumed to have
  been observed." A sequence that overlaps no other cannot be extended: it is
  a complete chain, and is reported.

``SEQ``, the paper's tree search written to save core memory, computes the
same answer as repeated OVERLAP and is not reproduced.
"""
from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from itertools import combinations_with_replacement, permutations

from .masses import H2O, MASSES


def canon(seq: str) -> str:
    """A sequence and its reverse are one entry.

    Keep the alphabetically smaller spelling.
    """
    reverse = seq[::-1]
    return seq if seq <= reverse else reverse


class Observed:
    """The observed mass list.

    Searchable within an absolute (Da) or relative (ppm) tolerance.
    """

    def __init__(self, masses, tol_da=None, tol_ppm=None):
        if (tol_da is None) == (tol_ppm is None):
            raise ValueError("give exactly one of tol_da and tol_ppm")
        self.masses = sorted(float(mass) for mass in masses)
        self.tol_da, self.tol_ppm = tol_da, tol_ppm

    def window(self, mass):
        if self.tol_da is not None:
            return self.tol_da
        return mass * self.tol_ppm * 1e-6

    def contains(self, mass):
        width = self.window(mass)
        index = bisect.bisect_left(self.masses, mass - width)
        return (index < len(self.masses)
                and self.masses[index] <= mass + width)

    def __len__(self):
        return len(self.masses)


def _longest_first(seq):
    """Sort key: longer sequences first, then alphabetical."""
    return (-len(seq), seq)


@dataclass
class Result:
    # (sequence, length) of chains nothing extends
    complete: list = field(default_factory=list)
    # (stage, length, count)
    history: list = field(default_factory=list)
    converged: bool = True
    reason: str = ""

    def chains(self, min_len=1):
        """Complete chains at least `min_len` long, longest first."""
        return sorted({seq for seq in self.complete if len(seq) >= min_len},
                      key=_longest_first)


def masfrg(observed: Observed, table=MASSES, water=H2O, max_len=4):
    """Every composition up to `max_len` residues matching an observed mass."""
    alphabet = sorted(table)
    comps = {}
    for length in range(1, max_len + 1):
        comps[length] = {
            comp
            for comp in combinations_with_replacement(alphabet, length)
            if observed.contains(sum(table[residue] for residue in comp)
                                 + water)
        }
    return comps


def _drop_one(comp):
    return [comp[:index] + comp[index + 1:] for index in range(len(comp))]


def prune(comps, min_subs=2):
    """PRUNE: keep a composition only if enough sub-compositions are kept.

    At least `min_subs` of its one-shorter sub-compositions must be kept.
    """
    out = {1: set(comps[1])}
    for length in range(2, max(comps) + 1):
        out[length] = {
            comp for comp in comps[length]
            if sum(sub in out[length - 1] for sub in _drop_one(comp))
            >= min_subs
        }
    return out


def weed(comps):
    """WEED: keep a composition only if a longer kept one contains it.

    The containing composition is one residue longer.
    """
    top = max(comps)
    out = {top: set(comps[top])}
    for length in range(top - 1, 0, -1):
        inside = {sub for comp in out[length + 1] for sub in _drop_one(comp)}
        out[length] = {comp for comp in comps[length] if comp in inside}
    return out


def cq(comps):
    """CQ: compositions to sequences whose prefix and suffix are both kept."""
    seqs = {1: {comp[0] for comp in comps[1]}}
    for length in range(2, max(comps) + 1):
        keep = set()
        for comp in comps[length]:
            for perm in set(permutations(comp)):
                seq = "".join(perm)
                if (canon(seq[:-1]) in seqs[length - 1]
                        and canon(seq[1:]) in seqs[length - 1]):
                    keep.add(canon(seq))
        seqs[length] = keep
    return seqs


def overlap_step(prev, length, accept):
    """OVERLAP: length-`length` candidates from pairs of shorter sequences.

    The pairs are of length `length`-1 sequences sharing `length`-2 residues.
    `accept(seq)` decides whether a candidate stands. Returns (kept,
    complete): the kept length-`length` sequences, and the length `length`-1
    sequences that took part in none of them.
    """
    oriented = set(prev) | {seq[::-1] for seq in prev}
    by_head = {}
    for tail in oriented:
        by_head.setdefault(tail[:-1], []).append(tail)
    kept, used = set(), set()
    for head in oriented:
        for tail in by_head.get(head[1:], ()):
            new = head + tail[-1]
            if accept(new):
                kept.add(canon(new))
                used.add(canon(head))
                used.add(canon(tail))
    return kept, set(prev) - used


def run(observed: Observed, table=MASSES, water=H2O, mass_limit=1000.0,
        max_len=4, min_subs=2, use_weed=False, max_protein=2000,
        max_list=200_000) -> Result:
    """The whole of MSPEC: MASFRG, PRUNE, WEED, CQ, then OVERLAP.

    OVERLAP repeats until nothing extends.
    """
    res = Result()
    comps = masfrg(observed, table, water, max_len)
    for length in comps:
        res.history.append(("MASFRG", length, len(comps[length])))
    comps = prune(comps, min_subs)
    for length in comps:
        res.history.append(("PRUNE", length, len(comps[length])))
    if use_weed:
        comps = weed(comps)
        for length in comps:
            res.history.append(("WEED", length, len(comps[length])))
    seqs = cq(comps)
    for length in seqs:
        res.history.append(("CQ", length, len(seqs[length])))
    # A sequence that nothing longer contains is already complete (a short
    # chain, or noise).
    for length in range(1, max_len):
        inside = ({canon(seq[:-1]) for seq in seqs[length + 1]}
                  | {canon(seq[1:]) for seq in seqs[length + 1]})
        res.complete += [seq for seq in seqs[length] if seq not in inside]

    def accept(seq):
        mass = sum(table[residue] for residue in seq) + water
        return mass > mass_limit or observed.contains(mass)

    cur, length = seqs[max_len], max_len + 1
    while cur and length <= max_protein:
        if len(cur) > max_list:
            res.converged, res.reason = False, (
                f"{len(cur)} candidate sequences of length {length - 1}: "
                "the list did not shrink")
            break
        cur, done = overlap_step(cur, length, accept)
        res.history.append(("OVERLAP", length, len(cur)))
        res.complete += sorted(done)
        length += 1
    return res
