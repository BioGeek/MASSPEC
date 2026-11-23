"""MASSPEC's reconstruction, as Dayhoff and Eck describe it (Comput. Biol. Med. 1, 5-28, 1970).

The input is a list of observed peptide masses: every fragment of a protein's near-random partial
hydrolysate, measured together in one high-resolution mass spectrum. The output is the protein's
sequence, found as the one set of overlapping fragments that stays mutually consistent among the
many spurious compositions the same masses also allow.

The programs, in the paper's order (its Appendix 1):

* ``MASFRG`` (reconstruction half): every amino acid COMPOSITION up to four residues long whose
  mass is within the tolerance of an observed mass. Order is unknown at this point, and many
  compositions share a mass (AV = GL), so the list is mostly spurious.
* ``PRUNE``: a composition of N residues is kept only if at least two of the compositions you get
  by removing one residue are in the list of length N-1. "If the peptide (F,G,H,I) is true in any of
  its possible permutations, there must also be at least two of the following: (F,G,H), (F,G,I),
  (F,H,I), (G,H,I)."
* ``WEED``: a composition is kept only if it is contained in at least one kept composition one
  residue longer, longest first.
* ``CQ``: the compositions become SEQUENCES; a sequence of length N is kept only if both of the
  length N-1 sequences inside it (its prefix and its suffix) are kept. Masses carry no direction,
  so a sequence and its reverse are the same entry ("AEJ means both AEJ and JEA").
* ``OVERLAP``: from length 5 on, two kept sequences of length N-1 sharing N-2 residues overlap into
  a candidate of length N, which is kept only if its mass was observed. "Masses above a present limit
  are all assumed to have been observed." A sequence that overlaps no other cannot be extended: it is
  a complete chain, and is reported.

``SEQ``, the paper's tree search written to save core memory, computes the same answer as repeated
OVERLAP and is not reproduced.
"""
from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from itertools import combinations_with_replacement, permutations

from .masses import MASSES, H2O


def canon(seq: str) -> str:
    """A sequence and its reverse are one entry: keep the alphabetically smaller spelling."""
    r = seq[::-1]
    return seq if seq <= r else r


class Observed:
    """The observed mass list, searchable within an absolute (Da) or relative (ppm) tolerance."""

    def __init__(self, masses, tol_da=None, tol_ppm=None):
        if (tol_da is None) == (tol_ppm is None):
            raise ValueError("give exactly one of tol_da and tol_ppm")
        self.masses = sorted(float(m) for m in masses)
        self.tol_da, self.tol_ppm = tol_da, tol_ppm

    def window(self, m):
        return self.tol_da if self.tol_da is not None else m * self.tol_ppm * 1e-6

    def contains(self, m):
        w = self.window(m)
        i = bisect.bisect_left(self.masses, m - w)
        return i < len(self.masses) and self.masses[i] <= m + w

    def __len__(self):
        return len(self.masses)


@dataclass
class Result:
    complete: list = field(default_factory=list)       # (sequence, length) of chains nothing extends
    history: list = field(default_factory=list)        # (stage, length, count)
    converged: bool = True
    reason: str = ""

    def chains(self, min_len=1):
        """Complete chains at least `min_len` long, longest first."""
        return sorted({s for s in self.complete if len(s) >= min_len}, key=lambda s: (-len(s), s))


def masfrg(observed: Observed, table=MASSES, water=H2O, max_len=4):
    """Every composition up to `max_len` residues whose mass matches an observed mass."""
    alphabet = sorted(table)
    comps = {}
    for n in range(1, max_len + 1):
        comps[n] = {c for c in combinations_with_replacement(alphabet, n)
                    if observed.contains(sum(table[a] for a in c) + water)}
    return comps


def _drop_one(c):
    return [c[:i] + c[i + 1:] for i in range(len(c))]


def prune(comps, min_subs=2):
    """PRUNE: keep a composition only if `min_subs` of its one-shorter sub-compositions are kept."""
    out = {1: set(comps[1])}
    for n in range(2, max(comps) + 1):
        out[n] = {c for c in comps[n] if sum(s in out[n - 1] for s in _drop_one(c)) >= min_subs}
    return out


def weed(comps):
    """WEED: keep a composition only if some kept composition one residue longer contains it."""
    top = max(comps)
    out = {top: set(comps[top])}
    for n in range(top - 1, 0, -1):
        inside = {s for c in out[n + 1] for s in _drop_one(c)}
        out[n] = {c for c in comps[n] if c in inside}
    return out


def cq(comps):
    """CQ: compositions to sequences whose prefix and suffix are both kept."""
    seqs = {1: {c[0] for c in comps[1]}}
    for n in range(2, max(comps) + 1):
        keep = set()
        for c in comps[n]:
            for p in set(permutations(c)):
                s = "".join(p)
                if canon(s[:-1]) in seqs[n - 1] and canon(s[1:]) in seqs[n - 1]:
                    keep.add(canon(s))
        seqs[n] = keep
    return seqs


def overlap_step(prev, n, accept):
    """OVERLAP: length-n candidates from pairs of length n-1 sequences sharing n-2 residues.

    `accept(seq)` decides whether a candidate stands. Returns (kept, complete): the kept length-n
    sequences, and the length n-1 sequences that took part in none of them.
    """
    oriented = set(prev) | {s[::-1] for s in prev}
    by_head = {}
    for t in oriented:
        by_head.setdefault(t[:-1], []).append(t)
    kept, used = set(), set()
    for s in oriented:
        for t in by_head.get(s[1:], ()):
            new = s + t[-1]
            if accept(new):
                kept.add(canon(new))
                used.add(canon(s))
                used.add(canon(t))
    return kept, set(prev) - used


def run(observed: Observed, table=MASSES, water=H2O, mass_limit=1000.0, max_len=4,
        min_subs=2, use_weed=False, max_protein=2000, max_list=200_000) -> Result:
    """The whole of MSPEC: MASFRG, PRUNE, WEED, CQ, then OVERLAP until nothing extends."""
    res = Result()
    comps = masfrg(observed, table, water, max_len)
    for n in comps:
        res.history.append(("MASFRG", n, len(comps[n])))
    comps = prune(comps, min_subs)
    for n in comps:
        res.history.append(("PRUNE", n, len(comps[n])))
    if use_weed:
        comps = weed(comps)
        for n in comps:
            res.history.append(("WEED", n, len(comps[n])))
    seqs = cq(comps)
    for n in seqs:
        res.history.append(("CQ", n, len(seqs[n])))
    # A sequence that nothing longer contains is already complete (a short chain, or noise).
    for n in range(1, max_len):
        inside = {canon(s[:-1]) for s in seqs[n + 1]} | {canon(s[1:]) for s in seqs[n + 1]}
        res.complete += [s for s in seqs[n] if s not in inside]

    def accept(seq):
        m = sum(table[a] for a in seq) + water
        return m > mass_limit or observed.contains(m)

    cur, n = seqs[max_len], max_len + 1
    while cur and n <= max_protein:
        if len(cur) > max_list:
            res.converged, res.reason = False, f"{len(cur)} candidate sequences of length {n - 1}: the list did not shrink"
            break
        cur, done = overlap_step(cur, n, accept)
        res.history.append(("OVERLAP", n, len(cur)))
        res.complete += sorted(done)
        n += 1
    return res
