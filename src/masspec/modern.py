"""MASSPEC taken to present-day data.

The 1970 method keeps a candidate only when EVERY fragment inside it was observed, and treats masses
above the instrument's limit as observed. Both assumptions fail on an LC-MS run of a non-specific
digest: fragments below ~350 Da are hardly ever seen, and many inside the window are missed.

This keeps the paper's phases and data structures, and changes three things, each a direct
generalisation of something in the paper:

1. **An observable window at both ends.** The paper assumed every mass above 1000 a.m.u. observed,
   because its instrument could not measure it. Here masses outside ``[lo, hi]`` are assumed
   observed for the same reason: below ``lo`` the fragment cannot be seen, so its absence says
   nothing. That lets the bottom-up phases (MASFRG, PRUNE, CQ) start without the short fragments.
2. **Support instead of completeness.** A candidate of length N is judged on the fragments that END
   at its last residue (its suffixes): exactly the fragments the extension introduces. The paper
   requires all of them observed; here it is enough that a fraction ``support`` of those with a mass
   in the window, and at most ``support_len`` residues long, were observed. Long fragments are rare in
   a real digest, so counting them would sink true sequences.
3. **A beam.** With support below 1, spurious candidates survive more easily, and the lists can grow
   without bound. The paper hit the same wall ("the storage is exceeded and the problem must be
   terminated"). Here the ``beam`` best-supported candidates of each length are kept.

Masses are matched within a ppm tolerance. Modifications enter as extra residue letters (see
`masses.MODIFIED`). The result is a set of contigs: with incomplete data a chain breaks where the
digest left no overlapping evidence, so the protein comes back in pieces whose order masses alone
cannot give.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations_with_replacement

from .masses import MASSES, H2O
from .mspec import Observed, canon, cq, prune


@dataclass
class ModernResult:
    contigs: list = field(default_factory=list)     # (sequence, support score), longest first
    history: list = field(default_factory=list)
    converged: bool = True
    reason: str = ""


class FastObserved:
    """Observed masses binned for O(1) lookup within a ppm tolerance."""

    def __init__(self, masses, tol_ppm):
        self.tol_ppm = tol_ppm
        self.masses = sorted(masses)
        self.width = max(self.masses) * tol_ppm * 1e-6 if self.masses else 1.0
        self.bins = {}
        for m in self.masses:
            self.bins.setdefault(int(m / self.width), []).append(m)

    def contains(self, m):
        w = m * self.tol_ppm * 1e-6
        k = int(m / self.width)
        return any(abs(x - m) <= w for b in (k - 1, k, k + 1) for x in self.bins.get(b, ()))

    def __len__(self):
        return len(self.masses)


def run_extend(observed, table=MASSES, water=H2O, lo=350.0, hi=6000.0, seed_len=4, support=0.5,
               support_len=20, min_hits=1, beam=2000, min_contig=6, max_protein=3000) -> ModernResult:
    """Grow chains one residue at a time instead of overlapping fixed-length windows.

    OVERLAP builds a length-N candidate only from two surviving length N-1 windows, so one missing
    window anywhere breaks every chain through it. Here a chain is extended by any residue whose
    suffixes (the fragments ending at the new residue) are supported: at least `min_hits` observed,
    and at least a fraction `support` of those in the window. Both directions are grown, since a
    seed may sit anywhere in the protein. The seeds are MASSPEC's own: CQ's sequences of
    `seed_len` residues whose mass, when in the window, was observed.
    """
    res = ModernResult()
    letters = sorted(table)
    mass = lambda s: sum(table[a] for a in s) + water
    seen = lambda m: m < lo or m > hi or observed.contains(m)
    comps = {n: {c for c in combinations_with_replacement(letters, n) if seen(mass(c))}
             for n in range(1, seed_len + 1)}
    seqs = cq(prune(comps, 2))
    seeds = {s for s in seqs[seed_len] if seen(mass(s))}
    cur = {}
    for s in seeds:
        cur[s] = 0
        cur[s[::-1]] = 0
    res.history.append(("SEEDS", seed_len, len(seeds)))

    def suffix_support(seq):
        o = t = 0
        m = water
        for k in range(1, min(len(seq), support_len) + 1):
            m += table[seq[-k]]
            if lo <= m <= hi:
                t += 1
                o += observed.contains(m)
        return o, t

    n = seed_len + 1
    while cur and n <= max_protein:
        nxt = {}
        for s, sc in cur.items():
            grown = False
            for a in letters:
                new = s + a
                o, t = suffix_support(new)
                if o >= min_hits and o >= support * t:
                    nxt[new] = max(nxt.get(new, -1), sc + o)
                    grown = True
            if not grown and len(s) >= min_contig:
                res.contigs.append((canon(s), sc))
        if len(nxt) > beam:
            res.converged = False
            res.reason = f"beam of {beam} reached"
            nxt = dict(sorted(nxt.items(), key=lambda kv: -kv[1])[:beam])
        res.history.append(("EXTEND", n, len(nxt)))
        cur = nxt
        n += 1
    best = {}
    for s, sc in res.contigs:
        best[s] = max(best.get(s, -1), sc)
    res.contigs = sorted(best.items(), key=lambda cs: (-cs[1], -len(cs[0])))
    return res


def run_modern(observed: Observed, table=MASSES, water=H2O, lo=350.0, hi=6000.0, seed_len=4,
               support=0.6, support_len=12, beam=5000, min_contig=6, max_protein=3000) -> ModernResult:
    res = ModernResult()
    letters = sorted(table)
    mass = lambda s: sum(table[a] for a in s) + water

    def seen(m):
        return m < lo or m > hi or observed.contains(m)

    # MASFRG with the window: compositions whose mass is observed, or cannot be.
    comps = {n: {c for c in combinations_with_replacement(letters, n) if seen(mass(c))}
             for n in range(1, seed_len + 1)}
    for n in comps:
        res.history.append(("MASFRG", n, len(comps[n])))
    comps = prune(comps, 2)
    seqs = cq(comps)
    for n in seqs:
        res.history.append(("CQ", n, len(seqs[n])))

    def score(seq):
        """(observed, in-window) over the suffixes of `seq` up to `support_len` residues."""
        o = t = 0
        m = water
        for k in range(1, min(len(seq), support_len) + 1):
            m += table[seq[-k]]
            if lo <= m <= hi:
                t += 1
                o += observed.contains(m)
        return o, t

    def accept(seq):
        o, t = score(seq)
        return t == 0 or o >= support * t

    cur = seqs[seed_len]
    # Seeds whose own mass is in the window but unobserved carry no evidence: drop them.
    cur = {s for s in cur if seen(mass(s))}
    best = {}
    n = seed_len + 1
    while cur and n <= max_protein:
        oriented = cur | {s[::-1] for s in cur}
        by_head = {}
        for t in oriented:
            by_head.setdefault(t[:-1], []).append(t)
        kept, used = {}, set()
        for s in oriented:
            for t in by_head.get(s[1:], ()):
                new = s + t[-1]
                if accept(new):
                    c = canon(new)
                    o, _ = score(new)
                    # carry the support accumulated along the chain
                    kept[c] = max(kept.get(c, 0), best.get(canon(s), 0) + o)
                    used.add(canon(s))
                    used.add(canon(t))
        for s in cur - used:
            if len(s) >= min_contig:
                res.contigs.append((s, best.get(s, 0)))
        if len(kept) > beam:
            res.converged = False
            res.reason = f"beam of {beam} reached at length {n}"
            kept = dict(sorted(kept.items(), key=lambda kv: -kv[1])[:beam])
        best = kept
        cur = set(kept)
        res.history.append(("OVERLAP", n, len(cur)))
        n += 1
    res.contigs.sort(key=lambda cs: (-len(cs[0]), -cs[1]))
    return res
