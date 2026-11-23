"""MASSPEC taken to present-day data.

The 1970 method keeps a candidate only when EVERY fragment inside it was
observed, and treats masses above the instrument's limit as observed. Both
assumptions fail on an LC-MS run of a non-specific digest: fragments below
~350 Da are hardly ever seen, and many inside the window are missed.

This keeps the paper's phases and data structures, and changes three things,
each a direct generalisation of something in the paper:

1. **An observable window at both ends.** The paper assumed every mass above
   1000 a.m.u. observed, because its instrument could not measure it. Here
   masses outside ``[lo, hi]`` are assumed observed for the same reason:
   below ``lo`` the fragment cannot be seen, so its absence says nothing.
   That lets the bottom-up phases (MASFRG, PRUNE, CQ) start without the short
   fragments.
2. **Support instead of completeness.** A candidate of length N is judged on
   the fragments that END at its last residue (its suffixes): exactly the
   fragments the extension introduces. The paper requires all of them
   observed; here it is enough that a fraction ``support`` of those with a
   mass in the window, and at most ``support_len`` residues long, were
   observed. Long fragments are rare in a real digest, so counting them would
   sink true sequences.
3. **A beam.** With support below 1, spurious candidates survive more easily,
   and the lists can grow without bound. The paper hit the same wall ("the
   storage is exceeded and the problem must be terminated"). Here the
   ``beam`` best-supported candidates of each length are kept.

Masses are matched within a ppm tolerance. Modifications enter as extra
residue letters (see `masses.MODIFIED`). The result is a set of contigs: with
incomplete data a chain breaks where the digest left no overlapping evidence,
so the protein comes back in pieces whose order masses alone cannot give.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations_with_replacement

from .masses import H2O, MASSES
from .mspec import Observed, canon, cq, prune


@dataclass
class ModernResult:
    # (sequence, support score), longest first
    contigs: list = field(default_factory=list)
    history: list = field(default_factory=list)
    converged: bool = True
    reason: str = ""


class FastObserved:
    """Observed masses binned for O(1) lookup within a ppm tolerance."""

    def __init__(self, masses, tol_ppm):
        self.tol_ppm = tol_ppm
        self.masses = sorted(masses)
        if self.masses:
            self.width = max(self.masses) * tol_ppm * 1e-6
        else:
            self.width = 1.0
        self.bins = {}
        for mass in self.masses:
            self.bins.setdefault(int(mass / self.width), []).append(mass)

    def contains(self, mass):
        return self.match(mass) is not None

    def match(self, mass):
        """The observed mass matching `mass` within the tolerance.

        The nearest one, or None.
        """
        tolerance = mass * self.tol_ppm * 1e-6
        centre = int(mass / self.width)
        best = None
        for bin_index in (centre - 1, centre, centre + 1):
            for candidate in self.bins.get(bin_index, ()):
                distance = abs(candidate - mass)
                if distance <= tolerance and (
                        best is None or distance < abs(best - mass)):
                    best = candidate
        return best

    def __len__(self):
        return len(self.masses)


def _by_score_descending(item):
    """Sort key for (chain, (score, used)) items: highest score first.

    Ties go to the alphabetically first chain, so the beam does not depend
    on set iteration order, which changes with Python's hash seed.
    """
    return (-item[1][0], item[0])


def _by_value_descending(item):
    """Sort key for (key, score) items: highest score, then key, first."""
    return (-item[1], item[0])


def _by_score_then_length(contig):
    """Sort key for (sequence, score): highest score, then longest, first.

    Remaining ties go to the alphabetically first sequence.
    """
    return (-contig[1], -len(contig[0]), contig[0])


def _by_length_then_score(contig):
    """Sort key for (sequence, score): longest, then highest score, first.

    Remaining ties go to the alphabetically first sequence.
    """
    return (-len(contig[0]), -contig[1], contig[0])


def run_extend(observed, table=MASSES, water=H2O, lo=350.0, hi=6000.0,
               seed_len=4, support=0.5, support_len=20, min_hits=1,
               beam=2000, min_contig=6, max_protein=3000) -> ModernResult:
    """Grow chains one residue at a time, not by overlapping windows.

    The 1970 method overlaps fixed-length windows. OVERLAP builds a length-N
    candidate only from two surviving length N-1 windows, so one missing
    window anywhere breaks every chain through it. Here a chain is extended
    by any residue whose suffixes (the fragments ending at the new residue)
    are supported: at least `min_hits` observed, and at least a fraction
    `support` of those in the window. Both directions are grown, since a seed
    may sit anywhere in the protein. The seeds are MASSPEC's own: CQ's
    sequences of `seed_len` residues whose mass, when in the window, was
    observed.
    """
    res = ModernResult()
    letters = sorted(table)

    def mass_of(seq):
        return sum(table[residue] for residue in seq) + water

    def seen(mass):
        return mass < lo or mass > hi or observed.contains(mass)

    comps = {
        length: {comp
                 for comp in combinations_with_replacement(letters, length)
                 if seen(mass_of(comp))}
        for length in range(1, seed_len + 1)
    }
    seqs = cq(prune(comps, 2))
    seeds = {seq for seq in seqs[seed_len] if seen(mass_of(seq))}
    cur = {}
    # Sorted: iterating the set directly would make the run depend on
    # Python's hash seed, which differs between processes.
    for seed in sorted(seeds):
        cur[seed] = 0
        cur[seed[::-1]] = 0
    res.history.append(("SEEDS", seed_len, len(seeds)))

    def suffix_support(seq, used):
        """Observed masses (not yet used by this chain) among the suffixes.

        Suffixes up to `support_len` residues, and the number of suffixes in
        the window.
        """
        hits, in_window = [], 0
        mass = water
        for depth in range(1, min(len(seq), support_len) + 1):
            mass += table[seq[-depth]]
            if lo <= mass <= hi:
                in_window += 1
                matched = observed.match(mass)
                if matched is not None and matched not in used:
                    hits.append(matched)
        return hits, in_window

    length = seed_len + 1
    # chain -> (score, masses it has used). A mass is ONE fragment of the
    # protein, so it may support a chain once: otherwise a low-complexity
    # chain (ALASALSAL...) whose suffix masses fall on the same observed
    # masses at every turn grows forever on a handful of them.
    state = {seq: (score, frozenset()) for seq, score in cur.items()}
    while state and length <= max_protein:
        nxt = {}
        for seq, (score, used) in state.items():
            grown = False
            for letter in letters:
                new = seq + letter
                hits, in_window = suffix_support(new, used)
                n_observed = len(set(hits))
                if (n_observed >= min_hits
                        and n_observed >= support * in_window):
                    if new not in nxt or nxt[new][0] < score + n_observed:
                        nxt[new] = (score + n_observed, used | set(hits))
                    grown = True
            if not grown and len(seq) >= min_contig:
                res.contigs.append((canon(seq), score))
        if len(nxt) > beam:
            res.converged = False
            res.reason = f"beam of {beam} reached"
            nxt = dict(sorted(nxt.items(), key=_by_score_descending)[:beam])
        res.history.append(("EXTEND", length, len(nxt)))
        state = nxt
        length += 1
    best = {}
    for seq, score in res.contigs:
        best[seq] = max(best.get(seq, -1), score)
    res.contigs = sorted(best.items(), key=_by_score_then_length)
    return res


def run_modern(observed: Observed, table=MASSES, water=H2O, lo=350.0,
               hi=6000.0, seed_len=4, support=0.6, support_len=12,
               beam=5000, min_contig=6, max_protein=3000) -> ModernResult:
    res = ModernResult()
    letters = sorted(table)

    def mass_of(seq):
        return sum(table[residue] for residue in seq) + water

    def seen(mass):
        return mass < lo or mass > hi or observed.contains(mass)

    # MASFRG with the window: compositions whose mass is observed, or cannot
    # be.
    comps = {
        length: {comp
                 for comp in combinations_with_replacement(letters, length)
                 if seen(mass_of(comp))}
        for length in range(1, seed_len + 1)
    }
    for length in comps:
        res.history.append(("MASFRG", length, len(comps[length])))
    comps = prune(comps, 2)
    seqs = cq(comps)
    for length in seqs:
        res.history.append(("CQ", length, len(seqs[length])))

    def score(seq):
        """(observed, in-window) over the suffixes of `seq`.

        Suffixes up to `support_len` residues.
        """
        n_observed = in_window = 0
        mass = water
        for depth in range(1, min(len(seq), support_len) + 1):
            mass += table[seq[-depth]]
            if lo <= mass <= hi:
                in_window += 1
                n_observed += observed.contains(mass)
        return n_observed, in_window

    def accept(seq):
        n_observed, in_window = score(seq)
        return in_window == 0 or n_observed >= support * in_window

    cur = seqs[seed_len]
    # Seeds whose own mass is in the window but unobserved carry no evidence:
    # drop them.
    cur = {seq for seq in cur if seen(mass_of(seq))}
    best = {}
    length = seed_len + 1
    while cur and length <= max_protein:
        oriented = cur | {seq[::-1] for seq in cur}
        by_head = {}
        for tail in oriented:
            by_head.setdefault(tail[:-1], []).append(tail)
        kept, used = {}, set()
        for head in oriented:
            for tail in by_head.get(head[1:], ()):
                new = head + tail[-1]
                if accept(new):
                    key = canon(new)
                    n_observed, _ = score(new)
                    # carry the support accumulated along the chain
                    kept[key] = max(kept.get(key, 0),
                                    best.get(canon(head), 0) + n_observed)
                    used.add(canon(head))
                    used.add(canon(tail))
        for seq in sorted(cur - used):
            if len(seq) >= min_contig:
                res.contigs.append((seq, best.get(seq, 0)))
        if len(kept) > beam:
            res.converged = False
            res.reason = f"beam of {beam} reached at length {length}"
            kept = dict(sorted(kept.items(),
                               key=_by_value_descending)[:beam])
        best = kept
        cur = set(kept)
        res.history.append(("OVERLAP", length, len(cur)))
        length += 1
    res.contigs.sort(key=_by_length_then_score)
    return res
