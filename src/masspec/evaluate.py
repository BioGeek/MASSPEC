"""Choosing an answer from the contigs, and scoring it against a known sequence.

A mass list cannot fix every residue. Some swaps are invisible to mass at any accuracy: Q is exactly
G+A, and N is exactly G+G (the same elemental formulas), and a run of residues can be reordered
when no fragment ends inside it. The paper says as much: "many different groups of amino acids ...
have the same chemical composition, for example, Ala-Val and Gly-Leu." So a reconstruction is scored
residue by residue after aligning it to the truth, and a mismatch that preserves mass over a short
stretch is counted apart from a plain error.
"""
from __future__ import annotations

from .masses import MASSES


def _similarity(a, b):
    """Gapped similarity of two contigs in either orientation, relative to the shorter one: the
    fraction of its residues in matching blocks. An indel (Q for GA shifts everything after it)
    does not hide a near-variant, as it would from an ungapped comparison."""
    from difflib import SequenceMatcher
    short = min(len(a), len(b))
    return max(sum(m.size for m in SequenceMatcher(None, x, b, autojunk=False).get_matching_blocks()) / short
               for x in (a, a[::-1]))


def assemble(contigs, redundant=0.6, min_score_frac=0.25, candidates=400, top=None):
    """Greedy answer: the best-supported contigs, skipping variants of ones already chosen.

    `contigs` is [(sequence, score)]. Taken best score first; a contig is skipped when it overlaps a
    chosen one with gapped similarity `redundant` or more (either orientation), which is what a
    near-variant (two residues swapped, Q for GA) looks like. Stops below `min_score_frac` of the best
    score: those contigs rest on too few masses to be worth reporting.
    """
    ranked = sorted(contigs, key=lambda cs: (-cs[1], -len(cs[0])))[:candidates]
    if not ranked:
        return []
    floor = ranked[0][1] * min_score_frac
    chosen = []
    for seq, score in ranked:
        if score < floor:
            break
        if any(_similarity(seq, c) >= redundant for c, _ in chosen):
            continue
        chosen.append((seq, score))
        if top and len(chosen) >= top:
            break
    return chosen


def _best_ungapped(contig, truth):
    """(identity, offset, orientation) of the best ungapped placement of `contig` on `truth`."""
    best = (-1.0, 0, contig)
    for oriented in (contig, contig[::-1]):
        n = len(oriented)
        for off in range(-n + 1, len(truth)):
            pairs = [(oriented[i], truth[off + i]) for i in range(n) if 0 <= off + i < len(truth)]
            if len(pairs) < n * 0.8:
                continue
            ident = sum(a == b for a, b in pairs) / n
            if ident > best[0]:
                best = (ident, off, oriented)
    return best


def score(contigs, truth, table=MASSES, min_identity=0.7):
    """Residue-level recall and precision of a set of contigs against the true sequence.

    Each contig is placed where it matches best. Its residues that equal the truth are CORRECT; a
    mismatching stretch of up to three residues whose total mass equals the truth's (0.02 Da) is
    MASS-EQUIVALENT (a reordering, or Q for GA); the rest are WRONG. A contig matching below
    `min_identity` everywhere counts as wrong throughout.
    """
    covered, eq_covered = set(), set()
    correct = equivalent = wrong = 0
    for seq, _ in contigs:
        ident, off, oriented = _best_ungapped(seq, truth)
        if ident < min_identity:
            wrong += len(seq)
            continue
        i = 0
        while i < len(oriented):
            p = off + i
            if 0 <= p < len(truth) and oriented[i] == truth[p]:
                correct += 1
                covered.add(p)
                i += 1
                continue
            # try to close a short mismatch by mass
            for w in (1, 2, 3):
                a = oriented[i:i + w]
                b = truth[max(p, 0):p + w] if p >= 0 else ""
                if len(a) == w and len(b) == w and abs(sum(table[x] for x in a) - sum(table[x] for x in b)) < 0.02:
                    equivalent += w
                    eq_covered.update(range(p, p + w))
                    i += w
                    break
            else:
                wrong += 1
                i += 1
    total = correct + equivalent + wrong
    return {
        "recall_exact": len(covered) / len(truth),
        "recall_with_equivalents": len(covered | eq_covered) / len(truth),
        "precision_exact": correct / total if total else 0.0,
        "precision_with_equivalents": (correct + equivalent) / total if total else 0.0,
        "contigs": len(contigs),
        "residues_output": total,
    }
