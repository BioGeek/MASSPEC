"""Choosing an answer from the contigs, and scoring it against a known
sequence.

A mass list cannot fix every residue. Some swaps are invisible to mass at any
accuracy: Q is exactly G+A, and N is exactly G+G (the same elemental
formulas), and a run of residues can be reordered when no fragment ends
inside it. The paper says as much: "many different groups of amino acids ...
have the same chemical composition, for example, Ala-Val and Gly-Leu." So a
reconstruction is scored residue by residue after aligning it to the truth,
and a mismatch that preserves mass over a short stretch is counted apart from
a plain error.
"""
from __future__ import annotations

from .masses import MASSES


def _similarity(first, second):
    """Gapped similarity of two contigs in either orientation, relative to the
    shorter one: the fraction of its residues in matching blocks. An indel (Q
    for GA shifts everything after it) does not hide a near-variant, as it
    would from an ungapped comparison."""
    from difflib import SequenceMatcher
    short = min(len(first), len(second))
    return max(
        sum(block.size
            for block in SequenceMatcher(
                None, oriented, second, autojunk=False
            ).get_matching_blocks()) / short
        for oriented in (first, first[::-1])
    )


def _rank_key(entry):
    """Best score first, and among equal scores the longest contig first."""
    return (-entry[1], -len(entry[0]))


def assemble(contigs, redundant=0.6, min_score_frac=0.25, candidates=400,
             top=None):
    """Greedy answer: the best-supported contigs, skipping variants of ones
    already chosen.

    `contigs` is [(sequence, score)]. Taken best score first; a contig is
    skipped when it overlaps a chosen one with gapped similarity `redundant`
    or more (either orientation), which is what a near-variant (two residues
    swapped, Q for GA) looks like. Stops below `min_score_frac` of the best
    score: those contigs rest on too few masses to be worth reporting.
    """
    ranked = sorted(contigs, key=_rank_key)[:candidates]
    if not ranked:
        return []
    floor = ranked[0][1] * min_score_frac
    chosen = []
    for seq, score in ranked:
        if score < floor:
            break
        if any(_similarity(seq, chosen_seq) >= redundant
               for chosen_seq, _ in chosen):
            continue
        chosen.append((seq, score))
        if top and len(chosen) >= top:
            break
    return chosen


def _best_ungapped(contig, truth):
    """(identity, offset, orientation) of the best ungapped placement of
    `contig` on `truth`."""
    best = (-1.0, 0, contig)
    for oriented in (contig, contig[::-1]):
        length = len(oriented)
        for off in range(-length + 1, len(truth)):
            pairs = [(oriented[index], truth[off + index])
                     for index in range(length)
                     if 0 <= off + index < len(truth)]
            if len(pairs) < length * 0.8:
                continue
            ident = sum(found == expected
                        for found, expected in pairs) / length
            if ident > best[0]:
                best = (ident, off, oriented)
    return best


def score(contigs, truth, table=MASSES, min_identity=0.7):
    """Residue-level recall and precision of a set of contigs against the true
    sequence.

    Each contig is placed where it matches best. Its residues that equal the
    truth are CORRECT; a mismatching stretch of up to three residues whose
    total mass equals the truth's (0.02 Da) is MASS-EQUIVALENT (a reordering,
    or Q for GA); the rest are WRONG. A contig matching below `min_identity`
    everywhere counts as wrong throughout.
    """
    covered, eq_covered = set(), set()
    correct = equivalent = wrong = 0
    for seq, _ in contigs:
        ident, off, oriented = _best_ungapped(seq, truth)
        if ident < min_identity:
            wrong += len(seq)
            continue
        index = 0
        while index < len(oriented):
            position = off + index
            if (0 <= position < len(truth)
                    and oriented[index] == truth[position]):
                correct += 1
                covered.add(position)
                index += 1
                continue
            # try to close a short mismatch by mass
            for width in (1, 2, 3):
                found = oriented[index:index + width]
                expected = (truth[max(position, 0):position + width]
                            if position >= 0 else "")
                if (len(found) == width and len(expected) == width
                        and abs(sum(table[residue] for residue in found)
                                - sum(table[residue]
                                      for residue in expected)) < 0.02):
                    equivalent += width
                    eq_covered.update(range(position, position + width))
                    index += width
                    break
            else:
                wrong += 1
                index += 1
    total = correct + equivalent + wrong
    return {
        "recall_exact": len(covered) / len(truth),
        "recall_with_equivalents": len(covered | eq_covered) / len(truth),
        "precision_exact": correct / total if total else 0.0,
        "precision_with_equivalents": ((correct + equivalent) / total
                                       if total else 0.0),
        "contigs": len(contigs),
        "residues_output": total,
    }
