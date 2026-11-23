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
    """Best score first, then the longest contig, then alphabetical.

    The last key makes the choice independent of the input order.
    """
    return (-entry[1], -len(entry[0]), entry[0])


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


_EXACT, _BLOCK, _MISMATCH, _OVERHANG = range(4)
_BLOCK_WIDTHS = [(found, expected) for found in (1, 2, 3)
                 for expected in (1, 2, 3)]


def _prefix_masses(seq, table):
    """Cumulative residue masses: element k is the mass of seq[:k]."""
    sums = [0.0]
    for residue in seq:
        sums.append(sums[-1] + table[residue])
    return sums


def _align(contig, truth, table, tolerance):
    """Best placement of `contig` on `truth`, allowing mass-equivalent blocks.

    The placement has no gaps. A contig residue either matches its truth
    residue (+1), mismatches it (-1), or belongs to a block of 1 to 3 contig
    residues standing for 1 to 3 different truth residues of the same total
    mass. That is how Q against GA, N against GG, or a reordering looks, and
    it is the only way the two sides can differ in length. A block earns 0.9
    per residue the two sides have in common and loses 1 per residue of
    length difference. A contig may run past either end of the truth (-1 per
    residue).

    General gaps are left out on purpose. With them, a short contig can be
    threaded through a 150-residue protein almost anywhere: scored against
    shuffled myoglobin, contigs from a real run reached 17% recall by chance,
    more than against myoglobin itself.

    Returns (score, steps), each step (kind, contig_start, contig_end,
    truth_start, truth_end).
    """
    rows, cols = len(contig) + 1, len(truth) + 1
    contig_mass = _prefix_masses(contig, table)
    truth_mass = _prefix_masses(truth, table)
    worst = float("-inf")
    best = [[worst] * cols for _ in range(rows)]
    came = [[None] * cols for _ in range(rows)]
    best[0] = [0.0] * cols
    for row in range(1, rows):
        for col in range(cols):
            options = []
            if col and contig[row - 1] == truth[col - 1]:
                options.append((best[row - 1][col - 1] + 1.0,
                                (_EXACT, 1, 1)))
            for width, span in _BLOCK_WIDTHS:
                if width > row or span > col:
                    continue
                found = contig[row - width:row]
                expected = truth[col - span:col]
                if found == expected:
                    continue
                gap = abs((contig_mass[row] - contig_mass[row - width])
                          - (truth_mass[col] - truth_mass[col - span]))
                if gap < tolerance:
                    value = 0.9 * min(width, span) - abs(width - span)
                    options.append((best[row - width][col - span] + value,
                                    (_BLOCK, width, span)))
            if col:
                options.append((best[row - 1][col - 1] - 1.0,
                                (_MISMATCH, 1, 1)))
            if col in (0, cols - 1):
                options.append((best[row - 1][col] - 1.0,
                                (_OVERHANG, 1, 0)))
            for value, step in options:
                if value > best[row][col]:
                    best[row][col], came[row][col] = value, step
    end = 0
    for col in range(1, cols):
        if best[-1][col] > best[-1][end]:
            end = col
    steps, row, col = [], rows - 1, end
    while row > 0:
        kind, width, span = came[row][col]
        steps.append((kind, row - width, row, col - span, col))
        row, col = row - width, col - span
    return best[-1][end], steps[::-1]


def score(contigs, truth, table=MASSES, min_identity=0.7, tolerance=0.02):
    """Residue-level recall and precision of a set of contigs against the true
    sequence.

    Each contig is placed where it fits best, in whichever orientation fits
    better (see `_align`). Its residues identical to the truth are CORRECT; a
    block of up to three residues standing for up to three different truth
    residues of the same total mass (within `tolerance` Da) is
    MASS-EQUIVALENT: a reordering, GA for Q, GG for N; the rest are WRONG. A
    contig with fewer than `min_identity` of its residues CORRECT counts as
    wrong throughout: mass-equivalent residues do not count towards it, as a
    block can be fitted to too many places.
    """
    covered, eq_covered = set(), set()
    correct = equivalent = wrong = 0
    for seq, _ in contigs:
        forward = _align(seq, truth, table, tolerance)
        backward = _align(seq[::-1], truth, table, tolerance)
        _, steps = forward if forward[0] >= backward[0] else backward
        exact = [step for step in steps if step[0] == _EXACT]
        blocks = [step for step in steps if step[0] == _BLOCK]
        n_equivalent = sum(step[2] - step[1] for step in blocks)
        if len(exact) / len(seq) < min_identity:
            wrong += len(seq)
            continue
        correct += len(exact)
        equivalent += n_equivalent
        wrong += len(seq) - len(exact) - n_equivalent
        covered.update(step[3] for step in exact)
        for step in blocks:
            eq_covered.update(range(step[3], step[4]))
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
