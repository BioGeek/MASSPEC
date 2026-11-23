from masspec.evaluate import score

TRUTH = "PEPTLDEQKLMNKVW"


def scored(contig):
    return score([(contig, 1.0)], TRUTH)


def test_exact_contig_in_either_orientation():
    for contig in ("PEPTLDEQKLMN", "NMLKQEDLTPEP"):
        result = scored(contig)
        assert result["precision_exact"] == 1.0
        assert result["recall_exact"] == 12 / 15


def test_ga_for_q_is_equivalent_and_keeps_the_rest_aligned():
    result = scored("PEPTLDEGAKLMN")
    assert result["precision_exact"] == 11 / 13
    assert result["precision_with_equivalents"] == 1.0
    assert result["recall_with_equivalents"] == 12 / 15


def test_reordering_is_equivalent():
    result = scored("PEPTLDEQLKMN")
    assert result["precision_exact"] == 10 / 12
    assert result["precision_with_equivalents"] == 1.0


def test_an_inserted_residue_is_an_error():
    # GG after the real N is not read as a second N: the placement has no
    # gaps, so an insertion that is not mass-equivalent stays an error
    result = scored("KLMNGGKVW")
    assert result["precision_with_equivalents"] < 1.0
    assert result["recall_exact"] <= 4 / 15


def test_shuffled_truth_scores_near_zero():
    from random import Random
    residues = list(TRUTH * 8)
    Random(0).shuffle(residues)
    decoy = "".join(residues)
    contigs = [(TRUTH[start:start + 8], 1.0) for start in range(0, 8, 2)]
    assert score(contigs, decoy)["recall_exact"] == 0.0


def test_unrelated_contig_is_wrong_throughout():
    result = scored("WWWWWWW")
    assert result["precision_with_equivalents"] == 0.0
    assert result["residues_output"] == 7
