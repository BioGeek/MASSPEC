"""The paper's insulin experiment, number for number where the paper gives
numbers."""
from masspec.fragments import ideal_masses
from masspec.masses import H2O_1970, MASSES_1970
from masspec.mspec import Observed, run

# Bovine insulin in the paper's alphabet: J = Leu/Ile, O = Tyr.
CHAIN_A = "GJVEQCCASVCSJOQJENOCN"
CHAIN_B = "FVNQHJCGSHJVEAJOJVCGERGFFOTPKA"


def _run(tol):
    obs = Observed(ideal_masses([CHAIN_A, CHAIN_B], 9, 1000.0, MASSES_1970,
                                H2O_1970),
                   tol_da=tol)
    return run(obs, MASSES_1970, H2O_1970, mass_limit=1000.0)


def _count(res, stage, length):
    return next(count for step, size, count in res.history
                if step == stage and size == length)


def _chain_order(chain):
    """Longest first, then alphabetical."""
    return (-len(chain), chain)


def test_masfrg_counts_match_fig6():
    result = _run(0.02)
    counts = (_count(result, "MASFRG", 1), _count(result, "MASFRG", 2),
              _count(result, "MASFRG", 3))
    assert counts == (16, 66, 197)


def test_prune_matches_fig7():
    assert _count(_run(0.02), "PRUNE", 3) == 93


def test_cq_matches_fig8():
    assert _count(_run(0.02), "CQ", 3) == 172


def test_ten_mers_match_text():
    # "Thirty-three fragments 10 long remain possible, which is the correct
    # number."
    assert _count(_run(0.02), "OVERLAP", 10) == 33


def test_insulin_solved_at_0_02():
    chains = _run(0.02).chains(15)
    assert chains == sorted([CHAIN_B[::-1], CHAIN_A], key=_chain_order)


def test_insulin_not_solved_at_0_05():
    # "the insulin problem could be solved with a resolution of +-0.02
    # a.m.u., but could not be solved at +-0.05": at 0.05 Q and K (0.036
    # apart) can no longer be told apart.
    chains = _run(0.05).chains(20)
    assert CHAIN_A in chains and len(
        [chain for chain in chains if len(chain) == 21]) > 1
