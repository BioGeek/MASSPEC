"""The paper's insulin experiment, number for number where the paper gives numbers."""
from masspec.fragments import ideal_masses
from masspec.masses import MASSES_1970, H2O_1970
from masspec.mspec import Observed, run

# Bovine insulin in the paper's alphabet: J = Leu/Ile, O = Tyr.
A = "GJVEQCCASVCSJOQJENOCN"
B = "FVNQHJCGSHJVEAJOJVCGERGFFOTPKA"


def _run(tol):
    obs = Observed(ideal_masses([A, B], 9, 1000.0, MASSES_1970, H2O_1970), tol_da=tol)
    return run(obs, MASSES_1970, H2O_1970, mass_limit=1000.0)


def _count(res, stage, n):
    return next(c for s, k, c in res.history if s == stage and k == n)


def test_masfrg_counts_match_fig6():
    r = _run(0.02)
    assert (_count(r, "MASFRG", 1), _count(r, "MASFRG", 2), _count(r, "MASFRG", 3)) == (16, 66, 197)


def test_prune_matches_fig7():
    assert _count(_run(0.02), "PRUNE", 3) == 93


def test_cq_matches_fig8():
    assert _count(_run(0.02), "CQ", 3) == 172


def test_ten_mers_match_text():
    # "Thirty-three fragments 10 long remain possible, which is the correct number."
    assert _count(_run(0.02), "OVERLAP", 10) == 33


def test_insulin_solved_at_0_02():
    chains = _run(0.02).chains(15)
    assert chains == sorted([B[::-1], A], key=lambda s: (-len(s), s))


def test_insulin_not_solved_at_0_05():
    # "the insulin problem could be solved with a resolution of +-0.02 a.m.u., but could not be
    # solved at +-0.05": at 0.05 Q and K (0.036 apart) can no longer be told apart.
    chains = _run(0.05).chains(20)
    assert A in chains and len([c for c in chains if len(c) == 21]) > 1
