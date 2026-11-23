import pytest

np = pytest.importorskip("numpy")
from masspec.features import envelopes, merge  # noqa: E402
from masspec.masses import ISOTOPE_SPACING, PROTON  # noqa: E402


def test_envelope_finds_charge_and_monoisotopic_mass():
    mass = 1234.5678
    charge = 2
    mono = mass / charge + PROTON
    step = ISOTOPE_SPACING / charge
    mz = [mono, mono + step, mono + 2 * step, 400.0]
    intensities = [100.0, 70.0, 30.0, 5.0]
    found = envelopes(mz, intensities, ppm=5)
    assert any(abs(found_mass - mass) < 1e-3 and found_charge == 2
               for found_mass, _, found_charge in found)


def test_merge_requires_two_scans():
    raw = [(1000.0, 1.0, 2, 0), (1000.001, 1.0, 2, 1), (1500.0, 1.0, 1, 0)]
    out = merge(raw, ppm=5, min_scans=2)
    assert len(out) == 1 and abs(out[0][0] - 1000.0005) < 1e-3


def test_merge_counts_scans_not_charges():
    one_scan = [(1000.0, 1.0, 2, 7), (1000.0001, 1.0, 3, 7)]
    assert merge(one_scan, ppm=10, min_scans=2) == []


def test_merge_does_not_chain_past_the_tolerance():
    # seven masses 9 ppm apart: each is within 10 ppm of the next, the ends
    # 54 ppm apart; a 2 x 10 ppm window from each cluster's first mass makes
    # three clusters, not one
    chain = [(1000.0 + step * 0.009, 1.0, 2, step) for step in range(7)]
    assert len(merge(chain, ppm=10, min_scans=1)) == 3
