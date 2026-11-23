import pytest

np = pytest.importorskip("numpy")
from masspec.features import C13, PROTON, envelopes, merge  # noqa: E402


def test_envelope_finds_charge_and_monoisotopic_mass():
    mass = 1234.5678
    z = 2
    mono = mass / z + PROTON
    mz = [mono, mono + C13 / z, mono + 2 * C13 / z, 400.0]
    it = [100.0, 70.0, 30.0, 5.0]
    found = envelopes(mz, it, ppm=5)
    assert any(abs(m - mass) < 1e-3 and c == 2 for m, _, c in found)


def test_merge_requires_two_scans():
    raw = [(1000.0, 1.0, 2), (1000.001, 1.0, 2), (1500.0, 1.0, 1)]
    out = merge(raw, ppm=5, min_scans=2)
    assert len(out) == 1 and abs(out[0][0] - 1000.0005) < 1e-3
