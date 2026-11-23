"""From a present-day LC-MS run to MASSPEC's input: one list of peptide
masses.

MASSPEC expects what Dayhoff and Eck imagined, "a single spectrum containing
all the fragment masses". An LC-MS run spreads the fragments over thousands
of MS1 scans, each peptide over several charge states and isotope peaks. This
module collapses that back to one list of neutral monoisotopic masses:

1. In every MS1 scan, find isotope envelopes: a peak followed by one or two
   peaks 1.00336/z apart (within `ppm`), for charges 1 to `max_charge`, with
   no comparably intense peak 1.00336/z BEFORE it (else it is not the
   monoisotopic peak). Each envelope gives a neutral mass
   (m/z - proton) * z.
2. Merge the masses of all scans within `ppm`, keeping those seen in at least
   `min_scans` scans, which removes most single-scan noise.

It is deliberately simple, so a reader can follow every step; a dedicated
feature finder (Dinosaur, OpenMS, MaxQuant) does better, and `read_mass_list`
accepts their output instead.

Needs the optional dependencies: ``uv sync --extra real`` (pyteomics, lxml,
numpy).
"""
from __future__ import annotations

import csv
from pathlib import Path

from .masses import ISOTOPE_SPACING, PROTON


def envelopes(mz, intensity, ppm=10.0, max_charge=6, min_rel=0.01,
              mono_ratio=0.6):
    """Neutral monoisotopic masses (with intensities) of the isotope
    envelopes in one centroided scan."""
    import numpy as np
    mz = np.asarray(mz, dtype=float)
    intensities = np.asarray(intensity, dtype=float)
    if mz.size == 0:
        return []
    order = np.argsort(mz)
    mz, intensities = mz[order], intensities[order]
    keep = intensities >= intensities.max() * min_rel
    mz, intensities = mz[keep], intensities[keep]

    def find(target):
        half_width = target * ppm * 1e-6
        first = np.searchsorted(mz, target - half_width)
        last = np.searchsorted(mz, target + half_width, side="right")
        if first == last:
            return -1
        return first + int(np.argmax(intensities[first:last]))

    out, used = [], set()
    for peak in np.argsort(-intensities):
        if peak in used:
            continue
        for charge in range(max_charge, 0, -1):
            step = ISOTOPE_SPACING / charge
            n1 = find(mz[peak] + step)
            if n1 < 0:
                continue
            n2 = find(mz[peak] + 2 * step)
            before = find(mz[peak] - step)
            # a comparably intense peak one isotope BELOW means `peak` is not
            # the monoisotopic peak
            if (before >= 0
                    and intensities[before] >= mono_ratio * intensities[peak]):
                continue
            members = [peak, n1] + ([n2] if n2 >= 0 else [])
            used.update(members)
            out.append(((mz[peak] - PROTON) * charge,
                        float(intensities[members].sum()), charge))
            break
    return out


def masses_from_mzml(path, ppm=10.0, max_charge=6, min_scans=2,
                     rt_range=None):
    """One merged list of neutral masses from every MS1 scan of a centroided
    mzML file.

    Returns [(mass, total_intensity, n_scans, charges)], sorted by mass.
    """
    from pyteomics import mzml
    raw = []
    with mzml.read(str(path)) as reader:
        for spec in reader:
            if spec.get("ms level") != 1:
                continue
            if rt_range:
                rt = spec["scanList"]["scan"][0].get("scan start time", 0.0)
                if not rt_range[0] <= float(rt) <= rt_range[1]:
                    continue
            for mass, inten, charge in envelopes(spec["m/z array"],
                                                 spec["intensity array"],
                                                 ppm, max_charge):
                raw.append((mass, inten, charge))
    return merge(raw, ppm, min_scans)


def merge(raw, ppm=10.0, min_scans=2):
    """Cluster per-scan masses within `ppm`: intensity-weighted mean, scan
    count, charges seen."""
    raw = sorted(raw)
    out, group = [], []

    def flush():
        if len(group) >= min_scans:
            weight = sum(member[1] for member in group)
            out.append((
                sum(member[0] * member[1] for member in group) / weight,
                weight,
                len(group),
                sorted({member[2] for member in group}),
            ))

    for record in raw:
        if group and record[0] - group[-1][0] > record[0] * ppm * 1e-6:
            flush()
            group = []
        group.append(record)
    if group:
        flush()
    return out


def calibrate(masses, calibrant, window_ppm=20.0, max_len=40):
    """Systematic mass error (ppm) estimated from a KNOWN protein in the same
    run.

    Like a lock mass, but from peptides: in a pepsin digest, pepsin digests
    itself, so its autolysis fragments are in every run and its sequence is
    public. Every mass is compared with every fragment of `calibrant` (up to
    `max_len` residues); true matches pile up at the instrument's offset while
    chance matches spread evenly, so the offset is the mode of the errors
    (2 ppm bins), refined to the median of the errors within 2 ppm of it. It
    never looks at the protein being sequenced. Returns (offset_ppm,
    matches_in_peak, matches_per_background_bin).
    """
    import bisect
    import statistics

    from .masses import H2O, MASSES
    seq = calibrant.upper().replace("I", "L")
    frags = set()
    for start in range(len(seq)):
        mass = H2O
        for end in range(start, min(len(seq), start + max_len)):
            if seq[end] not in MASSES:
                break
            mass += MASSES[seq[end]]
            frags.add(mass)
    frags = sorted(frags)
    errs = []
    for observed in masses:
        insert_at = bisect.bisect_left(frags, observed)
        for fragment in frags[max(0, insert_at - 1):insert_at + 1]:
            error = (observed - fragment) / fragment * 1e6
            if abs(error) <= window_ppm:
                errs.append(error)
    if not errs:
        return 0.0, 0, 0.0
    bins = {}
    for error in errs:
        bins.setdefault(int((error + window_ppm) // 2), []).append(error)

    def bin_size(key):
        return len(bins[key])

    peak = max(bins, key=bin_size)
    mode = statistics.median(bins[peak])
    near = [error for error in errs if abs(error - mode) <= 2]
    background = (len(errs) - len(bins[peak])) / max(1, len(bins) - 1)
    return statistics.median(near), len(bins[peak]), background


def read_mass_list(path, column=None):
    """Neutral masses from a table another tool wrote.

    Understands a plain list (one mass per line), a CSV/TSV with a mass column
    (`column`, or the first of: mass, Mass, neutral_mass, monoisotopic_mass,
    "Monoisotopic Mass"), MaxQuant's allPeptides.txt (column "Mass"), and
    Dinosaur's .features.tsv (column "mass").
    """
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace").splitlines()
    first = text[0] if text else ""
    if not any(char in first for char in ",\t;"):
        return sorted(float(line) for line in text
                      if line.strip() and not line.startswith("#"))
    delim = "\t" if "\t" in first else ("," if "," in first else ";")
    rows = csv.DictReader(text, delimiter=delim)
    names = rows.fieldnames or []
    candidates = ("mass", "Mass", "neutral_mass", "monoisotopic_mass",
                  "Monoisotopic Mass", "MonoisotopicMass")
    col = column or next(
        (candidate for candidate in candidates if candidate in names), None)
    if col is None:
        raise SystemExit(
            f"{path}: no mass column among {names}; pass --column")
    return sorted(float(row[col]) for row in rows
                  if row.get(col) not in (None, "", "NaN"))
