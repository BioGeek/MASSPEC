"""Real runs, packaged with the code: PRIDE PXD071026's undeuterated pepsin
digest of sperm whale myoglobin (run Mn_0s_1, LTQ Orbitrap), and the de novo
reads of two partial acid hydrolysates of horse myoglobin (PXD003804,
PXD030203; `acid_hydrolysates()`).

The raw file is 148 MB and needs ThermoRawFileParser; what the method uses
from it is one list of neutral masses, 5440 of them, which ships here
(scripts/package_pxd071026.py writes it from the run). So the real-data half
of the notebook runs anywhere the package is installed, Colab included.
Provenance and licences are in data/README.md beside the files.
"""
from __future__ import annotations

from importlib.resources import files
from typing import NamedTuple

from .masses import normalise

_DATA = files("masspec") / "data"


class Prediction(NamedTuple):
    """One MS/MS spectrum as InstaNovo read it."""

    scan: int
    precursor_mz: float
    charge: int
    peptide: str
    confidence: float
    ppm_error: float


class Run(NamedTuple):
    """The packaged run: uncalibrated masses, and the two sequences."""

    masses: list
    myoglobin: str
    pepsin: str


def _sequence(name):
    """A packaged FASTA's sequence, with I written as L."""
    lines = (_DATA / name).read_text().splitlines()
    return normalise("".join(line.strip() for line in lines
                             if not line.startswith(">")))


def pxd071026():
    """The run's merged MS1 masses, myoglobin's sequence and pepsin's.

    The masses are as measured: about 5 ppm low, which `calibrate()` finds
    from pepsin's autolysis fragments.
    """
    rows = (_DATA / "pxd071026_Mn_0s_1.tsv").read_text().splitlines()[1:]
    masses = [float(row.split("\t")[0]) for row in rows]
    return Run(masses, _sequence("myoglobin_6Z4T.fasta"),
               _sequence("pepsin_A_P00791.fasta"))


def pxd071026_instanovo():
    """InstaNovo's reading of every MS/MS spectrum in the same run.

    InstaNovo 1.2.2 with the instanovo-v1.2.0 model and InstaNovo+
    refinement, run on a GPU (scripts/package_instanovo_pxd071026.py has
    the command). `peptide` is ProForma with UNIMOD modifications;
    `confidence` is the model's probability for its prediction; `ppm_error`
    the gap between the prediction's mass and the precursor's.
    """
    rows = (_DATA / "pxd071026_Mn_0s_1_instanovo.tsv").read_text()
    out = []
    for row in rows.splitlines()[1:]:
        scan, mz, charge, peptide, confidence, ppm_error = row.split("\t")
        out.append(Prediction(int(scan), float(mz), int(charge), peptide,
                              float(confidence), float(ppm_error)))
    return out


class Read(NamedTuple):
    """One distinct peptide as PEAKS read it de novo from MS/MS."""

    peptide: str
    mass: float
    charge: int
    alc: float


class Hydrolysate(NamedTuple):
    """A partial acid hydrolysate: its de novo reads and its protein."""

    accession: str
    reads: list
    myoglobin: str


_HYDROLYSATES = {
    "PXD003804": "pxd003804_myo_rep1.tsv",
    "PXD030203": "pxd030203_apo_mb_maah_1.tsv",
}


def acid_hydrolysates():
    """Two microwave HCl hydrolysates of horse myoglobin, from PRIDE.

    PXD003804 (Q Exactive) and PXD030203 (Orbitrap Fusion Lumos), one run
    each. Each read is a distinct peptide PEAKS sequenced de novo with an
    average local confidence (`alc`) of at least 50%, written as PEAKS
    writes it (modifications in brackets, I and L as read), with the
    MEASURED neutral mass of its most confident spectrum's precursor.
    scripts/package_acid_hydrolysis.py writes them from the deposits.
    """
    myoglobin = _sequence("myoglobin_horse_P68082.fasta")
    out = []
    for accession, name in _HYDROLYSATES.items():
        reads = []
        for row in (_DATA / name).read_text().splitlines()[1:]:
            peptide, mass, charge, alc = row.split("\t")
            reads.append(Read(peptide, float(mass), int(charge), float(alc)))
        out.append(Hydrolysate(accession, reads, myoglobin))
    return out
