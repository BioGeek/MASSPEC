"""A real run, packaged with the code: PRIDE PXD071026's undeuterated pepsin
digest of sperm whale myoglobin (run Mn_0s_1, LTQ Orbitrap).

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
