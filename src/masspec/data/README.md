# Packaged data

What `masspec.example.pxd071026()` reads, so the notebook's real-data section
runs without the raw file.

| File | What | Source | Licence |
|---|---|---|---|
| `pxd071026_Mn_0s_1.tsv` | the run's MS1 scans merged into 5440 neutral masses (mass, summed intensity, distinct scans, charges seen); uncalibrated | run `Mn_0s_1` of PRIDE [PXD071026](https://www.ebi.ac.uk/pride/archive/projects/PXD071026), via `features.masses_from_mzml(ppm=10, min_scans=2)`; regenerate with `scripts/package_pxd071026.py` | CC0, as the deposit |
| `myoglobin_6Z4T.fasta` | sperm whale myoglobin, the protein in the run | RCSB PDB [6Z4T](https://www.rcsb.org/structure/6Z4T), as deposited with PXD071026 | CC0 (PDB data) |
| `pepsin_A_P00791.fasta` | mature porcine pepsin A (residues 63-388), the digesting enzyme, used only to calibrate | UniProt [P00791](https://www.uniprot.org/uniprotkb/P00791) | CC BY 4.0, (c) the UniProt Consortium |

The study behind the run: Gao H, Camargo EA, Ubi JN, Duan X, Tian X, Deng H,
Zheng G, Gao S. Hydrogen-deuterium exchange reveals catalytically linked
protein flexibility in myoglobin-mediated intramolecular C(sp3)-H activation.
*Protein Science* 35(1):e70410 (2026),
[doi:10.1002/pro.70410](https://doi.org/10.1002/pro.70410).
