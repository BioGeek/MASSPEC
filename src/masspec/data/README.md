# Packaged data

What `masspec.example.pxd071026()`, `pxd071026_instanovo()` and
`acid_hydrolysates()` read, so the notebook's real-data sections run without
the raw files or a GPU.

| File | What | Source | Licence |
|---|---|---|---|
| `pxd071026_Mn_0s_1.tsv` | the run's MS1 scans merged into 5440 neutral masses (mass, summed intensity, distinct scans, charges seen); uncalibrated | run `Mn_0s_1` of PRIDE [PXD071026](https://www.ebi.ac.uk/pride/archive/projects/PXD071026), via `features.masses_from_mzml(ppm=10, min_scans=2)`; regenerate with `scripts/package_pxd071026.py` | CC0, as the deposit |
| `pxd071026_Mn_0s_1_instanovo.tsv` | InstaNovo's prediction for each of the run's 3843 MS/MS spectra with a charge state: scan, precursor m/z and charge, peptide (ProForma), confidence, precursor mass error (ppm) | InstaNovo 1.2.2, model instanovo-v1.2.0 with InstaNovo+ instanovoplus-v1.1.0 refinement, on the run's mzML; the command is in `scripts/package_instanovo_pxd071026.py` | derived from CC0 data |
| `myoglobin_6Z4T.fasta` | sperm whale myoglobin, the protein in the run | RCSB PDB [6Z4T](https://www.rcsb.org/structure/6Z4T), as deposited with PXD071026 | CC0 (PDB data) |
| `pepsin_A_P00791.fasta` | mature porcine pepsin A (residues 63-388), the digesting enzyme, used only to calibrate | UniProt [P00791](https://www.uniprot.org/uniprotkb/P00791) | CC BY 4.0, (c) the UniProt Consortium |
| `pxd003804_myo_rep1.tsv` | 5449 distinct peptides PEAKS read *de novo* (ALC >= 50%) from a microwave HCl hydrolysate of horse myoglobin (Q Exactive): peptide as PEAKS wrote it, measured neutral precursor mass, charge, ALC | `Myo_Rep1_de-Novo_Peaks-output.csv` of PRIDE [PXD003804](https://www.ebi.ac.uk/pride/archive/projects/PXD003804); regenerate with `scripts/package_acid_hydrolysis.py` | [EMBL-EBI Terms of Use](https://www.ebi.ac.uk/about/terms-of-use/), as the deposit |
| `pxd030203_apo_mb_maah_1.tsv` | the same for 16462 peptides from a microwave HCl hydrolysate of horse apo-myoglobin (Orbitrap Fusion Lumos) | `apo_MB_MAAH_1_peaks.csv` of PRIDE [PXD030203](https://www.ebi.ac.uk/pride/archive/projects/PXD030203); same script | CC0, as the deposit |
| `myoglobin_horse_P68082.fasta` | horse myoglobin, mature chain (residues 2-154), the protein in both hydrolysates | UniProt [P68082](https://www.uniprot.org/uniprotkb/P68082) | CC BY 4.0, (c) the UniProt Consortium |

The study behind the run: Gao H, Camargo EA, Ubi JN, Duan X, Tian X, Deng H,
Zheng G, Gao S. Hydrogen-deuterium exchange reveals catalytically linked
protein flexibility in myoglobin-mediated intramolecular C(sp3)-H activation.
*Protein Science* 35(1):e70410 (2026),
[doi:10.1002/pro.70410](https://doi.org/10.1002/pro.70410).

The study behind PXD003804: Savidor A, Barzilay R, Elinger D, Yarden Y,
Lindzen M, Gabashvili A, Adiv Tal O, Levin Y. Database-independent Protein
Sequencing (DiPS) enables full-length de novo protein and antibody sequence
determination. *Mol Cell Proteomics* 16(6):1151-1161 (2017),
[doi:10.1074/mcp.O116.065417](https://doi.org/10.1074/mcp.O116.065417).

The study behind PXD030203: Mai Z-B, Zhou Z-H, He Q-Y, Zhang G. Highly
robust de novo full-length protein sequencing. *Anal Chem* 94(8):3467-3475
(2022), [doi:10.1021/acs.analchem.1c03718](https://doi.org/10.1021/acs.analchem.1c03718).
