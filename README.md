# Computational Analysis of TP53 Mutations in Colorectal Cancer

## Overview
Analysed somatic TP53 missense mutations and their predicted impact on
protein stability in colorectal cancer using sequence and structural
bioinformatics approaches.

## Tools & Technologies
- **Databases:** NCBI ClinVar, NCBI Protein, UniProt
- **Sequence Analysis:** BLASTp, Biopython (Entrez, SeqIO)
- **Structural Visualisation:** PyMOL (PDB: 2OCJ)
- **Data Analysis & Visualisation:** R, ggplot2

## What I Did
- Retrieved 40+ TP53 missense mutation records from NCBI ClinVar & UniProt
- Performed BLASTp alignment of wild-type vs mutant TP53 sequences
- Wrote Python (Biopython) scripts to automate FASTA retrieval and
  BLAST XML parsing
- Used R to analyse physicochemical properties of mutant residues
- Visualised TP53 DNA-binding domain structures in PyMOL to identify
  structural perturbations at mutation sites
- Identified 3 high-confidence loss-of-function variants:
  **R175H, R248W, R273H**

## Key Finding
These hotspot mutations cause significant predicted destabilisation of
the TP53 tetramerisation domain — consistent with published literature
and relevant to TP53-reactivation drug strategies.

## Skills Demonstrated
`Python` `Biopython` `BLASTp` `R` `PyMOL` `NCBI` `ClinVar` `Data Analysis`
