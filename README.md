# disulfide_bonds
data and scripts for the analysis of disulfide bonds in disulfide-rich peptides 

Files QueryResults_NMR.xls and QueryResults_XRD.xls 
These spreadsheets contain the datasets obtained from the search of the Protein Database. The datasets were produced using the ‘Advanced Search’ option in the PDB database (http://www.rcsb.org/, accessed April 2018) to find peptides with molecular mass <5.5 kDa that contain 1–6 disulfide bonds. The 5.5 kDa cut-off corresponds to molecules with ~50 residues, a commonly used definition to distinguish between peptides and proteins. Two separate queries for structures determined using NMR spectroscopy and X-ray crystallography (XRD) were carried out. For the latter, only structures with a resolution <2.5 Å were considered.

The spreadsheets contain the raw dataset as well as the 'clean' dataset in which  duplicate entries, de novo designed peptides, structures that are not peptides but domains of larger proteins, peptides with unnatural amino acids or chemical modification (e.g. glycosylation, fucosylation), structures of the same peptides solved at different temperatures, peptides with point-mutations of naturally occurring DRPs (e.g. the PDB ID 1RYW is a R29A mutant of the naturally occurring peptide hainantoxin-IV with PDB ID 1NIY) were removed. In addition there is a datasets for peptides that have disulfide bonds with CA-CA distances < 4.3 A (0.43 nm)

Files calc_ssbonds.py and calc_ssXtorsions.py
These files contain python scripts to calculate the geometries of disulfide bonds in PDB files. see comments in the scripts for details.
calc_ssbonds calculates the CA-CA distance in the disulfide bond 
calc_ssXtorsions calculates the five torsion angles in the disulfide bond 
