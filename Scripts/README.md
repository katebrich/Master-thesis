# Pipeline for statistical analysis of ligand-binding sites properties
This folder contains scripts for analysis of ligand-binding sites properties (called features). User can define a new feature and implement a method for its computation. The pipeline is able to download PDB and FASTA files from PDBe, label binding residues based on the ligands in the PDB file, compute values of defined features, perform statistical analysis to find out whether these features are interesting for the recognition of binding sites, and draw various histograms and plots. Furthermore, the extended pipeline can use obtained PDB files and features to train a new model in P2Rank predictor of binding sites (http://siret.ms.mff.cuni.cz/p2rank) to see how the new features help the prediction in reality

The only needed input is a dataset file with proteins.

## Prerequisities:
- xxx
- zzz
