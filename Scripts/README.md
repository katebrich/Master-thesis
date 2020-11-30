# Pipeline for statistical analysis of ligand binding sites properties
This folder contains scripts for analysis of ligand binding sites properties (called features). The pipeline is able to download PDB and FASTA files from PDBe, label binding residues based on the ligands in the PDB file, compute values of defined features, perform statistical analysis to find out whether these features are interesting for the recognition of binding sites, and draw various histograms and plots. Furthermore, the extended pipeline can use obtained PDB files and features to train a new model in [P2Rank](http://siret.ms.mff.cuni.cz/p2rank) predictor of binding sites to see if the new features improve the prediction.

Users can define a feature of their own by implementing a method for its computation.

The only needed input is a dataset file with proteins.

Below, you can learn how to:
 - [download data from database and run statistical analysis](#one)
 - [run analysis on custom data obtained from a different source](#two)
 - [train a P2Rank model with new features](#three)


### Prerequisities:
* Python 3.x
* Installed Python packages:
  * BioPython 1.76
  * NumPy
  * Matplotlib
  * Scipy
  * some modules from the Python Standard Library (e.g. logging, random, multiprocessing, collections etc.)
* P2Rank 2.2 or later (only for part 3 of the tutorial)

<a name="one"></a>
## 1. Downloading data and running analysis
This can be done using python script `analysis_pipeline.py`.

### Input - Dataset file
We need to specify the list of structures in the dataset file. It is a plain-text file where one row equals one structure. The columns are separated by whitespace. The first two columns are mandatory and they contain PDB ID and chain ID (pipeline can only work with single-chain structures). The third column is optional and it can define a list of specific ligands which will be used for ligand binding sites computation.

Example:
- Dataset containing 3 structures:
    1s69	A
    1do1	A
    1qgw	C
    
- Dataset with 2 structures where only ligands with specified codes will be considered:
    1loj	B	URI
    1ja1	A	FAD,NAP,FMN

### Options and Arguments

### Examples



<a name="two"></a>
## 2. Running analysis on custom data

### Directory structure

<a name="three"></a>
## 3. Training a P2Rank model
