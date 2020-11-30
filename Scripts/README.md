# Pipeline for statistical analysis of ligand-binding sites properties
This folder contains scripts for analysis of ligand-binding sites properties (called features). The pipeline is able to download PDB and FASTA files from PDBe, label binding residues based on the ligands in the PDB file, compute values of defined features, perform statistical analysis to find out whether these features are interesting for the recognition of binding sites, and draw various histograms and plots. Furthermore, the extended pipeline can use obtained PDB files and features to train a new model in [P2Rank](http://siret.ms.mff.cuni.cz/p2rank) predictor of binding sites to see if the new features improve the prediction.

Users can define a feature of their own by implementing a method for its computation.

The only needed input is a dataset file with proteins.

In this tutorial, it is described how to:
 - [download data from database and run statistical analysis](#one)
 - [run analysis on custom data obtained from a different source](#two)
 - [train a P2Rank model with new features](#three)


### Prerequisities:
- xxx
- zzz

<a name="one"></a>
## 1. Downloading data and running analysis

### Options and Arguments

### Dataset file

### Directory structure

### Examples
xxxxxxxxxxxxxxxxxxxxx
x
x
x
x
x
x
x
x
x
x
x
x
x
x
<a name="two"></a>
## 2. Running analysis on custom data


<a name="three"></a>
## 3. Training a P2Rank model
