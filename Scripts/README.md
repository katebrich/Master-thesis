# Pipeline for statistical analysis of ligand binding sites properties
This folder contains scripts for analysis of ligand binding sites properties (called features). The pipeline is able to download PDB and FASTA files from PDBe, label binding residues based on the ligands in the PDB file, compute values of defined features, perform statistical analysis to find out whether these features are interesting for the recognition of binding sites, and draw various histograms and plots. Furthermore, the extended pipeline can use obtained PDB files and features to train a new model in [P2Rank](http://siret.ms.mff.cuni.cz/p2rank) predictor of binding sites to see if the new features improve the prediction.

Users can define a feature of their own by implementing a method for its computation.

The only needed input is a dataset file with proteins.

Below, you can learn how to:
 - [download data from database and run statistical analysis with predefined features](#one)
 - [add a custom feature](#two)
 - [run analysis on custom data obtained from a different source](#three)
 - [train a P2Rank model with new features](#four)


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

    ```
    1s69	A
    1do1	A
    1qgw	C
    ```
    
- Dataset with 2 structures where only ligands with specified codes will be considered:

    ```
    1loj	B	URI
    1ja1	A	FAD,NAP,FMN
    ```
Several dataset can be found [here](../datasets)

### Output
All output files are located in an output root folder which was specified by argument `-o` or `--output_dir`. 

The root folder structure is as follows:

```
output_path
│   run.log  
│
└───PDB
│       1s69A.pdb
│       1do1A.pdb
|       ...  
└───FASTA
|       1s69A.fasta
|       1do1A.fasta
|        ...  
└───mappings    
|       1s69A.txt
|       1do1A.txt
|        ...
└───lbs    
|       1s69A.txt
|       1do1A.txt
|        ...
└───features    
|   |
│   └───feature1
│   |       1s69A.txt
|   |       1do1A.txt
│   |       ...
│   └───feature2
│           1s69A.txt
|           1do1A.txt
│           ...
└───analysis    
    |   binding_ratios.csv
    |   errors.txt
    |   means_difference.csv
    |   p_vals_perc.csv
    |   p_values.csv
    |   p_values_means.csv
    |
    └───feature1
    |       iterations.txt
    |       pairs.txt
    |       p_values.txt
    |       plots, histograms
    |       ...
    └───feature2
            iterations.txt
            pairs.txt
            p_values.txt
            plots, histograms
            ...

```

- File `run.log` is copied to the output directory when the program terminates. This is the detailed log; the brief log is printed to the console as the program runs.
- Folder `mappings` contains cached residue mappings; author residue number (plus insertion code, if any) from the PDB file is in the first column. In the second column, there is PDB molecule number for the residue.
- Folder `lbs` contains labeling of binding sites. For each structure, there is a file with one line per residue, where the first number is PDB molecule residue number and the second is label 0/1.
- The subfolders of `features` contain computed feature values. Again, the first number is PDB molecule residue number and the second number is the feature value for the residue.
- `analysis` folder contains several summary files. File `errors.txt` lists features where the analysis ended with error. This is often caused by lack of data (e.g. the sample size is bigger than number of rows or the data for a categorical feature are too sparse to meet the assumptions of Chi-squared test). Detailed information about the errors can be found in log.
- Furthermore, there are more statistics and graphs in separate folders for every feature. `pairs.txt` file contains paired ligand binding sites labels with feature values for all the structures in the dataset. It could be useful for more analysis or e.g. for training and testing a classifier.

### Options and Arguments

```
Usage: analysis_pipeline.py -d DATASET_FILE_PATH -o OUTPUT_DIR_PATH [OPTIONS]... 

Mandatory arguments to long options are mandatory for short options too.
  -d, --dataset                Mandatory; file with listed structures to process
  -o, --output_dir             Mandatory; root folder. Created if not exists
  -t, --tasks                  Default: 'A'. Comma-separated list of tasks to process. If data are missing in root folder for some task, they are computed even if their task is not in the list. Possible values: 'D' - download; 'L' - compute ligand binding sites; 'F' - compute features; 'A' - compute analysis
  -f, --force                  if an existing destination file cannot be





```

### Examples


<a name="two"></a>
## 2. Defining new features

<a name="three"></a>
## 2. Running analysis on custom data

pozor na pojmenovani filu...
### Directory structure

<a name="four"></a>
## 3. Training a P2Rank model
