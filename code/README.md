# PHARAOH tool

## Overview
This project offers a comprehensive pipeline for the analysis of orthologous proteins, creation of a replacement matrix, and conducting statistical tests to identify correlations between amino acid changes and lifespan. The pipeline is specifically tailored for Mus musculus or Homo sapiens PTM (Post-Translational Modification) data.

## Data Preparation
- **Proteome Data:** All mammalian proteome FASTA files are available in the `fastas` directory.
- **Ortholog Pairs:** "OrthoFinder" based Ortholog pairs are in either `Orthologues_Mus_musculus` or `Orthologues_Homo_sapiens` directories.
- **Species Tree:** The phylogenetic tree based on the orthologs can be found in `SpeciesTree_rooted_orthoFinder.txt`.
- **Evolutionary Distance Matrix:** The evolutionary distance matrix between each mammal is stored in `output_distance.xlsx`. (Refer to `read_trees_weights_by_distribution.py` for the code creating this distance matrix.)
- **Raw acetylome data:** both human and mouse acetylome data can be found in this directory.


## Usage
The primary script for executing the entire pipeline is `pharaoh.py`. Running this script initiates the complete workflow, including ortholog alignment, replacement matrix creation, and statistical tests.

## Prerequisites
- Python 3.9.16
- Required Python packages:
  - csv
  - sqlite3
  - os.path
  - glob
  - random
  - Bio
  - Phylo
  - matplotlib
  - numpy
  - pandas
  - openpyxl
  - statistics
  - statsmodels
  - math
  - subprocess
  - threading
  - swalign
  - concurrent.futures
  - time
  - multiprocessing
  - datetime
  - SeqIO
  - more_itertools
  - itertools
  - groupby
  - itemgetter
 
## Configuration
Key parameters can be adjusted in the script for customization:
- `ptm_file`: PTM file of your choise followint the example file (default: `CBS_example.csv`)
- `the_animal_tsv_dir`: Choose eather "Orthologues_Mus_musculus" or "Orthologues_Homo_sapiens"  (default: `Orthologues_Mus_musculus`)
- `mouse_or_human`: Organism selection - `Mus_musculus` or `Homo_sapiens` (default: `Mus_musculus`)
- `first_char`: List of first characters for statistical tests (default: `['K']`)
- `sec_char`: List of second characters for statistical tests (default: `['R']`)

## Output
The pipeline generates various output files, including ortholog alignment databases, the replacement matrix (`result.csv`), statistical test results (`results_Mann-whitney.csv`, `results_fdr.csv`), and tree validation results (`results_tree_dist.xlsx`, `PHARAOH_output.xlsx`). 
"PHARAOH_output.xlsx" is the final results file.

## Note
It is advisable to review the comments within the code for a better understanding of each section's purpose and functionality.



