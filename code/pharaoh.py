# ***********
# A note from Sarit:
# This code runs a small example of an acetylome file found in /data/CBS_example.csv to a small subgroup of the mammals found in /data/referenced_proteomes_and_life_spans.csv. Both files can be modified. all mammals can be found at /data/referenced_proteomes_and_life_spans_ALL_MAMMALS.csv , acetylomes are attached to our paper.
# ***********

from ortholog_aligment import ortholog_aligment_main
from ortho_matrix import matrix_main
from statistics import stat_main
from read_trees_weights_by_distribution import trees_main
import os 

# files for code ocean run
data_dir_code_ocean="/data/"
db_dir_code_ocean = "/results/DBs/"
run_files_dir_code_ocean = "/results/running_files/"
results = "/results/"

# changble parameters
ptm_file = data_dir_code_ocean+"CBS_example.csv" #choose the name of the PTM file
the_animal_tsv_dir = data_dir_code_ocean+"Orthologues_Mus_musculus/" # choose eather "Orthologues_Mus_musculus" or "Orthologues_Homo_sapiens" 
mouse_or_human = 'Mus_musculus' # choose eather "Mus_musculus" or "Homo_sapiens" 
# Choose 2 groups of AA for the statistical validarion
first_char = ['K']
sec_char = ['R']

# permenant parameters - please don't change
table_name = "PTM_table" 
lifespan_table = "orthologs_life_spans"
lifespans_file = data_dir_code_ocean+"referenced_proteomes_and_life_spans.csv"
code_folder = os.path.dirname(os.path.realpath(__file__)) # current directory
matrix_output = data_dir_code_ocean+"output_distance.xlsx"
dict_file = data_dir_code_ocean+"an_age_to_tree_names.xlsx"
fatas_dir = data_dir_code_ocean+"/fastas/"



# my advice - if ypu run this on your computer, to better understand, comment out each part and run them saparatly :)

# run the ortholog alignment algorithm
# in this part all kinds of .db files will be created
print("Starting ortholog alignment")
ortholog_aligment_main(table_name, the_animal_tsv_dir, lifespans_file, ptm_file,db_dir_code_ocean,run_files_dir_code_ocean, fatas_dir)
print("Done ortholog alignment")

# create the replacement matrix
# in this part the replacement matrix will be created as "result.csv"
print("Starting replacement matrix")
matrix_main(table_name, code_folder, mouse_or_human,db_dir_code_ocean,run_files_dir_code_ocean)
print("Done replacement matrix")

# first part of statistical test - mann whitney and FDR
# in this part 2 files are created: "results_Mann-whitney.csv" - for the mann-whitney test, "results_fdr.csv" - additional FDR correction
print("Starting part 1 statistics")
stat_main(db_dir_code_ocean+table_name+'.db', lifespan_table, code_folder ,first_char, sec_char, run_files_dir_code_ocean)
print("Done part 1 statistics")

# second part of statistical test - tree validation
# in this part 2 files are created: "results_tree_dist.xlsx" - the tree statistical test results, "PHARAOH_output.xlsx" - FINAL FILE!!
print("Starting part 2 statistics")
trees_main(matrix_output,lifespans_file,dict_file, results, run_files_dir_code_ocean)
print("Done part 2 statistics")

