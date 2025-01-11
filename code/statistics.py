import csv
import numpy as np
import os
from scipy.stats import mannwhitneyu, ttest_ind, wilcoxon
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import sem
from database import Database

"""
update a csv file

Parameters
---------
file : str
    csv file location
row : list
 list of data as row
"""


def update_file(file, row):
    with open(file, mode='a', newline='') as csv_file:
        csv_file = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        # if headers:
        #     csv_file.writerow(headers)
        csv_file.writerow(row)


"""
get a list of all lifespans from the db table

Parameters
---------
base_db_name : str
    db name (a.db)
lise_spans_table : str
    table name

Returns
---------
life_spans_ints : list
    all life spans as float values
"""


def get_lifespans(base_db_name, lise_spans_table, unique_folder):
    # open the db
    table_name1 = "orthologs_life_spans"
    db_path = os.path.join(unique_folder, base_db_name)
    with Database(db_path) as db:
        # iterate over the table
        animals = db.iterate_over_table(table_name1)
        animals
        life_spans_list = []
        names_list = []
        # create a list
        for animal in animals:
            names_list.append(animal[0])
            life_spans_list.append(animal[1])
        #
        #print(names_list)
        #print("-------------------------------------------")

        #print(life_spans_list)
        life_spans_ints = np.asarray([float(i) for i in life_spans_list])
    return life_spans_ints, names_list


"""
mann_whitney U test
"""


def mann_whitney(group_zero, group_one):
    w, p = mannwhitneyu(group_zero, group_one, use_continuity=True, alternative='two-sided')
    return p


"""
FDR correction
input file = results_after_mann_whitney
output file = results_after_mann_whitney_and_after_fdr
"""


def fdr_correction(input_file, output_file):
    pvals = []
    all_data = []
    rows = []
    # get all p values
    with open(input_file, newline='') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        headers = next(csv_reader)
        #headers = headers + ['mann_whitney_fdr', 'fdr<0.1']
        headers = headers + ['mann_whitney_fdr']
        update_file(output_file, headers)
        for row in csv_reader:
            if row:
                rows.append(row)
                p = float(row[-1])
                if np.isnan(p):
                    p = 0
                pvals.append(p)
                all_data.append(row)
    # plot p values histogram
    #n, bins, patches = plt.hist(x=pvals, bins='auto', color='#0504aa',
    #                            alpha=0.7, rwidth=0.85)
    #plt.grid(axis='y', alpha=0.75)
    #plt.xlabel('P Value')
    #plt.ylabel('Count')
    # plt.show()

    # FDR correction to p values list
    out = fdrcorrection(np.asarray(pvals), alpha=0.1)

    # update the FDR q values to the output file
    out_len = len(out[0])
    for i in range(out_len):
        rows[i].append(str(out[1][i]))
        #rows[i].append(str(out[0][i]))
        update_file(output_file, rows[i])


"""
statistical_test each row in the results file, and update to another file
"""


def statistical_test(file, life_spans,updated_file, first_char, sec_char):
    flag = False
    # open the results file
    with open(file, newline='') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        headers = next(csv_reader)
        headers = headers + ['mann_whitney_p-value']
        update_file(updated_file, headers)
        # get data by rows
        for row in csv_reader:
            if row:
                data = row[5:]
                all_r = []
                x_r = []
                all_k = []
                x_k = []
                no_otrho = 0
                # Divide into 2 groups by the replacement and original chars
                for i in range(len(data)):
                    if data[i] in first_char:
                        all_k.append(life_spans[i])
                        x_k.append(i)
                    elif data[i] in sec_char:
                        all_r.append(life_spans[i])
                        x_r.append(i)
                    elif data[i] is '#' or "No K in rel_pos" in data[i]:
                        no_otrho += 1
                # calculate the mean of each group
                r_mean = np.mean(all_r)
                k_mean = np.mean(all_k)

                # statistical test
                # if the each group size > 3 update the p value to file
                if r_mean > k_mean and len(all_r) >= 3 and len(all_k) >= 3:
                    #print("result found!!!!!!!!")
                    flag = True
                    stat = mann_whitney(all_k, all_r)
                    row.append(str(stat))
                    update_file(updated_file, row)

    return flag



def stat_main(db_name, animals_table, unique_folder, first_aa, second_aa, run_files_dir_code_ocean):
    res_file = run_files_dir_code_ocean + 'result.csv'
    updated_file = run_files_dir_code_ocean + "results_Mann-whitney.csv"
    fdr_file = run_files_dir_code_ocean + "results_fdr.csv"

    life_spans, names_list = get_lifespans(db_name, animals_table, unique_folder)
    x = statistical_test(res_file, life_spans, updated_file, first_aa, second_aa)
    if x:
        fdr_correction(updated_file, fdr_file)
        return True
    else:
        return False


#if __name__ == "__main__":
#    stat_main('CBS_example.db', 'orthologs_life_spans', '/home/stu/cohenlab/sarit/project10.05.2021/final_code/','K', 'R')