import random
from Bio import Phylo
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import openpyxl
from statistics import update_file
from statsmodels.stats.multitest import fdrcorrection

def lookup_by_names(tree):
    """ get the list of animals in a phylogenetic tree

    :param tree: a tree object
    :return: list of names of animals in tree
    """
    names = []
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names.append(clade.name)
    return names


def calculate_distance_to_matrix(tree, output):
    """ calculate a distance matrix using a tree input

    :param tree:
    :param output: a tree object
    :return: output file for the distance matrix
    """
    names = lookup_by_names(tree)
    matrix = np.zeros((len(names), len(names)))
    for i in range(len(names)):
        for j in range(len(names)):
            matrix[i][j] = tree.distance(names[i], names[j])

    matrix_df = pd.DataFrame(matrix, columns=names)
    matrix_df.to_excel(output)


def get_mean_lifespan_and_distance(matrix, animals, condition):
    """ calculate the mean of lifespans and distances

    :param matrix: file containing distance matrix
    :param animals: file containing spices and lifespans
    :return: mean of median for each file
    """
    matrix_df = pd.read_excel(matrix, index_col=0)
    matrix_df = matrix_df.div(matrix_df.to_numpy().max())
    matrix_mean = np.nanmean(matrix_df)
    matrix_median = np.median(matrix_df.values)

    animals_df = pd.read_csv(animals)

    animals_mean = animals_df['lifespan'].div(animals_df['lifespan'].to_numpy().max()).mean()
    animals_median = animals_df['lifespan'].div(animals_df['lifespan'].to_numpy().max()).median()

    if condition == "mean":
        return matrix_mean, animals_mean
    else:
        return matrix_median, animals_median


def normalize_df(df):
    """

    :param df: numeric dataframe
    :return: max val of df
    """

    return df.to_numpy().max()


def an_age_to_tree(file):
    dict = {}
    df = pd.read_excel(file)
    # i = 0
    for index, row in df.iterrows():
        dict[row[0].replace(" ", "_")] = row[1].replace(" ", "_")
    return dict


def effect_on_test(matrix_output, lifespans, results_file, dict_file, output):
    final_scores = []
    final_sig = []
    pvals = []
    # get mean life spans and distance
    matrix_mean, animals_mean = get_mean_lifespan_and_distance(matrix_output, lifespans, "mean")
    # get a normalizes distance matrix from tree
    distance_matrix = pd.read_excel(matrix_output, index_col=0)
    norm_distance_matrix = distance_matrix.div(distance_matrix.to_numpy().max())

    # tree_animals = pd.read_excel(matrix_output).columns.values
    animals_lifespans_df = pd.read_csv(lifespans)
    animals = animals_lifespans_df["species"]

    animals = [s.replace(' ', '_') for s in animals]

    norm_animals_lifespans = animals_lifespans_df
    norm_animals_lifespans["species"] = animals
    norm_animals_lifespans['lifespan'] = norm_animals_lifespans['lifespan'].div(
        norm_animals_lifespans['lifespan'].to_numpy().max()).transpose()
    animals_dict = an_age_to_tree(dict_file)

    # get weights
    weights = calculate_weights(norm_distance_matrix, norm_animals_lifespans, animals, animals_dict, matrix_mean,
                                animals_mean)

    results_df = pd.read_csv(results_file)
    proteins = results_df['Uniprot_ID'].to_list()
    pos = results_df['Pos'].to_list()

    for i in range(len(results_df.index)):
        #print(proteins[i] + "  " + str(pos[i]))
        found_df = results_df[(results_df['Uniprot_ID'] == proteins[i]) & (results_df['Pos'] == pos[i])]
        score, sig, persent = calculate_score(norm_distance_matrix, norm_animals_lifespans, animals, animals_dict,
                                              found_df, matrix_mean,
                                              animals_mean, weights)
        #print(sig)
        final_scores.append(score)
        pvals.append(persent)
        final_sig.append(str(sig))

    results_df["tree_scores"] = final_scores
    results_df["tree_p_val"] = pvals
    #results_df["new_sig"] = final_sig
    results_df.to_excel(output, index=False)


def calculate_score(norm_distance_matrix, norm_animals_lifespans, animals, animals_dict, found_df, matrix_mean,
                    animals_mean, weights):
    """

    :param weights: [big_small, small_big,small_small, big_big]
    :param norm_distance_matrix:
    :param norm_animals_lifespans:
    :param animals:
    :param animals_dict:
    :param found_df:
    :param matrix_mean:
    :param animals_mean:
    :return:
    """
    score = 0
    scores_sum = []
    firsr_mammal = []
    sec_mammal = []
    effects = []
    final_scores_chanes = []
    # for each pair of mammals
    for animal1 in animals:
        for animal2 in animals:
            # check if same aa
            aa = False
            if found_df[animal1].values[0] == found_df[animal2].values[0]:
                aa = True
            # get the name for the tree, and distances
            tree_animal1 = animals_dict[animal1]
            tree_animal2 = animals_dict[animal2]
            distance = norm_distance_matrix[tree_animal1][tree_animal2]
            # get delta lifespans
            lifespan1 = norm_animals_lifespans[norm_animals_lifespans["species"] == animal1]['lifespan'].values[0]
            lifespan2 = norm_animals_lifespans[norm_animals_lifespans["species"] == animal2]['lifespan'].values[0]
            dif_lifespans = np.abs(lifespan1 - lifespan2)

            # calculate the score
            # Strengthens claim
            if aa and distance >= matrix_mean and dif_lifespans < animals_mean:
                # big_small
                effect = weights[0]

            elif not aa and distance < matrix_mean and dif_lifespans >= animals_mean:
                # small_big
                effect = weights[1]
                # effect = 1/2

            # weakens the claim
            elif not aa and distance < matrix_mean and dif_lifespans < animals_mean:
                # small_small
                effect = -weights[2]
                # effect = -1/4
            elif not aa and distance >= matrix_mean and dif_lifespans < animals_mean:
                # big_small
                effect = -weights[0]
                # effect = -1/4

            elif aa and distance < matrix_mean and dif_lifespans >= animals_mean:
                # small_big
                effect = -weights[1]
                # effect = -1/4

            elif aa and distance >= matrix_mean and dif_lifespans > animals_mean:
                # big_big
                effect = -weights[3]
                # effect = -1/4

            else:
                effect = 0

            score += dif_lifespans * distance * effect
            scores_sum.append(dif_lifespans * distance)

            # ---------------for tests---------------------------------
            firsr_mammal.append(animal1)
            sec_mammal.append(animal2)
            effects.append(dif_lifespans * distance * effect)
            final_scores_chanes.append(score)

    df = pd.DataFrame()
    df['animals1'] = firsr_mammal
    df['animals2'] = sec_mammal
    df['score'] = effects
    df['total_score'] = final_scores_chanes
    # df.to_excel(R"C:\temp\significant_test.xlsx")

    # plotting
    # labels = []
    # for i in range(len(firsr_mammal)):
    #     labels.append(firsr_mammal[i] + sec_mammal[i])
    # y = range(0, len(effects))
    # fig, axs = plt.subplots(2, sharex=True)
    # axs[0].plot(y, effects)
    # axs[0].set_title("Score changes")
    # axs[1].plot(y, final_scores_chanes)
    # axs[1].set_title("Total score")
    # # axs[1].set_xlabel(labels)
    # # plt.xticks(rotation=90)
    # plt.show()
    # ---------------for tests---------------------------------

    curve, sig, persent = calculate_pval(scores_sum, score, 1000, weights)
    # plot histogram
    #plt.hist(curve)
    #plt.show()
    return score, sig, persent


def calculate_pval(scores_sum, original_score, permutation_size, all_weights):
    """
    statistical permutation test
    assign randome values from the list as strenghtening or weakning the test
    :param weights:
    :param scores_sum:
    :param original_score:
    :param permutation_size:
    :return:
    """
    # weights = [-1, 1]
    # weights = [0.5, 0.5, -0.25, -0.25, -0.25, -0.25]
    weights = [all_weights[0], all_weights[1], -all_weights[2], -all_weights[0], -all_weights[1], -all_weights[3]]
    curve = []
    more_then_score = 0
    sig = False
    for per in range(permutation_size):
        score = 0
        for i in range(len(scores_sum)):
            choise = random.choice(weights)
            score += scores_sum[i] * choise
        if score < original_score:
            more_then_score += 1
        curve.append(score)

    persent = 1 - more_then_score / permutation_size
    if persent <= 0.05:
        sig = True
    return curve, sig, persent


def calculate_weights(norm_distance_matrix, norm_animals_lifespans, animals, animals_dict, matrix_mean, animals_mean):
    """
    get the
    :param norm_distance_matrix:
    :param norm_animals_lifespans:
    :param animals:
    :param animals_dict:
    :param matrix_mean:
    :param animals_mean:
    :return:
    """
    big_small = 0
    small_big = 0
    small_small = 0
    big_big = 0
    other = 0

    # for each pair of mammals
    for animal1 in animals:
        for animal2 in animals:
            # get the name for the tree, and distances
            tree_animal1 = animals_dict[animal1]
            tree_animal2 = animals_dict[animal2]
            distance = norm_distance_matrix[tree_animal1][tree_animal2]
            # get delta lifespans
            lifespan1 = norm_animals_lifespans[norm_animals_lifespans["species"] == animal1]['lifespan'].values[0]
            lifespan2 = norm_animals_lifespans[norm_animals_lifespans["species"] == animal2]['lifespan'].values[0]
            dif_lifespans = np.abs(lifespan1 - lifespan2)

            # calculate the score
            if animal1 != animal2:
                if distance >= matrix_mean and dif_lifespans < animals_mean:
                    big_small += 1
                elif distance < matrix_mean and dif_lifespans >= animals_mean:
                    small_big += 1
                elif distance < matrix_mean and dif_lifespans < animals_mean:
                    small_small += 1
                elif distance >= matrix_mean and dif_lifespans > animals_mean:
                    big_big += 1
                else:
                    other += 1

    all = big_small + small_big + small_small + big_big + other
    return [big_small / all, small_big / all, small_small / all, big_big / all]






"""
FDR correction
input file = results_after_mann_whitney
output file = results_after_mann_whitney_and_after_fdr
"""


def fdr_correction(input_file, output_file):
    all_data = []
    rows = []
    # get all p values
    df = pd.read_excel(input_file)
    pvals = df["tree_p_val"]  
    
    # FDR correction to p values list
    out = fdrcorrection(np.asarray(pvals), alpha=0.1)
    #print(out)

    # update the FDR q values to the output file
    out_len = len(out[0])
    all_q = []
    all_sig = []
    for i in range(out_len):
        all_q.append(str(out[1][i]))
        all_sig.append(str(out[0][i]))
    #print(all_q)
    df["q_val"] = all_q
    #df["fdr"] = all_sig
    df.to_excel(output_file, index=False)
    
        
        
def trees_main(matrix_output,lifespans,dict_file, result, run_files_dir_code_ocean):
    # caculate distance matrix - already done for the orthofinder tree 
    # matrix_output = "output_distance.xlsx"
    # calculate_distance_to_matrix(tree, matrix_output)

    results_file = run_files_dir_code_ocean + "results_fdr.csv"
    dist_res = run_files_dir_code_ocean + "results_tree_dist.xlsx"
    final_res = result + "PHARAOH_output.xlsx"
    effect_on_test(matrix_output, lifespans, results_file, dict_file, dist_res)
    fdr_correction(dist_res, final_res)


#if __name__ == '__main__':
#    trees_main(matrix_output,lifespans,dict_file)
