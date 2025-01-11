import csv
import datetime
import os
import xml.dom.minidom
import more_itertools as mit
from itertools import groupby
from operator import itemgetter
from database import Database
import pandas as pd
import numpy as np
import time
from multiprocessing import Process


def write_csv(csv_filename, fields, data):
    with open(csv_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        writer.writerows(data)


def parse_xml_species_data(xml_filename):
    xml_doc = xml.dom.minidom.parse(xml_filename)
    # xml_doc = xml.dom.minidom.parse(xml_filename)
    # parse_xml_species_data(xml_doc)

    species1_dict = {}
    species2_dict = {}

    # print('\n[' + datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S") +
    #       '] Handling XML species configuration')
    # f_log.write('[' + datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S") +
    #             '] Handling XML species configuration.\n')

    id_to_protein = {}
    species_items = xml_doc.getElementsByTagName('species')
    species1_name = species_items[0].getAttribute('name')
    species1_genes = species_items[0].getElementsByTagName('gene')
    species2_name = species_items[1].getAttribute('name')
    species2_genes = species_items[1].getElementsByTagName('gene')
    for gene in species1_genes:
        id = gene.getAttribute('id')
        protId = gene.getAttribute('protId')
        geneId = gene.getAttribute('geneId')

        gene_item = (species1_name, id, protId, geneId)
        # print(gene_item)
        species1_dict[id] = gene_item
        id_to_protein[id] = protId

    for gene in species2_genes:
        id = gene.getAttribute('id')
        protId = gene.getAttribute('protId')
        geneId = gene.getAttribute('geneId')

        gene_item = (species2_name, id, protId, geneId)
        # print(gene_item)
        species2_dict[id] = gene_item
        id_to_protein[id] = protId

    # put here you logic to combine both species from dictionaries and sent to csv writer

    # print('\n[' + datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S") +
    #       '] Create CSV file ' + csv_filename)

    species1 = [int(i) for i in list(species1_dict.keys())]
    species2 = [int(i) for i in list(species2_dict.keys())]

    grplist1 = [list(group) for group in mit.consecutive_groups((species1))]
    for lis in grplist1:
        if len(lis) > 1:
            lis
    grplist2 = [list(group) for group in mit.consecutive_groups((species2))]

    orthologs = {}
    ziped = list(zip(grplist1, grplist2))
    for pair in ziped:
        s1 = pair[0]
        s2 = pair[1]
        if len(s1) > 1 and len(s2) == 1:
            for p in s1:
                orthologs[str(p)] = str(s2[0])
        elif len(s2) > 1 and len(s1) == 1:
            for p in s2:
                orthologs[str(s1[0])] = str(p)
        elif len(s1) == 1 and len(s2) == 1:
            orthologs[str(s1[0])] = str(s2[0])
        else:
            for i in s1:
                for j in s2:
                    orthologs[str(i)] = str(j)

    ortho_pairs = find_orthologs(orthologs, id_to_protein, species1_name, species2_name)

    return ortho_pairs

    # create_all_dict(species1_dict, species2_dict)
    # write_csv(csv_filename, fields, species_list)


def find_orthologs(orthologs, id_to_protein, species1_name, species2_name):
    ortho_pairs = {}
    for s1, s2 in orthologs.items():
        id1 = id_to_protein[str(s1)]
        id2 = id_to_protein[str(s2)]
        if "Homo sapiens" in species1_name:
            ortho_pairs[id1] = id2
        else:
            ortho_pairs[id2] = id1
    return ortho_pairs


def create_all_dict(species1_dict, species2_dict):
    all_dict = {}
    key_prot_dict = {}
    i = 1
    j = 1
    # id1 = species1_dict[str(i)][1]
    # id2 = species2_dict[str(j)][1]
    size = len(species1_dict) + len(species2_dict) + 1
    for id in range(1, size, 1):
        if str(id) in species1_dict:

            all_dict[str(id)] = species1_dict[str(id)]
            key_prot_dict[species1_dict[str(id)][2]] = str(id)
        else:
            all_dict[str(id)] = species2_dict[str(id)]
            key_prot_dict[species2_dict[str(id)][2]] = str(id)

    return all_dict, key_prot_dict



def human_to_mouse(mouse_ortho_pairs, other_ortho_pairs, human_or_mouse):
    all_orthologs = {}
    if human_or_mouse == "Mus_musculus":

        for key, val in mouse_ortho_pairs.items():
            if key in other_ortho_pairs:
                all_orthologs[val] = other_ortho_pairs[key]
    else:
        for key, val in mouse_ortho_pairs.items():
            if key in other_ortho_pairs:
                all_orthologs[key] = other_ortho_pairs[val]

    return all_orthologs


def find_all_orthologs(ortho_dict, db_dir, mouse_ortho_pairs, file, human_or_mouse):
    if "Mus musculus" not in file:
        if "sapiens" not in file:
            name = file.split('.')[0].replace(" ", "_")
            fila_path = ortho_dict + file
            animal_ortho_pairs = parse_xml_species_data(fila_path)
            all_orthologs = human_to_mouse(mouse_ortho_pairs, animal_ortho_pairs, human_or_mouse)
            # get acetylome
            db_path = os.path.join(db_dir, name + '.db')
            with Database(db_path) as db:
                table = name + "_all_rep"
                rows = db.iterate_over_table(table)
                ortho_table = name + "_orthologs_results"
                db.create_ac_rep_table(ortho_table)
                find_orthologs_in_db(rows, all_orthologs, db, ortho_table)
                db.close()


def find_orthologs_in_db(rows, animal_ortho_pairs, db, ortho_table):
    for key, value in animal_ortho_pairs.items():
        matching_key = [s for s in rows if key in s]
        if matching_key:
            matching_value = [s for s in matching_key if value in s]
            if matching_value:
                # new_rows.extend(matching_value)
                for m in matching_value:
                    to_db = list(m)
                    to_db = [str(i) for i in to_db]
                    db.insert_row_into_rep_table(ortho_table, [to_db])
    # return new_rows


def find_all_orthologs_in_human(ortho_dict, mouse_xml, db_dir, human_or_mouse):
    # get human-mouse
    mouse_ortho_pairs = parse_xml_species_data(mouse_xml)
    # iterate over the directories
    for subdir, dirs, files in os.walk(ortho_dict):
        if human_or_mouse == "Mus_musculus":
            name = "Homo_sapiens"
        else:
            name = "Mus_musculus"
        all_orthologs = {}
        for key, val in mouse_ortho_pairs.items():
            all_orthologs[val] = key
        db_path = os.path.join(db_dir, name + '.db')
        with Database(db_path) as db:
            table = name + "_all_rep"
            rows = db.iterate_over_table(table)
            ortho_table = name + "_orthologs_results"
            db.create_ac_rep_table(ortho_table)
            find_orthologs_in_db(rows, all_orthologs, db, ortho_table)
            db.close()


def ortho_matrix(base_db_name, base_table, db_dir, m_file, animals, human_or_mouse, db_dir_code_ocean):
    final_results = []
    
    # get the PRM file rows and animals list
    with Database(base_db_name) as db:
        ac_rows = db.iterate_over_table(base_table)
        animals = db.iterate_over_table(animals)
        db.close()
        
    # iterate over the animals, for each animal get the results
    for animal in animals:
        try:
            name = animal[0].replace(' ', '_')
            #print("--------------------------------------")
            #print(name)
            #print(human_or_mouse)
            #print("--------------------------------------")

            if name == human_or_mouse:
                results = [human_or_mouse]

                for row in ac_rows:
                    seq = row[3]
                    rel_pos = row[5]
                    res = seq[int(rel_pos)].upper()
                    # print("----------------------- im res@@@@@@@@@@@@@@@@@@@ ------------------")
                    # print(res)
                    results.append(res)
                    # print("----------------------- im results!!!!! ------------------")
                    # print(results)
                final_results.append(results)
                #print("----------------> risults are ready!!!! <-------------------")
                #print(final_results)

            else:
                name = animal[0].replace(' ', '_')
                db_path = os.path.join(db_dir_code_ocean, name + '.db')
                with Database(db_path) as db:
                    ortho_table = name + "_all_rep"
                    rows = db.iterate_over_table(ortho_table)
                    if rows:
                        animal_results = [name]
                        animal_results.extend(find_orthologs_to_matrix(rows, ac_rows, db, ortho_table))
                        final_results.append(animal_results)
                        db.close()

        except Exception as e:
            print(e)


    # create the data structure of the final results
    numpy_array = np.array(final_results)
    transposed = numpy_array.T
    transposed_list = transposed.tolist()
    temp_list = []

    for i in range(len(transposed_list)-1):
        row = list(transposed_list[i+1])
        ac_row = list(ac_rows[i])
        n = [ac_row[0], ac_row[6], ac_row[2], ac_row[4], ac_row[3]]
        n.extend(row)
        temp_list.append(n)
        
    # create a results file
    first_row = ['Key', 'Uniprot_ID', 'Ref', 'Pos', 'Seq']
    first_row.extend(transposed_list[0])
    update_file(m_file, first_row)

    #write all to a csv results file
    for row in temp_list:
        update_file(m_file, row)




def update_file(m_file, row):
    with open(m_file, mode='a', newline='') as csv_file:
        csv_file = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csv_file.writerow(row)


def find_orthologs_to_matrix(rows, ac_rows, db, ortho_table):
    results = []
    data = pd.DataFrame(rows,
                        columns=["Animal", "Acetylome_id", "Protein_id", "Pos", "Identities", "Positives", "GN",
                                 "Seq", "Acetylome_seq", "Protein_seq", "Replacement"])
    for ac_site in ac_rows:
        id = ac_site[6]
        seq = ac_site[3]
        site = ac_site[4]
        found = data.loc[(data['Acetylome_id'] == id) & (data['Pos'] == site) & (data['Seq'] == seq)]
        # find the result code for this row
        if found.empty:
            results.append("#")
        else:
            m_res = found['Replacement'].to_string(index=False)
            results.append(m_res)
    return results


def find_in_rep_code(found, original_char, rep):
    i = 0
    best_pos = 60
    rep_found = False
    original_char_found = False
    m_rep = ""
    # if a result was found

    # sort the results by the positive values
    found['Positives'] = found['Positives'].astype(int)
    sorted_results = found.sort_values(by=['Positives'], ascending=False)
    # iter over the result
    for index, m_row in sorted_results.iterrows():
        m_positives = int(m_row['Positives'])
        # if the positive value is above 60:
        if m_positives >= 60:
            # get the highest positive value
            if i == 0:
                best_pos = m_positives
            # if its in the top five or has the best positive
            if m_positives == best_pos:
                m_rep = m_row['Replacement']
                # if we found the rep char
                for r in rep:
                    if m_rep == r:
                        rep_found = True
                    # if we found the original char
                    elif m_rep == original_char:
                        original_char_found = True
                        # print(original_char)
            else:
                break
        i += 1
    # if we found the rep char - 1
    # if we didn't found the rep char but found the original char - 0
    # else - we didnt found a result or one of the chars - *
    return m_rep.upper()


def matrix_main(base_table, unique_folder, human_or_mouse, db_dir_code_ocean,run_files_dir_code_ocean):
    base_db_name = db_dir_code_ocean + base_table +".db"
    matrix_file = run_files_dir_code_ocean + 'result.csv'
    animals = "orthologs_life_spans"

    ortho_matrix(base_db_name, base_table, unique_folder, matrix_file, animals, human_or_mouse, db_dir_code_ocean)
    #print("All done creating results.csv")

#if __name__ == "__main__":
#  matrix_main('CBS_example', "/home/stu/cohenlab/sarit/project10.05.2021/final_code/", 'Mus_musculus')
