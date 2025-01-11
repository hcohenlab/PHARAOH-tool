import math
import os
import sys
from os import walk
import glob
from database import Database
import subprocess
import csv
from Bio.Blast import NCBIXML
import threading
from threading import Thread, Lock
import swalign

import concurrent.futures
import time
import multiprocessing as mp
from multiprocessing import Pool
from multiprocessing import Process
import cProfile
from datetime import datetime
import pandas as pd
from Bio import SeqIO


mutex = Lock()
results_list = []

"""
create a file
write to the file - name, content of file: seq
use mutex for protection of the file

Parameters
---------
name : str
  file location and name
seq : str
  file content
"""


def create_seq_file(name, seq):
  mutex.acquire()
  try:
      f = open(name, "w")
      f.write(seq)
      f.close()
  finally:
      mutex.release()


"""
the main loop of the algorithm
iterate over all proteins in DB
find for each protein matching proteins from blast
for each match look at the acetylation site for replacements
write all results into a file and db table

Parameters
---------
db_name : str
  name of the DB
db_root_loc : str
  location of the dit that contains all the dbs
input_loc, out_loc : str
  file locations for the sequences, used in blast
results_file : str
  file location for the blast results
table : str
  table name from the db
"""


def find_ac_site(db, rows, res_table, df, animal, fatas_dir):
  # get a row from the db
  row = rows.pop(0)
  # get all the data about the protein from the DB
  ac_name = row[2].upper()
  seq = row[3]
  pos = row[4]
  rel_pos = row[5]
  id = row[6]
  biggest_score = -math.inf

  found = df[df["Mus_musculus"].str.contains(id)][animal]
  temp_list1 = found.to_string(index=False).split(',')
  if not found.empty:
    try:
      for element in temp_list1:
        sequences = [i for i in SeqIO.parse(fatas_dir + animal + '.fasta','fasta')]
        for seq1 in sequences:
          #print("IM element->seq1 ------------------> " + element + "-----------------------------------------------> " + seq1.name)
          if element.replace(" ", "") == seq1.name.replace(" ", ""):
            #print("We Are In!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            match = 2
            mismatch = -1
            gap = -2
            scoring = swalign.NucleotideScoringMatrix(match, mismatch)
            sw = swalign.LocalAlignment(scoring, gap)
            aln = sw.align(seq1.seq, seq)
            score1 = aln.score
            if score1 >= biggest_score:
              seq_to_save = seq1
              biggest_score = score1
              best_aln = aln.identity
              all_matches = aln.matches
              position = abs(int(rel_pos) - aln.q_pos)
              new_k = aln.ref[aln.r_pos + position]
            break

      temp_list3 = [animal, id, pos, seq_to_save.id.split("|")[1], best_aln, all_matches, ac_name, seq, seq, str(seq_to_save.seq) ,new_k]
      update_rep_row_to_db(db, [temp_list3], res_table)

    # catch any exceptions and add them to the results file
    except Exception as e:
      print(e)

"""
update a row into a database

Parameters
---------
db : object
  db object for the update
m_list : list
  a row as a list
table : str
  table name for update
results_file : str
"""
def update_rep_row_to_db(db, m_list, table):
  db.insert_row_into_rep_table(table, m_list)
  
  
"""
run the loop of the algorithm for each aminal

Parameters
---------
animal : str
  name of the animal
table_name : str
  name of the PTM table in the DB
the_animal_tsv_dir : str
  directory containing all orthologs files
"""
def run_loop(animal, table_name, the_animal_tsv_dir, db_dir_code_ocean, fatas_dir):
  print("running alignment to:")
  print(animal)
  res_table = animal + "_all_rep"
  db_name = db_dir_code_ocean + animal + ".db"
  df = pd.DataFrame()

  # get the rows from the PTM file
  with Database(db_dir_code_ocean + table_name+'.db') as db:
      rows = db.iterate_over_table(table_name)
      db.close()
      
  # for this specific animal, get the orthologus proteins
  for subdir, dirs, files in os.walk(the_animal_tsv_dir):
    for file in files:
      #print(file)
      if animal in file:
        df = pd.read_csv(the_animal_tsv_dir + file, sep='\t')
        break

  # find the replacements for each site
  with Database(db_name) as db:
      db.create_ac_rep_table(res_table)
      while len(rows) > 0:
        find_ac_site(db, rows, res_table, df, animal, fatas_dir)


# *****************************************MAIN******************************8
def ortholog_aligment_main(table_name, the_animal_tsv_dir, lifespans_file, ptm_file,db_dir_code_ocean,run_files_dir_code_ocean, fatas_dir):
  procs = []
  db_name = db_dir_code_ocean+table_name+'.db'
  with Database(db_name) as db:
    db.create_species_life_spans_table("orthologs_life_spans")
    db.insert_species_life_spans_csv_to_db(lifespans_file, "orthologs_life_spans")
    
    
    db.create_acetylom_table(table_name)
    db.insert_acetylom_csv_to_db(ptm_file, table_name)
    
    rows = db.iterate_over_table("orthologs_life_spans")
    db.close()

  animals_list = []

  for row in rows:
    animals_list.append(row[0].replace(' ', '_'))


  # give for any animal it own process to run
  for animal in animals_list:
    #sequencial run
    run_loop(animal, table_name, the_animal_tsv_dir,db_dir_code_ocean, fatas_dir)
    
    # Multyprocess run
    #proc = Process(target=run_loop, args=(animal, table_name, the_animal_tsv_dir))
    #procs.append(proc)
    #proc.start()

  #for proc in procs:
    #proc.join()




#if __name__ == "__main__":
#  ortholog_aligment_main('CBS_example', "/home/stu/cohenlab/sarit/project10.05.2021/final_code/Orthologues_Mus_musculus/", 'referenced_proteomes_and_life_spans.csv', 'CBS_example.csv')
