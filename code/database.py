import csv
import sqlite3
import os.path
import glob

"""
a DB class with general sqlite3 functions
"""


class Database:
    def __init__(self, DBname):
        """
        connection = sqlite3.connect(self.DBname)
        cursor = connection.cursor()
        """
        self.DBname = DBname
        self._conn = sqlite3.connect(DBname, timeout=100.0)
        self._cursor = self._conn.cursor()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close

    @property
    def get_cursor(self):
        return self._cursor

    @property
    def get_connection(self):
        return self._conn

    def commit(self):
        self.get_connection.commit()

    def close(self, commit=True):
        if commit:
            self.commit
        self.get_connection.close()

    def execute(self, sql, params=None):
        self.get_cursor.fetchall()

    def executemany(self, cmd, to_db):
        self.get_cursor.executemany(cmd, to_db)

    def fetchall(self):
        return self.get_cursor.fetchone()

    def query(self, sql, params=None):
        self.get_cursor.execute(sql, params or ())
        return self.fetchall()

    def query_fetch_all(self, sql, params=None):
        self.get_cursor.execute(sql, params or ())
        return self.get_cursor.fetchall()

    """
    check if a table exist in the DB
    print the result
    
    Parameters
    ---------
    table : str
        a table name from the db
    """

    def check_if_table_exist(self, table):
        cmd = "SELECT count(name) FROM sqlite_master WHERE type='table' AND name = ?"
        val = self.query(cmd, (table,))
        if val[0] == 1:
            print(table + ' table exist')
        else:
            print(table + ' table does not exist')
        self.commit()

    """
    delete a row from a table based on a given value 
        
    Parameters
    ---------
    table : str
        table name from the db
    val : str
        name of column from the table
    id : str
        the value to delete
    """
    def delete_from_table(self, table, val, id):
        cmd = "DELETE FROM " + table + " WHERE " + val + "=?"
        self.executemany(cmd, id)
        self.commit()

    # ****************************PROTEOME********************************************
    """ 
    create a DB table for protein data:
    the table:
    Ref: short name -  SIRT6_human
    Uniprot_ID: protein id by uniprot
    Length: length of protein sequence
    Seq: sequence itself
    
    Parameters
    ---------
    table : str
        a table name from the db
    """

    def create_protein_table(self, table):
        print(self.DBname)
        cmd = "CREATE TABLE IF NOT EXISTS " + table + " (Ref VARCHAR(30), Uniprot_ID VARCHAR(30), Length, Seq)"
        # cmd = """CREATE TABLE human_proteome(Ref VARCHAR(30), Length, Seq); """
        self.query(cmd)
        # self.execute(cmd)
        self.commit()
        print("finished creating protein DB table " + table)
        # example:
        # createDB("lab_project.db", "human_proteome")

    """ 
    insert a protein csv file to a table in the db 
    
    Parameters
    ---------
    csv_file : str
        file name and location in which the data is from 
    table : str
        a table name from the db
    """

    def insert_protein_csv_to_db(self, csv_file, table):
        # connect to DB
        with open(csv_file, 'r') as fin:
            # csv.DictReader uses first line in file for column headings by default
            dr = csv.DictReader(fin)  # comma is default delimiter
            to_db = [(i['Ref'], i['Uniprot_ID'], i['Length'], i['Seq']) for i in dr]

        cmd = ("INSERT INTO " + table + " (Ref, Uniprot_ID, Length, Seq) VALUES (?, ?, ?, ?);")
        self.executemany(cmd, to_db)
        self.commit()
        print("finished uploading csv file into " + table)

    """
    find a protein in the db by name(Ref)
    
    Parameters
    ---------
    protein_name :  str
        the protein short name - reference
    table : str
        table name from the db
    """

    def find_in_protein_db_by_ref(self, protein_name, table):
        cmd = ("SELECT Ref, Length, Seq FROM " + table + " WHERE Ref=?")
        rows = self.query(cmd, (protein_name,))
        if rows is not None:
            # for row in rows:
            #     print(row)
            self.commit()
            return rows

    """
    find a protein in the proteim db by ID
     
    Parameters
    ---------
    ID :  str
        the protein uniprot id
    table : str
        table name from the db
    """

    def find_in_protein_db_by_id(self, ID, table):
        cmd = ("SELECT Ref, Length, Seq FROM " + table + " WHERE Uniprot_ID=?")
        rows = self.query(cmd, (ID,))
        if rows is not None:
            # for row in rows:
            #     print(row)
            self.commit()
            return rows

    # ****************************ACETYLOME********************************************
    """ 
    create a DB table for acetylom data
    
    the table:
    full name: the proteins full name
    Ref: short name
    Seq:sequence itself
    Pos: position of acetylated Lysine(k)
    
    Parameters
    ---------
    table : str
        table name from the db
    """

    def create_acetylom_table(self, table):
        cmd = "CREATE TABLE IF NOT EXISTS " + table + "(Key ,Name VARCHAR(30), Ref VARCHAR(30), Seq, Pos, Rel_pos, " \
                                                      "Uniprot_ID) "
        self.query(cmd)
        self.commit()
        print("finished creating protein DB table " + table)

    """ 
    insert a acetylom csv file to the db 
    
    Parameters
    ---------
    csv_file : str
        file name and location in which the data is from 
    table : str
        a table name from the db
    """

    def insert_acetylom_csv_to_db(self, csv_file, table):
        # connect to DB
        with open(csv_file, 'r') as fin:
            # csv.DictReader uses first line in file for column headings by default
            dr = csv.DictReader(fin)  # comma is default delimiter
            to_db = [(i['Key'], i['Name'], i['Ref'], i['Seq'], i['Pos'], i['Rel_pos'], i['Uniprot_ID']) for i in dr]

        cmd = (
                "INSERT INTO " + table + " (Key, Name, Ref, Seq, Pos, Rel_pos, Uniprot_ID) VALUES (?, ?, ?, ?, ?, ?, ?);")
        self.executemany(cmd, to_db)
        self.commit()
        print("finished uploading csv file into " + table)

    """
    find a acetylome in db by ref
    print the result
    
    Parameters
    ---------
    protein_name :  str
        the protein short name - reference
    table : str
        table name from the db
    """

    def find_in_acetylome_db_by_ref(self, protein_name, table):
        cmd = ("SELECT * FROM " + table + " WHERE Ref=?")
        rows = self.query_fetch_all(cmd, (protein_name,))
        for row in rows:
            print(row)
        self.commit()
        return rows

    """
    find a acetulome in db by ID
    print the result
       
    Parameters
    ---------
    ID :  str
        the protein uniprot id
    table : str
        table name from the db
    """

    def find_in_acetylome_db_by_id(self, ID, table):
        cmd = ("SELECT * FROM " + table + " WHERE Uniprot_ID=?")
        rows = self.query_fetch_all(cmd, (ID,))
        for row in rows:
            print(row)
        self.commit()
        return rows

    """
    get all values from a table
    
    Parameters
    ---------
    table : str
        table name from the db
    """

    def iterate_over_table(self, table):
        rows = self.query_fetch_all("SELECT * FROM " + table)
        return rows

    """
    update a value by table, protein reference, column and value
    
    Parameters
    ---------
    table : str
        table name from the db
    key : str
        key from the db table
    col : str
        column name from the table
    val : str
        value to insert to the table
    """

    def update_value_by_protein_key(self, table, key, col, val):
        try:
            cmd = ("UPDATE " + table + " SET " + col + " = " + str(val) + " WHERE Key=?")
            # print("command is: " + cmd)
            rows = self.query_fetch_all(cmd, (key,))
            self.commit()
            # return rows
        except Exception as e:
            print("exectipn during running the command + " + cmd + " is: " + e)

    # *******************************Homology**************************************
    """ 
    create a DB table for acetylom data
    
    the table: 
    full name: the proteins full name
    Ref: short name
    Seq:sequence itself
    Pos: position of acetylated Lysine(k)
    
    Parameters
    ---------
    table : str
        table name from the db
    """

    def create_homology_mouse_human_table(self, table):
        cmd = "CREATE TABLE IF NOT EXISTS " + table + "(Mouse_name ,Mouse_ID ,Human_name, Human_ID)"
        self.query(cmd)
        self.commit()
        print("finished creating protein DB table " + table)

    """ 
    insert a data csv file to the db 
    
    Parameters
    ---------
    csv_file : str
        file name and location in which the data is from 
    table : str
        a table name from the db
    """

    def insert_hom_csv_to_db(self, csv_file, table):
        # connect to DB
        with open(csv_file, 'r') as fin:
            # csv.DictReader uses first line in file for column headings by default
            dr = csv.DictReader(fin)  # comma is default delimiter
            to_db = [(i['Mouse_name'], i['Mouse_ID'], i['Human_name'], i['Human_ID']) for i in dr]

        cmd = ("INSERT INTO " + table + " (Mouse_name, Mouse_ID, Human_name, Human_ID) VALUES (?, ?, ?, ?);")
        self.executemany(cmd, to_db)
        self.commit()
        print("finished uploading csv file into " + table)

    """
    find row from the homology table by uniprot id
    
    Parameters
    ---------
    ID :  str
        the protein uniprot id
    table : str
        table name from the db
    """
    def find_in_hom_db_by_mouse_id(self, ID, table):
        cmd = ("SELECT * FROM " + table + " WHERE Mouse_ID=?")
        rows = self.query_fetch_all(cmd, (ID,))
        for row in rows:
            print(row)
        self.commit()
        return rows

        # *******************************Acetylome vs Animals Replacements**************************************

    """ 
    create a DB table for all replacements between a specific a.a data
    
    the table:
    aminal: aminal name
    animal id: id from the acetylome
    protein_id: id from protein in animal
    identities, positives: from blast
    GN: gene name

    Parameters
    ---------
    table : str
        table name from the db
    """

    def create_ac_rep_table(self, table):
        # row = ["Animal", "Acetylome_id","Pos", "Protein_id", "Identities", "Positives", "GN","Seq" "Acetylome_seq",
        #        "Protein_seq","Replacement"]
        cmd = "CREATE TABLE IF NOT EXISTS " + table + "(Animal ,Acetylome_id, Pos ,Protein_id, Identities, Positives, GN, " \
                                                      "Seq, Acetylome_seq , Protein_seq, Replacement)"
        self.query(cmd)
        self.commit()
        print("finished creating protein DB table " + table)

    """ 
    insert a data csv file to the db 
    
    Parameters
    ---------
    csv_file : str
        file name and location in which the data is from 
    table : str
        a table name from the db
    """

    def insert_ac_rep_csv_to_db(self, csv_file, table):
        # connect to DB
        with open(csv_file, 'r') as fin:
            # csv.DictReader uses first line in file for column headings by default
            dr = csv.DictReader(fin)  # comma is default delimiter
            to_db = [(i['Animal'], i['Acetylome_id'], i['Protein_id'], i['Identities'], i['Positives'],
                      i['GN'], i['Acetylome_seq'], i['Protein_seq'], i['Replacement']) for i in dr]

        cmd = ("INSERT INTO " + table + " (Animal ,Acetylome_id ,Protein_id, Identities, Positives, GN, "
                                        "Acetylome_seq , Protein_seq, Replacement) VALUES (?, ?, ?, ?, ?, ?, ?, "
                                        "?, ?);")
        self.executemany(cmd, to_db)
        self.commit()
        print("finished uploading csv file into " + table)
        
    """
    insert a row of data to the table
    
    Parameters
    ---------
    table : str
        the name of the table
    to_db : list
        a list of values
    """
    def insert_row_into_rep_table(self, table, to_db):
        cmd = ("INSERT INTO " + table + " (Animal ,Acetylome_id ,Protein_id, Pos, Identities, Positives, GN, "
                                        "Seq, Acetylome_seq , Protein_seq, Replacement) VALUES (?, ?, ?, ?, ?, ?, ?, "
                                        "?, ?, ?, ?);")
        self.executemany(cmd, to_db)
        self.commit()
        
        
    """
    insert all files from a dir to all_rep tables
    
    Parameters
    ---------
    dir_path : str
        path to a dir that contains all the csv files
    """

    def create_ac_rep_files_from_dir(self, dir_path):
        os.chdir(dir_path)
        all_files = dir_path + "*.csv"
        for file in glob.glob(all_files):
            print("in: " + file)
            name1 = file.split('/')[5]
            name = name1.split('_results')[0]
            table = name + '_all_ac_replacements'
            self.create_ac_rep_table(table)
            self.insert_ac_rep_csv_to_db(file, table)

    """
    find rows in table where id=id and replacement = rep
    
    Parameters
    ---------
    table : str
        table name from the db
    ID :  str
        the protein uniprot id
    rep : char
        the char that was found in the location
    """

    def find_ac_rep_in_db(self, table, ID, rep):
        cmd = ("SELECT * FROM " + table + " WHERE Acetylome_id =? AND Replacement =?")
        rows = self.query_fetch_all(cmd, (ID, rep,))
        # for row in rows:
        #     print(row)
        # self.commit()
        return rows

    # *******************************ALL Replacements**************************************
    """ 
    create a DB table for all replacements between a specific a.a data
    
    the table:
    aminal: aminal name
    animal id: uniprot id from the acetylome
    protein_id: uniprot id from protein in animal
    identities, positives: from blast
    GN: gene name
    Query_seq: peptide sequence from the acetylome
    Sbjct_seq: protein sequence to compare
    Replacements: a char that was found in the location
    Error: errors that are found
    
    Parameters
    ---------
    table : str
        table name from the db
    """

    def create_all_rep_table(self, table):
        # row = ["Animal", "Animal_id", "Protein_id", "Identities", "Positives", "GN", "Query_seq", "Sbjct_seq",
        #        "Replacements", "Error"]
        cmd = "CREATE TABLE IF NOT EXISTS " + table + "(Animal ,Animal_id ,Protein_id, Identities, Positives, GN, " \
                                                      "Seq, Query_seq , Sbjct_seq, Replacements)"
        self.query(cmd)
        self.commit()
        print("finished creating protein DB table " + table)

    """ 
    insert a data csv file to the db 
    
    Parameters
    ---------
    csv_file : str
        file name and location in which the data is from 
    table : str
        a table name from the db
    """

    def insert_all_rep_csv_to_db(self, csv_file, table):
        # connect to DB
        with open(csv_file, 'r') as fin:
            # csv.DictReader uses first line in file for column headings by default
            dr = csv.DictReader(fin)  # comma is default delimiter
            to_db = [(i['Animal'], i['Aminal_id'], i['Protein_id'], i['Identities'], i['Positives'],
                      i['GN'], i['Query_seq'], i['Sbjct_seq'], i['Replacements'], i['Error']) for i in dr]

        cmd = ("INSERT INTO " + table + " (Animal ,Aminal_id ,Protein_id, Identities, Positives, GN, "
                                        "Query_seq , Sbjct_seq, Replacements, Error ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, "
                                        "?, ?);")
        self.executemany(cmd, to_db)
        self.commit()
        print("finished uploading csv file into " + table)

    """
    insert a row of data to the table

    Parameters
    ---------
    table : str
        the name of the table
    to_db : list
        a list of values
    """

    def insert_row_into_all_rep_table(self, table, to_db):
        cmd = ("INSERT INTO " + table + " (Animal ,Animal_id ,Protein_id, Identities, Positives, GN, "
                                        "Seq, Query_seq , Sbjct_seq, Replacements ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, "
                                        "?, ?);")

        self.executemany(cmd, to_db)
        self.commit()
    """
    insert all files from a dir to all_rep tables
    
    Parameters
    ---------
    dir_path : str
        path to a dir that contains all the csv files
    """

    def create_all_rep_files_from_dir(self, dir_path, rep):
        os.chdir(dir_path)
        all_files = dir_path + "*.csv"
        for file in glob.glob(all_files):
            print("in: " + file)
            name1 = file.split('/')[5]
            name = name1.split('_all_rep_results')[0]
            table = name + '_' + rep
            self.create_all_rep_table(table)
            self.insert_all_rep_csv_to_db(file, table)




# **************************************************SIRT6 acetylation*************************************************
    """ 
    create a DB table for all replacements between a specific a.a data

    the table:
    Uniprot_id: uniprot id from protein in animal
    Pos: position of the ac site
    GN: gene name
    Seq: peptide sequence from the acetylome
    
    Parameters
    ---------
    table : str
        table name from the db
    """

    def create_sirt6_table(self, table):
        cmd = "CREATE TABLE IF NOT EXISTS " + table + "(Uniprot_ID ,Pos ,GN, Seq)"
        self.query(cmd)
        self.commit()
        print("finished creating protein DB table " + table)

    """ 
    insert a data csv file to the db 

    Parameters
    ---------
    csv_file : str
        file name and location in which the data is from 
    table : str
        a table name from the db
    """

    def insert_sirt6_csv_to_db(self, csv_file, table):
        # connect to DB
        with open(csv_file, 'r') as fin:
            # csv.DictReader uses first line in file for column headings by default
            dr = csv.DictReader(fin)  # comma is default delimiter
            to_db = [(i['Uniprot_ID'], i['Pos'], i['GN'], i['Seq']) for i in dr]

        cmd = ("INSERT INTO " + table + " (Uniprot_ID ,Pos ,GN, Seq) VALUES (?, ?, ?, ?);")
        self.executemany(cmd, to_db)
        self.commit()
        print("finished uploading csv file into " + table)

    """
    find a acetulome in db by ID
    print the result

    Parameters
    ---------
    ID :  str
        the protein uniprot id
    table : str
        table name from the db
    """

    def find_in_sirt6_db_by_ID(self, id, table):
        cmd = ("SELECT * FROM " + table + " WHERE Uniprot_ID=?")
        rows = self.query_fetch_all(cmd, (id,))
        for row in rows:
            print(row)
        self.commit()
        return rows


    """
    find a acetulome in db by gene name
    print the result

    Parameters
    ---------
    ID :  str
        the protein gene
    table : str
        table name from the db
    """

    def find_in_sirt6_db_by_GN(self, gn, table):
        cmd = ("SELECT * FROM " + table + " WHERE GN=?")
        rows = self.query_fetch_all(cmd, (gn,))
        for row in rows:
            print(row)
        self.commit()
        return rows

    # **************************************************Animals LifeSpans*************************************************
    """ 
    create a DB table for all replacements between a specific a.a data

    the table:
    Uniprot_id: uniprot id from protein in animal
    Pos: position of the ac site
    GN: gene name
    Seq: peptide sequence from the acetylome

    Parameters
    ---------
    table : str
        table name from the db
    """

    def create_species_life_spans_table(self, table):
        cmd = "CREATE TABLE IF NOT EXISTS " + table + "(species, lifespan)"
        self.query(cmd)
        self.commit()
        print("finished creating protein DB table " + table)

    """ 
    insert a data csv file to the db 

    Parameters
    ---------
    csv_file : str
        file name and location in which the data is from 
    table : str
        a table name from the db
    """

    def insert_species_life_spans_csv_to_db(self, csv_file, table):
        # connect to DB
        with open(csv_file, 'r') as fin:
            # csv.DictReader uses first line in file for column headings by default
            dr = csv.DictReader(fin)  # comma is default delimiter
            to_db = [(i['species'], i['lifespan']) for i in dr]

        cmd = ("INSERT INTO " + table + " (species, lifespan) VALUES (?, ?);")
        self.executemany(cmd, to_db)
        self.commit()
        print("finished uploading csv file into " + table)
        
    


