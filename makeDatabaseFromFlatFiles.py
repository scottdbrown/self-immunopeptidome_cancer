'''
Make Database from Flat Files
Using the flat files of the lists of peptides presented by each HLA, create sqlite3 database

Date: July 6, 2018
@author: sbrown
'''

## Import Libraries
import sys
import argparse
import os
import time
import sqlite3

DEBUG = False
VERB = False



def connectAndWriteDB(database, data, query):
    db = sqlite3.connect(database)
    db.executemany(query, data)
    db.commit()
    db.close()

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    ENDC = '\033[0m'

def log_print(msg_type, msg):
    print("{}{}{} {}[{}]{}: {}".format(bcolors.BOLD, msg_type, bcolors.ENDC, bcolors.OKBLUE, time.strftime("%Y/%m/%d %T"), bcolors.ENDC, msg))


if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Make Database from Flat Files")
    ## add_argument("name", "(names)", metavar="exampleOfValue - best for optional", type=int, nargs="+", choices=[allowed,values], dest="nameOfVariableInArgsToSaveAs")
    parser.add_argument("flatfile_dir", help = "Directory containing flat files", type = str)
    parser.add_argument("new_database", help = "Database to create", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    print("=======================================================")
    print("Python version: {}".format(sys.version))
    print("Server: {}".format(os.uname()[1]))
    print("Current directory: {}".format(os.getcwd()))
    print("Command: {}".format(" ".join(sys.argv)))
    print("Time: {}".format(time.strftime("%Y/%m/%d %T")))
    print("=======================================================\n")


    ## Check that database file does not already exist.
    if os.path.exists(args.new_database):
        print("Database file {} already exists. Please provide a non-existant database.".format(args.database_file))
        sys.exit()


    ## Get list of HLA
    hlas = os.listdir(args.flatfile_dir)
    hla_ids = {}
    i = 1
    hla_toWrite = []
    for hla in hlas:
        ## cleave off trailing ".txt"
        hla = hla.split(".")[0]
        ## add to id dictionary
        hla_ids[hla] = i
        hla_toWrite.append((i, hla))
        i += 1

    ## write HLA to database
    log_print("STATUS", "Writing HLA data to database...")
    db = sqlite3.connect(args.new_database)
    db.execute("CREATE TABLE hla(id INT, allele TEXT)")
    db.executemany("INSERT INTO hla(id, allele) VALUES (?, ?)", hla_toWrite)
    db.commit()
    
    ## create peptide table and binders table
    db.execute("CREATE TABLE peptide(id INT, sequence TEXT)")
    db.execute("CREATE TABLE binders(hla_id INT, pep_id INT)")
    db.commit()
    db.close()

    ## Parse through files
    log_print("STATUS", "Parsing through files...")
    pep_ids = {}
    i = 1
    pep_toWrite = []
    bind_toWrite = []
    for hla_file in hlas:
        hla = hla_file.split(".")[0]
        ## open the file
        for line in open(os.path.join(args.flatfile_dir, hla_file), "r"):
            pep = line.rstrip()
            ## check if we have seen this peptide before
            if pep not in pep_ids:
                pep_ids[pep] = i
                pep_toWrite.append((i, pep))
                i += 1
            bind_toWrite.append((hla_ids[hla], pep_ids[pep]))

    ## write to database
    log_print("STATUS", "Writing peptides and binders to database...")
    db = sqlite3.connect(args.new_database)
    db.executemany("INSERT INTO peptide(id, sequence) VALUES (?, ?)", pep_toWrite)
    db.executemany("INSERT INTO binders(hla_id, pep_id) VALUES (?, ?)", bind_toWrite)
    db.commit()
    db.close()

    print("\n=================================================")
    print("To improve performance of database, connect to database using '$ sqlite3 {}'".format(args.new_database))
    print("Then run indexing by typing:")
    print("> CREATE INDEX peptide_ind ON peptide(id);")
    print("> CREATE INDEX peptide_seq ON peptide(sequence);")
    print("> CREATE INDEX binder_hla_ind ON binders(hla_id);")
    print("> CREATE INDEX binder_pep_ind ON binders(pep_id);")
    print("=================================================\n")
    log_print("STATUS", "Done.")