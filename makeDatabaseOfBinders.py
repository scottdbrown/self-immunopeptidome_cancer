'''
Make Database of Binders
Parse through NetMHCpan results and build an SQlite3 database of the binding peptides.

Date: February 1, 2017
@author: sbrown
'''

## Import Libraries
import sys
import argparse
import os
import sqlite3
import multiprocessing as mp
import time
import scandir
import traceback

DEBUG = False
VERB = False

IC50_THRESH = 500
maxBufferSize = 100000

SENTINEL = None

hlaID = {}
pepID = {}


def connectAndWriteDB(database, data, query):
    db = sqlite3.connect(database)
    db.executemany(query, data)
    db.commit()
    db.close()

def processPeptideFile(q, oq, hids, pids, i):
    numProcessed = 0
    resHolder = []
    numInHolder = 0
    while True:
        dat = q.get()
        #print("DEBUG q.get: {}".format(dat))
        if dat is SENTINEL:
            break
        hla, pepLen, scoreFile, refFile = dat
        numProcessed += 1

        for lineScore, lineRef in zip(open(scoreFile, "r"), open(refFile, "r")):
            if float(lineScore.rstrip()) <= IC50_THRESH:
                resHolder.append((hids[hla], pids[lineRef.rstrip()], lineScore.rstrip()))
                numInHolder += 1

        if numInHolder > maxBufferSize:
            oq.put(resHolder)
            resHolder = []
            numInHolder = 0

    oq.put(resHolder)
    resHolder = []
    numInHolder = 0

    oq.put(SENTINEL)
    print("\nProcessing complete in slave {}...".format(i+1))


if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Make Database of Binders")
    parser.add_argument("species_code", help = "embl species oscode", type = str)
    parser.add_argument("root_dir", help = "Root directory for species analysis", type = str)
    parser.add_argument("hla_list", help = "File with HLA alleles to use (each only once)", type = str)
    parser.add_argument("database_file", help = "Database file to create", type = str)
    parser.add_argument("maxNumberProcesses", help = "Maximum number of processes to start", type = int)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB


    ## Check that database file does not already exist.
    if os.path.exists(args.database_file):
        print("Database file {} already exists. Please provide a non-existant database.".format(args.database_file))
        sys.exit()


    ## Make peptide mapping

    db = sqlite3.connect(args.database_file)
    db.execute("CREATE TABLE peptide(id INT, sequence TEXT)")
    db.commit()
    db.close()

    print("Processing all peptide sequences...")
    timecheck = time.time()
    
    pep_i = 1
    pep_toWrite = []

    for f in os.listdir(args.root_dir):
        if f.endswith("_peptides.txt"):
            for line in open(os.path.join(args.root_dir,f), "r"):
                pep = line.rstrip()

                pepID[pep] = pep_i
                pep_toWrite.append((pep_i, pep))
                pep_i += 1

                if len(pepID) % maxBufferSize == 0:
                    print("                                                 ", end="\r")
                    print("{:,} peptides processed...".format(len(pepID)), end="", flush=True)
                    print("Writing...", end="", flush=True)
                    qry = "INSERT INTO peptide(id, sequence) VALUES (?, ?)"
                    connectAndWriteDB(args.database_file, pep_toWrite, qry)
                    pep_toWrite = []
                    print("done.", end="\r", flush=True)

    ## clear the writing buffer
    print("                                                 ", end="\r")
    print("{:,} peptides processed...".format(len(pepID)), end="", flush=True)
    print("Writing...", end="", flush=True)
    qry = "INSERT INTO peptide(id, sequence) VALUES (?, ?)"
    connectAndWriteDB(args.database_file, pep_toWrite, qry)
    pep_toWrite = []
    print("done.", end="\r", flush=True)


    print("\nTook {:.2f} seconds...".format(time.time() - timecheck))

    ## Make peptide table

    ## Extract all HLA, make table
    print("Reading in HLA file...")
    timecheck = time.time()

    hla_i = 1
    hla_toWrite = []

    for line in open(args.hla_list, "r"):
        hla = line.rstrip().replace(":","-")
        
        hlaID[hla] = hla_i
        hla_toWrite.append((hla_i, hla))
        hla_i += 1

    print("Took {:.2f} seconds...".format(time.time() - timecheck))

    print("Writing HLA to database...")
    timecheck = time.time()

    db = sqlite3.connect(args.database_file)
    db.execute("CREATE TABLE hla(id INT, allele TEXT)")
    db.executemany("INSERT INTO hla(id, allele) VALUES (?, ?)", hla_toWrite)
    db.commit()
    db.close()

    print("Took {:.2f} seconds...".format(time.time() - timecheck))

    ## Distribute result file parsing
    print("Finding and parsing all results files...")
    print("Will print files being submitted to slave 1...")
    timecheck = time.time()

    db = sqlite3.connect(args.database_file)
    db.execute("CREATE TABLE binders(hla_id INT, pep_id INT, ic50 REAL)")
    db.commit()
    db.close()

    ## read in all peptides, keep those meeting threshold.
    ## most species will have multiple contig files/runs for each index.
    try:
        procs = []
        pqs = []
        out_qs = []
        for i in range(args.maxNumberProcesses):
            q = mp.Queue()
            oq = mp.Queue()
            p = mp.Process(target=processPeptideFile, args=(q, oq, hlaID, pepID, i))
            procs.append(p)
            p.start()
            pqs.append(q)
            out_qs.append(oq)


        ## assign HLA alleles to each process.


        proc_ind = 0
        for root, dirs, files in scandir.walk(args.root_dir):
            for f in files:
                if f.endswith(".pMHC.parsed"):
                    ind = "_".join(f.split(".")[0].split("_")[0:3])
                    hla = f.split(".")[0].split("_")[1]

                    if proc_ind == 0: print("Submitting {} to slave {}...     ".format(f, proc_ind+1), end="\r")
                    
                    pepLen = int(f.split(".")[0].split("_")[2])
                    contigFileNum = int(f.split(".")[0].split("_")[3])

                    pepScoreFile = os.path.join(root,f)
                    pepRefFile = os.path.join(args.root_dir, "prot{}_{}_{}_peptides.txt".format(pepLen,contigFileNum,args.species_code))

                    ## add to queue
                    pqs[proc_ind].put([hla, pepLen, pepScoreFile, pepRefFile])

                    proc_ind += 1
                    if proc_ind >= args.maxNumberProcesses:
                        proc_ind = 0
                    

        ## Send END signal to each queue
        for pq in pqs:
            pq.put(SENTINEL)

        print("\nTook {:.2f} seconds...".format(time.time() - timecheck))

        print("Writing to database as files are processed...")
        timecheck = time.time()
        
        ## Collect chunks from Processes in out_qs, write to database
        proc_ind = 0
        completed = [False for x in range(args.maxNumberProcesses)]
        
        while not all(completed):
            ## cycle through output queues and slurp results
            if not out_qs[proc_ind].empty():
                res = out_qs[proc_ind].get()
                if res is SENTINEL:
                    ## that proc is done.
                    completed[proc_ind] = True
                else:
                    ## write to database
                    print("                                          ", end="\r", flush=True)
                    print("Writing chunk from slave {}...".format(proc_ind+1), end="\r", flush=True)
                    ## If res is a list of lists, .get() only grabs the first.
                    qry = "INSERT INTO binders(hla_id, pep_id, ic50) VALUES (?, ?, ?)"
                    connectAndWriteDB(args.database_file, res, qry)
                    print("Writing chunk from slave {}...done. ".format(proc_ind+1), end="\r", flush=True)
            proc_ind += 1
            if proc_ind >= args.maxNumberProcesses:
                proc_ind = 0

        ## make sure all procs have finished.
        print("\nShutting down all slaves...")
        for p in procs:
            p.join()

    except KeyboardInterrupt:
        print("Keyboard Interruption: ".format(sys.exc_info()[0]))
        print(traceback.format_exc())
        print("Please wait while processes close...")
        for q in pqs:
            while not q.empty():
                foo = q.get()
            q.put(SENTINEL)
        for p in procs:
            p.join()

        sys.exit()

    print("Done.")