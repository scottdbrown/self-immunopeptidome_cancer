'''
Lookup HLA genotypes
Using multiprocessing, query database to get sets of peptides presented by random HLA genotypes.

Date: February 3, 2017
@author: sbrown
'''

## Import Libraries
import sys
import argparse
import os
import sqlite3
import time
import multiprocessing as mp

DEBUG = False
VERB = False

maxBufferSize = 100

SENTINEL = None


def lookupGenotype(dbp, q, out_q, i):
    resHolder = []
    numInHolder = 0

    ## get HLA ids
    if VERB: print("Getting HLA ids in slave {}...".format(i+1))
    hlaID = {}
    db = sqlite3.connect(dbp)
    res = db.execute("SELECT * FROM hla")
    for hid, allele in res:
        hlaID[allele] = hid
    
    ## Allows return of non-tuples:
    db.row_factory = lambda cursor, row: row[0]

    ## Process genotype queue
    if VERB: print("Processing queue in slave {}...".format(i+1))
    while True:
        hla_geno_list = q.get()
        if hla_geno_list is SENTINEL:
            break
        sid, a1, a2, b1, b2, c1, c2 = hla_geno_list

        genotype = "{}_{}_{}_{}_{}_{}".format(a1, a2, b1, b2, c1, c2)

        if DEBUG: print("Looking up genotype {}...".format(genotype))

        ## only add to set if HLA does not end in "N"
        hlas = set()
        for h in [a1, a2, b1, b2, c1, c2]:
            if not h.endswith("N"):
                hlas.add(hlaID[h])

        #hlas = set([hlaID[a1], hlaID[a2], hlaID[b1], hlaID[b2], hlaID[c1], hlaID[c2]])

        hlaqry = "({})".format(",".join([str(x) for x in hlas]))

        ## Test using sqlite DISTINCT vs. pulling all and doing set()
        numPepBind = len(set(db.execute("SELECT pep_id FROM binders WHERE hla_id IN {}".format(hlaqry)).fetchall()))
        ## this is slower:
        #setOfBindingPeps = db.execute("SELECT DISTINCT(pep_id) FROM binders WHERE hla_id IN {}".format(hlaqry)).fetchall()

        #resHolder.append([sid, genotype, numPepBind])
        resHolder.append([sid, numPepBind])
        numInHolder += 1
        if numInHolder == maxBufferSize:
            out_q.put(resHolder)
            resHolder = []
            numInHolder = 0
    out_q.put(resHolder)
    resHolder = []
    numInHolder = 0

    db.close()

    print("\nGenotype lookup complete in slave {}...".format(i+1))


if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Construct and Lookup Random HLA genotypes")
    parser.add_argument("hla_genotype_list", help = "File with HLA genotypes to lookup", type = str)
    parser.add_argument("database_file", help = "Database file to read", type = str)
    parser.add_argument("outFile", help = "file to write output to", type = str)
    parser.add_argument("maxNumberProcesses", help = "Maximum number of processes to start", type = int)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    ## Load HLA-id mapping from database
    print("Getting HLA ids from database...")
    timecheck = time.time()

    hlaID = {}
    db = sqlite3.connect(args.database_file)
    res = db.execute("SELECT * FROM hla")
    for hid, allele in res:
        hlaID[allele] = hid
    db.close()

    print("Took {:.2f} seconds...".format(time.time() - timecheck))


    ## load all HLA genotypes from list
    geno = []
    
    print("Reading in HLA file...")
    timecheck = time.time()
    for line in open(args.hla_genotype_list, "r"):
        ## expected format: sub_id \t hlaa_hlaa_hlab_hlab_hlac_hlac \n
        sid, gt = line.rstrip().split("\t")
        
        ## initialize hla_peptides for each hla
        UNKNOWN_HLA = False
        for hla in gt.split("_"):
            #if DEBUG: print("hla: {}\n".format(hla))
            ## check for HLA ending with N, meaning not expressed. ignore this for now, handle it later.
            if not hla.endswith("N") and hla not in hlaID:
                print("Unknown HLA allele: {} - dropping subject {}".format(hla, sid))
                UNKNOWN_HLA = True


        if not UNKNOWN_HLA:
            geno.append([sid,gt])

        

    print("Took {:.2f} seconds...".format(time.time() - timecheck))

    ## Spin up processes.
    print("Making genotype processing slaves...")

    try:
        ## Make output queue
        out_qs = []

        ## Start processes.
        procs = []
        pqs = []
        for i in range(args.maxNumberProcesses):
            ## make queue to hold output
            oq = mp.Queue()
            ## make queue to hold genotypes to process
            q = mp.Queue()
            ## make process to process genotype
            p = mp.Process(target=lookupGenotype, args=(args.database_file, q, oq, i+1))
            ## add process to list of processes.
            procs.append(p)
            ## start the process
            p.start()
            ## add queue to list of queues.
            pqs.append(q)
            out_qs.append(oq)

        ## Build genotypes and send to queues
        print("Submitting genotypes...")
        timecheck = time.time()

        proc_ind = 0
        geno_num = 0
        
        while geno_num < len(geno):
            if geno_num % 100 == 0:
                print("{:,} genotypes submitted to slaves...".format(geno_num), end="\r")

            sid = geno[geno_num][0]
            a1, a2, b1, b2, c1, c2 = geno[geno_num][1].split("_")

            ## add this genotype to a queue.
            pqs[proc_ind].put([sid,a1,a2,b1,b2,c1,c2])
            ## rotate through to next proc
            proc_ind += 1
            if proc_ind >= args.maxNumberProcesses:
                proc_ind = 0

            geno_num += 1

        ## done all genotypes, now send the end value
        for pq in pqs:
            pq.put(SENTINEL)

        print("\nTook {:.2f} seconds...".format(time.time() - timecheck))


        print("All genotypes submitted, beginning to write in {:,} line blocks...".format(maxBufferSize))


        timecheck = time.time()
        out = open(args.outFile, "w")
        #out.write("sample\tgenotype\tnumBinders\n")
        out.write("sample\tnumBinders\n")

        resultBuffer = []
        #numInBuffer = 0
        numBlocks = 0
        blockTime = time.time()


        numGeno = 0
        proc_ind = 0
        while numGeno < len(geno):
            ## cycle through output queues and slurp results
            if not out_qs[proc_ind].empty():
                res = out_qs[proc_ind].get()
                for trip in res:
                    resultBuffer.append("\t".join(map(str,trip)))
                    #numInBuffer += 1
                    numGeno += 1

            proc_ind += 1
            if proc_ind >= args.maxNumberProcesses:
                proc_ind = 0


            if len(resultBuffer) > 0:
                numBlocks += 1
                print("\nBlock {} took {:.2f} seconds...writing...".format(numBlocks, time.time() - blockTime))
                print("{} genotypes written.".format(numGeno))
                out.write("\n".join(resultBuffer))
                out.write("\n")
                resultBuffer = []
                #numInBuffer = 0
                blockTime = time.time()

        if len(resultBuffer) > 0:
            print("\nWriting final block...")
            out.write("\n".join(resultBuffer))
            out.write("\n")
            out.close()

        print("Took {:.2f} seconds...".format(time.time() - timecheck))

        ## make sure all procs have finished.
        print("Making sure all genotype slaves have shutdown...")
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

    print("done.")