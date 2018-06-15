'''
Lookup Mutation Read Support
Using multiprocessing, check for read support from bam files for mutations

Date: May 7, 2018
@author: sbrown
'''

## Import Libraries
import sys
import argparse
import os
import sqlite3
import time
import multiprocessing as mp
import subprocess
import traceback

DEBUG = False
VERB = False

maxBufferSize = 100

MAX_ATTEMPTS = 5

SENTINEL = None


SAM_BIN = "/gsc/software/linux-x86_64-centos5/samtools-0.1.8/samtools"


def baseComplement(base):
    if base == "A":
        return "T"
    elif base == "T":
        return "A"
    elif base == "G":
        return "C"
    elif base == "C":
        return "G"
    else:
        return "X"

def lookup(in_q, out_q, i):
    resHolder = []
    numInHolder = 0

    ## Process queue
    if VERB: print("Processing queue in slave {}...".format(i+1))
    while True:
        ## do lookup
        mut_info = in_q.get()
        if mut_info is SENTINEL:
            break
        
        barcode, chrom, pos, mut, wild, bam = mut_info

        mutCount = 0
        wildCount = 0
        otherCount = 0

        ## determine if chr needs to be appended to chromosome (inconsistent file naming)
        cmd = "{} view {} | head -1".format(SAM_BIN, bam)
        VALID_COM = False
        while not VALID_COM:
            call = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            (res, err) = call.communicate()
            attempts = 0
            if err.decode('ascii') == "":
                VALID_COM = True
            elif attempts < MAX_ATTEMPTS:
                print("ERROR IN RUNNING COMMAND: {}\nError: {}".format(cmd, str(err.decode('ascii'))))
                print("Waiting 5 seconds to try again...")
                attempts += 1
                time.sleep(5)
            else:
                sys.exit("Unable to run command.")

        reads = res.decode('ascii')

        chr_prefix = ""
        if reads.split("\t")[2].startswith("chr"):
            chr_prefix = "chr"

        ## build cmd
        cmd = "{} mpileup -r {}{}:{}-{} {} | cut -f 5".format(SAM_BIN, chr_prefix, chrom, pos, pos, bam)

        #print("cmd: {}".format(cmd))
        VALID_COM = False
        while not VALID_COM:
            call = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            (res, err) = call.communicate()
            if err.decode('ascii') == "":
                VALID_COM = True
            else:
                print("ERROR IN RUNNING COMMAND: {}\nError: {}".format(cmd, str(err.decode('ascii'))))
                print("Waiting 2 seconds to try again...")
                time.sleep(2)
        
        ## capture the output, parse to check if read support variant, write to file       
        bases = res.decode('ascii')
        bases = bases.rstrip()

        basecounts = {"A":0, "T":0, "C":0, "G":0}
        totalcounts = 0
        if len(bases) > 0:
            for b in bases:
                b = b.upper()
                if b in basecounts:
                    basecounts[b] += 1
                    totalcounts += 1
        

        wildCount = basecounts[wild]
        mutCount = basecounts[mut]
        otherCount = totalcounts - wildCount - mutCount


        #resHolder.append([sid, genotype, numPepBind])
        resHolder.append([barcode, chrom, pos, wild, mut, wildCount, mutCount, otherCount])
        numInHolder += 1
        if numInHolder == maxBufferSize:
            out_q.put(resHolder)
            resHolder = []
            numInHolder = 0
    out_q.put(resHolder)
    resHolder = []
    numInHolder = 0

    print("\nLookup complete in slave {}...".format(i+1))


if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Check for mutation read support")
    parser.add_argument("list_of_muts", help = "File with mutations and bams to lookup", type = str)
    parser.add_argument("outFile", help = "file to write output to", type = str)
    parser.add_argument("maxNumberProcesses", help = "Maximum number of processes to start", type = int)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB


    ## Spin up processes.
    print("Making processing slaves...")

    try:
        ## Make output queue
        out_qs = []

        ## Start processes.
        procs = []
        pqs = []
        for i in range(args.maxNumberProcesses):
            ## make queue to hold output
            oq = mp.Queue()
            ## make queue to hold mutations to process
            q = mp.Queue()
            ## make process to process genotype
            p = mp.Process(target=lookup, args=(q, oq, i+1))
            ## add process to list of processes.
            procs.append(p)
            ## start the process
            p.start()
            ## add queue to list of queues.
            pqs.append(q)
            out_qs.append(oq)

        ## Build mutations and send to queues
        print("Submitting mutations...")
        
        proc_ind = 0
        mut_num = 0

        HEADER = True
        for line in open(args.list_of_muts, "r"):
            if HEADER:
                HEADER = False
            else:
                line = line.rstrip().split("\t")
                barcode = line[1]
                chrom = line[3]
                pos = int(line[4])
                strand = line[6]
                hgvsc = line[7]
                bam = line[8]


                hgvsc = hgvsc.split(">")
                wild = hgvsc[0][-1]
                mut = hgvsc[1]

                if strand == "-1":
                    ## complement the nucleotides
                    wild = baseComplement(wild)
                    mut = baseComplement(mut)

                if mut_num % 1000 == 0:
                    print("{:,} mutations submitted to slaves...".format(mut_num), end="\r")

                
                ## add this mutation to a queue.
                pqs[proc_ind].put([barcode, chrom, pos, mut, wild, bam])
                ## rotate through to next proc
                proc_ind += 1
                if proc_ind >= args.maxNumberProcesses:
                    proc_ind = 0

                mut_num += 1

        ## done all mutations, now send the end value
        for pq in pqs:
            pq.put(SENTINEL)


        print("All mutations submitted, beginning to parse in {:,} line blocks...".format(maxBufferSize))


        toWrite = {}

        numMut = 0
        proc_ind = 0
        while numMut < mut_num:
            ## cycle through output queues and slurp results
            if not out_qs[proc_ind].empty():
                res = out_qs[proc_ind].get()
                for trip in res:
                    barcode, chrom, pos, wild, mut, wildCount, mutCount, otherCount = trip
                    if (barcode, chrom, pos, wild, mut) not in toWrite:
                        toWrite[(barcode, chrom, pos, wild, mut)] = [wildCount, mutCount, otherCount]
                    else:
                        ## already have entry for this mutation from sample with two bam files (>1 lane, unmerged)
                        toWrite[(barcode, chrom, pos, wild, mut)][0] += wildCount
                        toWrite[(barcode, chrom, pos, wild, mut)][1] += mutCount
                        toWrite[(barcode, chrom, pos, wild, mut)][2] += otherCount
                    numMut += 1

            proc_ind += 1
            if proc_ind >= args.maxNumberProcesses:
                proc_ind = 0


        ## make sure all procs have finished.
        print("Making sure all lookup slaves have shutdown...")
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



    out = open(args.outFile, "w")
    #out.write("sample\tgenotype\tnumBinders\n")
    out.write("barcode\tchrom\tpos\twild\tmut\twildCount\tmutCount\totherCount\n")

    resultBuffer = []
    #numInBuffer = 0
    numBlocks = 0
    numMut = 0

    for k in toWrite:
        resultBuffer.append("\t".join(map(str, list(k)+toWrite[k])))
        numMut += 1
    

        if len(resultBuffer) > 0 & len(resultBuffer)%100==0:
            numBlocks += 1
            print("{} mutations written.".format(numMut))
            out.write("\n".join(resultBuffer))
            out.write("\n")
            resultBuffer = []
            #numInBuffer = 0

    if len(resultBuffer) > 0:
        print("\nWriting final block...")
        out.write("\n".join(resultBuffer))
        out.write("\n")
        out.close()


    print("done.")
