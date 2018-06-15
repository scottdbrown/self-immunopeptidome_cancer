'''
Tally Parsed Data
Get summary stats for each HLA regarding number of peptides predicted to bind at different levels

Date: September 23, 2016
@author: sbrown

Edited January 11, 2017:
    - Write incompleted job names to a file
    - Add time parsing data too
Edited January 15, 2017:
    - Use multiprocessing to speed up.
'''

## Import Libraries
import sys
import argparse
import os
import re
import datetime
import multiprocessing as mp
import time

DEBUG = False
VERB = False

SENTINEL = None

maxBufferSize = 10000
maxSlaveSize = 100

def processJob(q, out_q, i):
    resHolder = []
    ## will hold [[jobname, isClean], [jobname,hla,pepLen,date,dt,duration], [hla, {}]]
    numInHolder = 0

    while True:
        dat = q.get()
        if DEBUG: print("\ndat in slave {} is: {}".format(i, dat))
        if dat is SENTINEL:
            break

        jobname, direc = dat

        timeRes = []
        hlaRes = []
        
        errFile = None
        stdFile = None
        resFile = None

        if os.path.exists(direc):

            for file in os.listdir(direc):
                if re.match(".*\.sh\.e.*", file):
                    errFile = file
                elif re.match(".*\.sh\.o.*", file):
                    stdFile = file
                elif re.match(".*\.parsed", file):
                    resFile = file

            jobIsClean = True

        ## check that shell err file empty
        if errFile:
            for line in open(os.path.join(direc, errFile), "r"):
                if line:
                    ## file should be empty, there is an error.
                    jobIsClean = False
        else:
            ## job did not run
            jobIsClean = False

        ## check the time
        if jobIsClean:
            if stdFile:
                times = [0,0]
                lineNum = 0
                for line in open(os.path.join(direc, stdFile), "r"):
                    times[lineNum] = datetime.datetime.strptime(line.rstrip(), "%a %b %d %H:%M:%S %Z %Y")
                    lineNum += 1
                #print(times)

                if times[0] == 0 or times[1] == 0:
                    if DEBUG: print("Job {} did not complete.".format(f.split(".")[0]))
                    jobIsClean = False

                else:
                    duration = times[1] - times[0]
                    duration = duration.total_seconds()

                    date = str(times[0].date())
                    dt = str(times[0])
                    hla = jobname.split("_")[1]
                    pepLen = jobname.split("_")[2]

                    timeRes = [jobname, hla, pepLen, date, dt, duration]
            else:
                jobIsClean = False

        ## parse results
        if jobIsClean:
            ## resFile like CHLTR_HLA-A01-01_10_1.pMHC.parsed
            (species, hla, pepLen, fnum) = resFile.split(".")[0].split("_")
            pepLen = int(pepLen)
            hlaResDict = {50: [0 for x in range(0,4)], 100: [0 for x in range(0,4)], 500: [0 for x in range(0,4)]}
            for line in open(os.path.join(direc,resFile), "r"):
                line = float(line.rstrip())
                for cutoff in hlaResDict:
                    if line <= cutoff:
                        hlaResDict[cutoff][pepLen - 8] += 1
                hlaRes = [hla, hlaResDict]

        resHolder.append([[jobname,jobIsClean], hlaRes, timeRes])
        
        numInHolder += 1
        if numInHolder == maxSlaveSize:
            out_q.put(resHolder)
            resHolder = []
            numInHolder = 0
    out_q.put(resHolder)
    resHolder = []
    numInHolder = 0

    print("\nProcessing in slave {} complete...".format(i))



if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Tally Parsed Data")
    ## add_argument("name", "(names)", metavar="exampleOfValue - best for optional", type=int, nargs="+", choices=[allowed,values], dest="nameOfVariableInArgsToSaveAs")
    parser.add_argument("scriptReferenceFile", help = "Script file used by clusterTAS, containing job names", type = str)
    parser.add_argument("resultsDir", help = "Directory with parsed result files", type = str)
    parser.add_argument("resultsOutputFile", help = "File to write results output to", type = str)
    parser.add_argument("timeOutputFile", help = "File to write time output to", type = str)
    parser.add_argument("incompleteJobsFile", help = "File to write incomplete job names to", type = str)
    parser.add_argument("maxNumberProcessess", help = "Maximum number of processes to start", type = int)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    res = {}    ## keys are HLA alleles
    jobs = set()
    failedJobs = set()

    timeRes = []


    ## get list of jobs that were submitted.
    timecheck = time.time()
    print("Reading in list of all jobs that were submitted for processing...")
    for line in open(args.scriptReferenceFile, "r"):
        jobs.add(line.split("\t")[0])

    print("Complete. Took {:.2f} seconds...".format(time.time()-timecheck))


    ## go through jobs and tally data.
    timecheck = time.time()
    print("Checking each job...")
    try:
        print("Making job processing slaves...")
        procs = []
        pqs = []
        out_qs = []
        for i in range(args.maxNumberProcessess):
            q = mp.Queue()
            oq = mp.Queue()
            p = mp.Process(target=processJob, args=(q, oq, i+1))
            procs.append(p)
            p.start()
            pqs.append(q)
            out_qs.append(oq)

        proc_ind = 0
        job_num = 0
        print("Submitting jobs to slaves...")
        for jobname in jobs:
            if job_num % 100 == 0:
                print("{:,} jobs submitted to slaves...".format(job_num), end="\r")


            #if DEBUG: print("Processing job {}".format(jobname))
            direc = os.path.join(os.path.abspath(args.resultsDir),jobname)

            ## add to queue
            pqs[proc_ind].put([jobname,direc])
            proc_ind += 1
            if proc_ind >= args.maxNumberProcessess:
                proc_ind = 0

            job_num += 1



        for pq in pqs:
            pq.put(SENTINEL)

        print("\nTook {:.2f} seconds...".format(time.time() - timecheck))

        print("All jobs submitted, beginning to write timing information in {:,} line blocks...".format(maxBufferSize))

        timecheck = time.time()

        timeOut = open(args.timeOutputFile, "w")
        timeOut.write("jobname\thla\tpepLen\tdate\tdatetime\tduration\n")

        failedOut = open(args.incompleteJobsFile, "w")

        resultBuffer = []
        numInBuffer = 0
        numBlocks = 0
        blockTime = time.time()

        numJobs = 0
        proc_ind = 0
        while numJobs < len(jobs):
            if not out_qs[proc_ind].empty():
                resultChunk = out_qs[proc_ind].get()
                for item in resultChunk:
                    cleanStat, bindStat, timeStat = item

                    ## if job failed:
                    if not cleanStat[1]:
                        failedOut.write("{}\n".format(cleanStat[0]))

                    else:
                        ## hla binding
                        hla, hlaNums = bindStat
                        if hla not in res:
                            res[hla] = {50: [0 for x in range(0,4)], 100: [0 for x in range(0,4)], 500: [0 for x in range(0,4)]}
                        for cutoff in hlaNums:
                            for i in range(len(hlaNums[cutoff])):
                                res[hla][cutoff][i] += hlaNums[cutoff][i]


                        ## time
                        resultBuffer.append("\t".join(str(x) for x in timeStat))
                        numInBuffer += 1
                        if numInBuffer % 100 == 0:
                            print("Buffer contains {:,} lines...".format(numInBuffer), end="\r")

                    numJobs += 1

            proc_ind += 1
            if proc_ind >= args.maxNumberProcessess:
                proc_ind = 0

            if numInBuffer >= maxBufferSize:
                numBlocks += 1
                print("\nBlock {} took {:.2f} seconds...writing...".format(numBlocks, time.time() - blockTime))
                if DEBUG: print(resultBuffer)
                timeOut.write("\n".join(resultBuffer))
                timeOut.write("\n")
                resultBuffer = []
                numInBuffer = 0
                blockTime = time.time()

        print("\nWriting final time block...")
        timeOut.write("\n".join(resultBuffer))
        timeOut.write("\n")
        timeOut.close()

        failedOut.close()


        print("Took {:.2f} seconds...".format(time.time() - timecheck))

    except KeyboardInterrupt:
        print("Keyboard Interruption: ".format(sys.exc_info()[0]))
        print(traceback.format_exc())
        print("Please wait while processes close...")
        for q in pqs:
            q.empty()
            q.put(SENTINEL)
        for p in procs:
            p.join()

        sys.exit()



    timecheck = time.time()
    print("Writing HLA data...")

    out = open(args.resultsOutputFile, "w")
    out.write("hla\tco50_8mer\tco50_9mer\tco50_10mer\tco50_11mer\tco100_8mer\tco100_9mer\tco100_10mer\tco100_11mer\tco500_8mer\tco500_9mer\tco500_10mer\tco500_11mer\n")
    for hla in res:
        out.write("{}\t{}\t{}\t{}\n".format(hla, "\t".join(str(x) for x in res[hla][50]), "\t".join(str(x) for x in res[hla][100]), "\t".join(str(x) for x in res[hla][500])))
    out.close()

    print("Took {:.2f} seconds...".format(time.time() - timecheck))

    print("done.")
