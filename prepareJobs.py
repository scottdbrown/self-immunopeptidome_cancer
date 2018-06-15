'''
Prepare Jobs
Makes files required for clusterTAS to run.

Date: September 13, 2016
@author: sbrown

Edited November 4, 2016:
    - Include result parsing (parseNetMHCpanOutput.py)
'''

## Import Libraries
import sys
import argparse
import os
import math

DEBUG = False
VERB = False

## CHANGE THESE PATHS FOR YOUR SYSTEM
PYTHON3ENV = "/home/sbrown/bin/pythonvenv/python3/bin/activate"
NETMHCPAN = "/home/sbrown/bin/netMHCpan-3.0/netMHCpan"
RESPARSER = "/home/sbrown/scripts/parseNetMHCpanOutput.py"



def processContigsWriteFiles(contigFile, n):
    #scripts = ""
    #fof = ""

    contigs = open(contigFile, "r").readlines()
    ## need to sort and then distribute.
    contigs.sort(key = len, reverse = True)
    numJobs = math.ceil(len(contigs) / args.contigsPerJob)
    if VERB: print("{} jobs will be created for {}.".format(numJobs, contigFile))
    files = ["" for x in range(0, numJobs)]

    ## add sequences to files back and forth (so that total length in each is similar)
    i = 0
    ASCENDING = True
    for seq in contigs:
        seq = seq.rstrip()
        files[i] += ">sim\n{}\n".format(seq)

        if ASCENDING:
            i += 1
        else:
            i -= 1
        ## catch edges
        if i == numJobs:
            ASCENDING = False
            i -= 1
        elif i == -1:
            ASCENDING = True
            i += 1

    ## write files and add line to script holder.
    ## need to write protein file for this job, and then add it to the fof and scripts files.
    fnum = 0

    for seqs in files:
        fnum += 1
        out = open(os.path.join(args.destDir, "prot{}_{}_{}.fa".format(n, fnum, args.species)), "w")
        out.write(seqs)
        out.close()

        if VERB: print("Going through each HLA for file {}...".format(fnum))
        hscript = ""
        hfof = ""
        for h in hlas:
            if DEBUG: print("Setting up jobs for {}...".format(hlas[h][0]))

            ## add to script and fof files.
            ## script line like "TCGA-A6-6781-01A-22D-A270-10    source /home/sbrown/bin/pythonvenv/python3/bin/activate;/home/sbrown/bin/netMHCpan-3.0/netMHCpan -tdir tmpdirXXXXXX -a HLA-C07:01 -f peptides.fa > bindingRes.pMHC;"
            hscript += "{}_{}_{}_{}\tsource {}; {} -tdir tmpdirXXXXXX -a {} -l {} -f {} > {}_{}_{}_{}.pMHC;".format(args.species, h, n, fnum, PYTHON3ENV, NETMHCPAN, hlas[h][0], n, "prot{}_{}_{}.fa".format(n, fnum, args.species), args.species, h, n, fnum)
            hscript += "python {} {}_{}_{}_{}.pMHC {}_{}_{}_{}.pMHC.parsed; rm {}_{}_{}_{}.pMHC;".format(RESPARSER, args.species, h, n, fnum, args.species, h, n, fnum, args.species, h, n, fnum)
            hscript += "\n"

            hfof += "{}\n".format(os.path.join(args.destDir, "prot{}_{}_{}.fa".format(n, fnum, args.species)))

            hlas[h][1] += hscript
            hlas[h][2] += hfof

            hscript = ""
            hfof = ""

    '''
    scriptFile = open(os.path.join(args.destDir,  "scripts.sh"), "a")
    scriptFile.write(scripts)
    scriptFile.close()

    fofFile = open(os.path.join(args.destDir, "files.fof"), "a")
    fofFile.write(fof)
    fofFile.close()
    '''


if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Prepare jobs for clusterTAS")
    ## add_argument("name", "(names)", metavar="exampleOfValue - best for optional", type=int, nargs="+", choices=[allowed,values], dest="nameOfVariableInArgsToSaveAs")
    parser.add_argument("--species", metavar = "species_name", help = "Name of species", type = str, default = None)
    parser.add_argument("--contig8mer", metavar = "file", help = "8mer contig file", type = str, default = None)
    parser.add_argument("--contig9mer", metavar = "file", help = "9mer contig file", type = str, default = None)
    parser.add_argument("--contig10mer", metavar = "file", help = "10mer contig file", type = str, default = None)
    parser.add_argument("--contig11mer", metavar = "file", help = "11mer contig file", type = str, default = None)
    parser.add_argument("--contigsPerJob", metavar = "N", help = "Number of contigs per job", type = int, default = None)
    parser.add_argument("--hlaAlleleList", metavar = "file", help = "File of HLA alleles to use", type = str, default = None)
    parser.add_argument("--destDir", metavar = "directory", help = "Directory to write output files to", type = str, default = None)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB


    hlas = {}
    for line in open(args.hlaAlleleList, "r"):
        hlas[line.rstrip().replace(":","-")] = [line.rstrip(),"",""]


    ## flush output files
    scriptFile = open(os.path.join(args.destDir,  "scripts.sh"), "w").close()
    fofFile = open(os.path.join(args.destDir, "files.fof"), "w").close()


    ## Process 8mer contigs
    if args.contig8mer:
        processContigsWriteFiles(args.contig8mer, 8)

    ## Process 9mer contigs
    if args.contig9mer:
        processContigsWriteFiles(args.contig9mer, 9)

    ## Process 10mer contigs
    if args.contig10mer:
        processContigsWriteFiles(args.contig10mer, 10)

    ## Process 11mer contigs
    if args.contig11mer:
        processContigsWriteFiles(args.contig11mer, 11)


    ## write by HLA
    scriptFile = open(os.path.join(args.destDir,  "scripts.sh"), "a")
    fofFile = open(os.path.join(args.destDir, "files.fof"), "a")

    for h in hlas:
        scriptFile.write(hlas[h][1])
        fofFile.write(hlas[h][2])

    scriptFile.close()
    fofFile.close()

    print("done.")
