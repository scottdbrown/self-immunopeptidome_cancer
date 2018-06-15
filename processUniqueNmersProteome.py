'''
Process unique n-mers from the proteome
Goes through proteome fasta files and saves all unique n-mers

Date: August 18, 2016
@author: sbrown
'''

## Import Libraries
import sys
import argparse

DEBUG = False
VERB = False

nmers = set()

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Process unique n-mers from the proteome")
    ## add_argument("name", "(names)", metavar="exampleOfValue - best for optional", type=int, nargs="+", choices=[allowed,values], dest="nameOfVariableInArgsToSaveAs")
    parser.add_argument("--fasta", nargs="+", help = "List of fasta files", type = str)
    parser.add_argument("-n", "--nMerLength", help = "Length of n-mer to process", dest = "nmer", type = int)
    parser.add_argument("outputFile", help = "File to write results to", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    #print(args)
    ## arguments accessible as args.cmdline_arg or args.cmdflg or args.destName
    ## can test using parser.parse_args("-cmdflg value other_value".split())

    ## read fasta and cut into nmers.

    protSeq = ""
    totalCount = 0
    uniqueCount = 0

    print("Reading through fasta files...")
    for f in args.fasta:
        for line in open(f, "r"):
            if line.startswith(">"):
                ## process previous sequence
                for i in range(0,len(protSeq) - args.nmer + 1):
                    nmers.add(protSeq[i : i + args.nmer])
                    totalCount += 1
                protSeq = ""
            else:
                protSeq += line.rstrip()
    
    ## process the final protein in the file.
    for i in range(0,len(protSeq) - args.nmer + 1):
        nmers.add(protSeq[i : i + args.nmer])
        totalCount += 1
    protSeq = ""
    


    ## write file.
    print("Writing unique {}mers...".format(args.nmer))
    out = open(args.outputFile, "w")
    while len(nmers) > 0:
        out.write("{}\n".format(nmers.pop()))
        uniqueCount += 1
    out.close()

    print("There are {} total {}mers found.".format(totalCount, args.nmer))
    print("{} are unique.".format(uniqueCount))
    print("done.")
