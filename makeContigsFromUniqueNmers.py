'''
Make Contigs From Unique N-mers
Join together all unique n-mers with n-1 overlaps. Increases protein length to be submitted to netMHCpan, decreasing compute time.

Date: August 23, 2016
@author: sbrown
'''

## Import Libraries
import sys
import argparse

DEBUG = False
VERB = False




## For pseudocode, assuming 8mers

## for each 8mer:
    ## findMateOrAdd(8mer)

## def findMateOrAdd(8mer):
    ## if nterm7mer in ctermDict:
        ## seq = ctermDict[7mer] + 8mer[7:]
        ## remove ctermDict[7mer] and ntermDict[shifted7mer]
        ## recurse on seq
    ## else if cterm7mer in ntermDict:
        ## seq = 8mer[:1] + ntermDict[7mer]
        ## remove ntermDict[7mer] and ctermDict[shifted7mer]
        ## recurse on seq
    ## else
        ## ntermDict[nterm7mer] = 8mer[7:]
        ## ctermDict[cterm7mer] = 8mer[:1]

def findMateOrAdd(seq):
    ntermMer = seq[:n-1]
    ntermRemain = seq[n-1:]
    ctermMer = seq[len(seq)-n+1:]
    ctermRemain = seq[:len(seq)-n+1]
    
    if ntermMer in ctermDict:
        ## join the sequence
        joinSeq = ctermDict[ntermMer][0] + ntermMer + ntermRemain
        ## get the seq of the joining partner
        joinPair = ctermDict[ntermMer][0] + ntermMer
        ## remove the partner from dictionaries
        del ctermDict[ntermMer][0]
        if len(ctermDict[ntermMer]) == 0:
            del ctermDict[ntermMer]
        ntermDict[joinPair[:n-1]].remove(joinPair[n-1:])
        if len(ntermDict[joinPair[:n-1]]) == 0:
            del ntermDict[joinPair[:n-1]]
        ## recurse on joined sequence
        findMateOrAdd(joinSeq)
    elif ctermMer in ntermDict:
        ## join the sequence
        joinSeq = ctermRemain + ctermMer + ntermDict[ctermMer][0]
        ## get the seq of the joining partner
        joinPair = ctermMer + ntermDict[ctermMer][0]
        ## remove the partner from dictionaries
        del ntermDict[ctermMer][0]
        if len(ntermDict[ctermMer]) == 0:
            del ntermDict[ctermMer]
        ctermDict[joinPair[len(joinPair)-n+1:]].remove(joinPair[:len(joinPair)-n+1])
        if len(ctermDict[joinPair[len(joinPair)-n+1:]]) == 0:
            del ctermDict[joinPair[len(joinPair)-n+1:]]
        ## recurse on joined sequence
        findMateOrAdd(joinSeq)
    else:
        ## add seq to dictionaries
        ## there might already be the ntermMer in the ntermDict...
        if ntermMer not in ntermDict:
            ntermDict[ntermMer] = []
        ntermDict[ntermMer].append(ntermRemain)
        if ctermMer not in ctermDict:
            ctermDict[ctermMer] = []
        ctermDict[ctermMer].append(ctermRemain)

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Make contigs from unique n-mers")
    ## add_argument("name", "(names)", metavar="exampleOfValue - best for optional", type=int, nargs="+", choices=[allowed,values], dest="nameOfVariableInArgsToSaveAs")
    parser.add_argument("uniqueNmerFile", help = "File containing unique nmers", type = str)
    parser.add_argument("N", help = "Value of N", type = int)
    parser.add_argument("outputFile", help = "File to write output to", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB
    n = args.N

    ntermDict = {}
    ctermDict = {}
    
    print("Joining nmers...")
    
    numNmers = 0

    for line in open(args.uniqueNmerFile, "r"):
        numNmers += 1
        findMateOrAdd(line.rstrip())

    print("Writing output...")

    numContigs = 0
    listOfLengths = []

    out = open(args.outputFile, "w")
    for ntermMer in ntermDict:
        for resid in ntermDict[ntermMer]:
            out.write("{}{}\n".format(ntermMer, resid))
            numContigs += 1
            listOfLengths.append(len("{}{}".format(ntermMer, resid)))
    out.close()

    print("There were {} unique {}mers.".format(numNmers, n))    
    print("After joining, there are {} unique sequences.".format(numContigs))
    print("The average length after joining is {}, with a min of {} and a max of {}.".format(sum(listOfLengths)/len(listOfLengths), min(listOfLengths), max(listOfLengths)))
    print("done.")

