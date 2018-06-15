'''
Parse NetMHCpan output
Parse the output from NetMHCpan 3.0

Date: September 19, 2016
@author: sbrown

Edited September 20, 2016:
 - Can distill the data even more. Do not need to store HLA and peptide...
 - HLA is in the filename. Peptide can be inferred from protein file.
   - ordering of peptides is the same as they occur in the contigs.
'''

## Import Libraries
import sys
import argparse

DEBUG = False
VERB = False


if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Parse NetMHCpan 3.0")
    ## add_argument("name", "(names)", metavar="exampleOfValue - best for optional", type=int, nargs="+", choices=[allowed,values], dest="nameOfVariableInArgsToSaveAs")
    parser.add_argument("netMHCpan_file", help = "File to parse", type = str)
    parser.add_argument("output_file", help = "File to write output to", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    ## will output HLA, peptide, ic50.
    ## HLA is a bit redundant because for self-immunopeptidome work each file is all one HLA, but to make this generalizable will keep it.

    #print("Parsing {}".format(args.netMHCpan_file))

    out = open(args.output_file, "w")

    INPREDICTIONS = False
    for line in open(args.netMHCpan_file, "r"):
        if INPREDICTIONS:
            if line.strip() == "":
                ## finished this chunk
                INPREDICTIONS = False
            elif not line.startswith("---"):
                line = line.strip().rstrip().split()
                #out.write("{}\t{}\t{}\n".format(line[1], line[2], line[12]))
                out.write("{}\n".format(line[12]))

        elif line.strip().startswith("Pos"):
            INPREDICTIONS = True

    #print("done.")
