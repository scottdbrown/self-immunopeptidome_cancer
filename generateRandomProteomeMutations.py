'''
Generate Random Proteome Mutations
Take the reference proteome and generate a random mutation list along with peptide lists

Date: November 30, 2017
@author: sbrown
'''

## Import Libraries
import sys
import argparse
import os
import time
import random

DEBUG = False
VERB = False

AMINO_ACIDS_FOR_COUNTS = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
COUNTS = {}



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


def weighted_choice(items):
    '''items is a list of tuples in the form (item, weight)'''
    '''returns random item for list given the weighted probability'''
    weight_total = sum((item[1] for item in items))
    n = random.uniform(0, weight_total)
    for item, weight in items:
        if n < weight:
            return item
        n = n - weight
    return item

def processMutation(MUTI, PROTI):
    ## Get the position (0-based) on this sequence:
    prot_pos = mut_positions[MUTI] - PROTI - 1
    ref_aa = prot_seq[prot_pos:prot_pos+1]
    
    if ref_aa not in AMINO_ACIDS_FOR_COUNTS:
        mut_aa = "X"
    else:
        if args.transition_counts:
            ## pick random different amino acid based on transition frequencies
            mut_aa = weighted_choice(TRANS_FREQS[ref_aa])
        else:
            ## Pick a random different amino acid
            mut_aa = random.sample(AMINO_ACIDS^set(ref_aa), 1)[0]
            ## Note: The above behaves unexpectedly if ref_aa not in AMINO_ACIDS - resulting set is the union of the two sets. This resulted in an "X1X" mutation. Will leave for now since I want to mutate the same sites, but keep this in mind. Also generated a Z.

        ## Add aa change to COUNTS
        COUNTS[ref_aa][mut_aa] += 1
    
    ## Get minimal peptide reference sequence
    minPepStart = prot_pos - 11 + 1
    minPepEnd = prot_pos + 11 - 1
    if minPepStart < 0:
        minPepStart = 0
    if minPepEnd >= len(prot_seq):
        minPepEnd = len(prot_seq) - 1
    refMinPep = prot_seq[minPepStart:minPepEnd + 1]

    ## Get variant position (0-based)
    varPos = prot_pos - minPepStart

    ## mutate reference sequence
    mutPep = refMinPep[:varPos] + mut_aa + refMinPep[varPos+1:]

    ## Add entry to mutation list
    return {"proteome_position":mut_positions[MUTI], "protein":cur_prot[1:], "mutation": "{}{}{}".format(ref_aa, prot_pos+1, mut_aa), "minimal_peptide": mutPep, "variant_position": varPos+1}


if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "Generate Random Proteome Mutations")
    ## add_argument("name", "(names)", metavar="exampleOfValue - best for optional", type=int, nargs="+", choices=[allowed,values], dest="nameOfVariableInArgsToSaveAs")
    parser.add_argument("reference_fasta", help = "Fasta file with proteome reference", type = str)
    parser.add_argument("output_dir", help = "Directory to write output to", type = str)
    parser.add_argument("count_output_file", help = "File to write amino acid change counts to", type = str)
    parser.add_argument("--seed", help = "Random seed", type = int, default=None)
    parser.add_argument("--num_mutations", help = "Number of mutations to generate", type = int, default=None)
    parser.add_argument("--transition_counts", help = "File with amino acid transition counts.")
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

    log_print("STATUS","Building array...")
    for aa1 in AMINO_ACIDS_FOR_COUNTS:
        COUNTS[aa1] = {}
        for aa2 in AMINO_ACIDS_FOR_COUNTS:
            COUNTS[aa1][aa2] = 0


    ## process transition frequencies, if present
    log_print("STATUS", "Processing transition counts...")
    if args.transition_counts:
        ## first pass through to get counts for each reference amino acid
        ref_aa_counts = {}
        HEADER = True
        for line in open(args.transition_counts, "r"):
            if HEADER:
                HEADER = False
            else:
                line = line.rstrip().split("\t")
                ref = line[0]
                mut = line[1]
                count = int(line[2])

                if ref not in ref_aa_counts:
                    ref_aa_counts[ref] = 0
                ref_aa_counts[ref] += count

        ## then another pass through to create arrays
        TRANS_FREQS = {}
        HEADER = True
        for line in open(args.transition_counts, "r"):
            if HEADER:
                HEADER = False
            else:
                line = line.rstrip().split("\t")
                ref = line[0]
                mut = line[1]
                count = int(line[2])

                if ref not in TRANS_FREQS:
                    TRANS_FREQS[ref] = []

                TRANS_FREQS[ref].append([mut,count/ref_aa_counts[ref]])


    random.seed(args.seed)

    ## determine proteome length
    log_print("STATUS","Determining total proteome length...")
    TOTAL_PROTEOME_LENGTH = 0
    for line in open(args.reference_fasta, "r"):
        if not line.startswith(">"):
            TOTAL_PROTEOME_LENGTH += len(line.rstrip())


    ## make sure that number of mutations requested < total proteome length
    if args.num_mutations > TOTAL_PROTEOME_LENGTH:
        log_print("ERROR", "Requested {} mutations, but total proteome length is only {}.".format(args.num_mutations, TOTAL_PROTEOME_LENGTH))
        sys.exit()


    ## generate random positions to mutate.
    log_print("STATUS", "Generating set of positions to mutate...")
    mut_positions = random.sample(set(range(1,TOTAL_PROTEOME_LENGTH+1)), args.num_mutations)
    log_print("STATUS", "Sorting list of mutation positions...")
    mut_positions.sort()

    ## add position to the end of mut_positions as a "dummy flag"
    ## Make it a position that is not in this proteome.
    ## hacky solution to checking if the MUT_IND goes past the end.
    mut_positions.append(TOTAL_PROTEOME_LENGTH+1)

    ## walk through proteome and mutate as you go.
    log_print("STATUS", "Walking through proteome and producing mutations...")
    MUTATIONS = [{} for x in range(args.num_mutations)]
    AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
    PROT_IND = 0
    MUT_IND = 0
    PEPTIDES = {8:[], 9:[], 10:[], 11:[]}

    cur_prot = ""
    prot_seq = ""
    for line in open(args.reference_fasta, "r"):
        line = line.rstrip()
        if line.startswith(">"):
            ## see if next mutation(s) position is in the cur_prot
            while mut_positions[MUT_IND] <= PROT_IND + len(prot_seq):
                ## it is on this line
                MUTATIONS[MUT_IND] = processMutation(MUT_IND, PROT_IND)

                ## Generate all 8-11mers
                varPos = MUTATIONS[MUT_IND]["variant_position"] - 1 ## conver to 0-based
                for n in [8,9,10,11]:
                    for i in range(0,len(MUTATIONS[MUT_IND]["minimal_peptide"])+1-n):
                        ## if mutation is within the peptide:
                        if varPos >= i and varPos <= i+n-1:
                            PEPTIDES[n].append(MUTATIONS[MUT_IND]["minimal_peptide"][i:i+n])

                ## Move on to next mutation
                MUT_IND += 1

            cur_prot = line
            PROT_IND += len(prot_seq)
            prot_seq = ""
        else:
            prot_seq += line
            

    ## handle the last sequence (that may contain a mutation)
    while mut_positions[MUT_IND] <= PROT_IND + len(prot_seq):
        ## it is on this line
        MUTATIONS[MUT_IND] = processMutation()

        ## Generate all 8-11mers
        varPos = MUTATIONS[MUT_IND]["variant_position"] - 1 ## conver to 0-based
        for n in [8,9,10,11]:
            for i in range(0,len(MUTATIONS[MUT_IND]["minimal_peptide"])+1-n):
                ## if mutation is within the peptide:
                if varPos >= i and varPos <= i+n-1:
                    PEPTIDES[n].append(MUTATIONS[MUT_IND]["minimal_peptide"][i:i+n])

        ## Move on to next mutation
        MUT_IND += 1

    ## print out mutation list
    log_print("STATUS", "Writing mutation metadata...")
    out = open(os.path.join(args.output_dir, "mutations.tsv"), "w")
    out.write("index\tproteome_position\tprotein\tmutation\tminimal_peptide\tvariant_position\n")
    mut_ind = 0
    for mut in MUTATIONS:
        mut_ind += 1
        out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(mut_ind, mut["proteome_position"], mut["protein"], mut["mutation"], mut["minimal_peptide"], mut["variant_position"]))
    out.close()


    ## print out the 8-11mer peptides
    log_print("STATUS", "Writing 8-11mer peptides...")
    outs = {8: open(os.path.join(args.output_dir,"8mer_peptides.txt"), "w"), 9: open(os.path.join(args.output_dir,"9mer_peptides.txt"), "w"), 10: open(os.path.join(args.output_dir,"10mer_peptides.txt"), "w"), 11: open(os.path.join(args.output_dir,"11mer_peptides.txt"), "w")}

    for n in PEPTIDES:
        outs[n].write("\n".join(PEPTIDES[n]))
        outs[n].close()

    ## print out aa change counts
    log_print("STATUS", "Writing output...")
    out = open(args.count_output_file, "w")
    out.write("reference\talternate\tcount\n")

    for aa1 in COUNTS:
        for aa2 in COUNTS[aa1]:
            out.write("{}\t{}\t{}\n".format(aa1, aa2, COUNTS[aa1][aa2]))

    out.close()

    ## done
    log_print("STATUS", "Complete.")
