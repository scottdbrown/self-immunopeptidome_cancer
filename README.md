# self-immunopeptidome_cancer
Scripts relevant to Self-immunopeptidome and Cancer manuscript


## Environment
* Python v3.4.5 (and python_requirements.txt)
* Samtools v0.1.8
* sqlite3 v3.6.20
* NetMHCpan v3.0

## Tasks described here
* [Condensing proteomes](#procedure-to-condense-a-proteome)
* [Prepare predictions](#prepare-predictions)
* [Generating random amino acid changes](#generating-random-proteome-mutations)
* [Getting variant RNA-seq support](#getting-variant-rna-seq-read-support)

## Procedure to condense a proteome
This will condense a given proteome (.fasta(s)) into sets of artificial proteins ("contigs")


#### 8mers:
Extract all unique n-mers from the proteome file(s)
```bash
$ python processUniqueNmersProteome.py output_directory/8mers.txt -n 8 --fasta proteome_file1.fasta [proteome_file2.fasta ...]
```

Combine unique n-mers into a set of artificial proteins ("contigs") containing only each unique n-mer exactly once (when parsing with a sliding window of size n)
```bash
$ python makeContigsFromUniqueNmers.py output_directory/8mers.txt 8 output_directory/8mers_contigs.txt
```

Record lengths of the resulting contigs
```bash
$ awk '{{print length($0);}}' output_directory/8mers_contigs.txt > output_directory/8mers_lengths.txt
```

#### Repeat for 9, 10, and 11mers


## Prepare predictions

First, get list of all HLA available for prediction in NetMHCpan
```bash
$ /path/to/netMHCpan-3.0/netMHCpan -listMHC | awk '{print $1}' | grep -e "HLA-[ABC]" > allHLAI.txt
```

Note, this was designed to create input files for [clusterTAS](https://github.com/scottdbrown/bcgsc-scripts/blob/master/clusterTAS) cluster submission and management script locally at the BC Genome Sciences Centre. This will likely need modification to run on your system. Format of each line of scripts.sh is `job_name bash_command;`.

Prepare NetMHCpan invocations using the condensed proteomes. First edit lines 22-24 to paths on your system. Then run as:
```bash
$ python prepareJobs.py --species HUMAN --contig8mer output_directory/8mers_contigs.txt --contig9mer .output_directory/8mers_contigs.txt --contig9mer output_directory/10mers_contigs.txt --contig11mer output_directory/11mers_contigs.txt --contigsPerJob 1000 --hlaAlleleList allHLAI.txt --destDir /path/to/output/jobs_dir/
```
This breaks the input into many smaller individual jobs, depending on the proteome size.
Example of one line of the scripts.sh file is:
```bash
HUMAN_HLA-B13-23_8_32	source /home/sbrown/bin/pythonvenv/python3/bin/activate; /path/to/netMHCpan-3.0/netMHCpan -tdir tmpdirXXXXXX -a HLA-B13:23 -l 8 -f prot8_32_HUMAN.fa > HUMAN_HLA-B13-23_8_32.pMHC; python /home/sbrown/scripts/parseNetMHCpanOutput.py HUMAN_HLA-B13-23_8_32.pMHC HUMAN_HLA-B13-23_8_32.pMHC.parsed; rm HUMAN_HLA-B13-23_8_32.pMHC;
```
Note, since the HLA is in the file name, and the prot8_32_HUMAN.fa is the lest of peptides, this is parsed down to just be the IC50 scores for each peptide (in the same order as prot8_32_HUMAN.fa) to save space.

Within the folder holding the results of all the predictions, we will check to see that all jobs completed successfully, and get simple summaries
```bash
$ python tallyParsedData_multiProc.py -v ../scripts.sh analysis_results/ singleHLAdata.tsv timeCharacteristics.tsv failedJobs.txt 16
```

Make sure that failedJobs.txt is an empty file before continuing.

Now we will build an sqlite3 (v3.6.20) database to hold information on pMHC binding.
```bash
$ python makeDatabaseOfBinders.py HUMAN /path/to/results/ allHLAI.txt HUMAN_binders.db 16
```
Note: Database holds all peptides and hla, but only pMHC interactions (binders) with IC50 < 500 nM.

Details on the schema of the created database:
```bash
$ sqlite3 HUMAN_binders.db
```
```sql
SQLite version 3.6.20
Enter ".help" for instructions
Enter SQL statements terminated with a ";"
sqlite> .tables
binders  hla      peptide
sqlite> .schema peptide
CREATE TABLE peptide(id INT, sequence TEXT);
CREATE INDEX peptide_ind ON peptide(id);
CREATE INDEX peptide_seq ON peptide(sequence);
sqlite> .schema hla
CREATE TABLE hla(id INT, allele TEXT);
sqlite> .schema binders
CREATE TABLE binders(hla_id INT, pep_id INT, ic50 REAL);
CREATE INDEX binder_hla_ind ON binders(hla_id);
CREATE INDEX binder_pep_ind ON binders(pep_id);
```
Note that indices need to be created manually after running `makeDatabaseOfBinders.py`.


## Creating SQLite3 database from flat files

If you are using our human immunopeptidome data (downloaded from: __________), you can generate a comparable SQLite3 database by running:
```bash
$ python makeDatabaseFromFlatFiles.py human_immunopeptidome_database_flat/ HUMAN_binders.db
```

Note that this requires loading all data into memory, and is quite RAM-intensive, requiring > 50GB of RAM.


Given a list of HLA genotypes, we can calculate the size of the self-immunopeptidome
```bash
$ python lookupHLAgenotypesSQL.py samples_hlaGeno.tsv HUMAN_binders.db output_selfimmunopeptidome_sizes.tsv 12
```
Where samples_hlaGeno.tsv is like:
```
sample_id HLA-A02-01_HLA-A28-01_HLA-B57-01_HLA-B27-05_HLA-C03-03_HLA-C03-03
...
```

## Generating Random Proteome Mutations

To generate random proteome mutations (note that TCGA_aaChange_counts.tsv is included in this repository, but can be substituted with your own frequencies)
```bash
$ python generateRandomProteomeMutations.py proteome_reference.fasta /path/for/output/ generated_aaChange_counts.tsv --seed 171201 --num_mutations 50000 --transition_counts TCGA_aaChange_counts.tsv
```


## Getting variant RNA-seq read support

To get transcriptome read support for SNVs:
Edit line 29 to be the path to your Samtools v0.1.8 binary.
```bash
$ python lookupMutationReadSupport.py variants_and_bams_toQuery.tsv variants_readCounts.tsv 36
```
Where variants_and_bams_toQuery.tsv is like: 
```
subject     sample    NCBI_Build      Chromosome      Start_Position  End_Position    STRAND  HGVSc   file
TCGA-VQ-A91D    TCGA-VQ-A91D-01A-11D-A410-08    GRCh37  19      15234051        15234051        -1      c.341C>T        /path/to/TCGA-VQ-A91D-01A-11R-A414-31_diseased_HG19_GSC_merge_RNA-Seq.bam
...
```