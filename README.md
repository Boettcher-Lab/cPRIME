# 🧬 cPRIME
We introduce cPRIME (prime editing with cDNA readout), a scalable method for assessing variant effects on RNA metabolism. cPRIME employs prime editing to insert genetic variants into their endogenous genomic locus, and quantifies changes in RNA abundance by comparing editing frequencies between genomic DNA and complementary DNA.

This repository contains code to reproduce results of our publication "cPRIME enables the systematic mapping of genetic variant effects on RNA metabolism".

## ⭐ Highlights
- scripts to compute synonymous marker variants
- synonymous marker filtering
- code to select pegRNAs from PrimeDesign output
- fast and easy to use pipeline for read counting that contain variants introduced by prime editing

## ℹ️ Overview
This repository is structured in two sections:
- ex27_screen
  
  -> contains code and associated data to reproduce results of the intial exon 27 screen of the NF1 gene
- multi_exon_screen
  
  -> contains code and associated data to reproduce results of the multi exon screen of the NF1 gene

## ex27_screen


### Count_edits_in_fastq
Contains Bash and Python scripts to count edits in FASTQ files.

Count_edits.sh -> you have to change in line 119 the reference file depending on if you want to search for hits in:
 - gDNA (with synonymous marker)
 - cDNA (with synonymous marker)
 - gDNA (with wild type sequence (WT); no intended or synonymous marker edit)
 - gDNA (without synonymous marker; only intended edit)

Line 119:  awk '{ print $2 }' Searching_sequences_gDNA.txt > sequences.txt

### search_sequences
Contains the search sequences that were used to count the number of reads containing combined variants.
Use these search sequences as described in Count_edits_in_fastq/Count_edits.sh


### Design_pegRNAs_additional_STOP_codon
Contains R script to design PrimeDesign input files for nonsense mutations and synonymous mutations on codon positions in exon 27 that are not described in cBioPortal. Also includes R scripts to generate the reference files with the search sequences, to determine the number of each intended edit in the fastq files.

- Design_input_sequences_for_PrimeDesign_STOP_and_synonymous_edit.R -> design input sequences for PrimeDesign
- Design_searching_sequences_for_STOP_and_synonymous_edit_cDNA.R -> design searching sequences for additional nonsense (STOP) and synonymous mutations for cDNA
- Design_searching_sequences_for_STOP_and_synonymous_edit_gDNA.R -> design searching sequences for additional nonsense (STOP) and synonymous mutations for gDNA
- Design_WT_sequences_for_STOP_and_synonymous_edit_cDNA.R -> design WT sequences for additional nonsense (STOP) and synonymous mutations for cDNA
- Design_WT_sequences_for_STOP_and_synonymous_edit_gDNA.R -> design WT sequences for additional nonsense (STOP) and synonymous mutations for gDNA

- Chr17_full_sequence.txt -> Sequence of chromosome 17 from hg19 genome
- NF1-Exons-position.txt -> Information about start and end position of each exon from NF1 in chromosome 17

Additional nonsense mutations or synonymous edits that are missing in the final Overview_names.xlsx list are removed manually for example all InDels are removed.


### Design_pegRNAs_cBioPortal_mutations
Contains R script to design PrimeDesign input files for the mutations listed in cBioPortal for exon 27. Also includes R scripts to generate the reference files with the search sequences, to determine the number of each intended edit in the fastq files.

- cBioPortal_27Jan2022.tsv -> cBioPortal mutations for NF1 (Download: 27 January 2022)
- Chr17_full_sequence.txt -> Sequence of chromosome 17 from hg19 genome
- NF1-Exons-position.txt -> Information about start and end position of each exon from NF1 in chromosome 17
- Exon_splice_information.txt -> Information about number of bases of the first and last codon of each exon. 
For example: Start_codon_bp_count 3 means that the first codon of the exon starts with three bases; Start_codon_bp_count 2 means that the exon starts with two bases while the "third" base for this codon is located on the exon before (her the last base).
- Design_input_sequences_mutations_for_PrimeDesign.R -> design input sequences for PrimeDesign for mutations
- Design_input_sequences_synonymous_edit_for_PrimeDesign.R -> design input sequences for PrimeDesign for synonymous edit
- Design_searching_sequences_for mutations_cDNA.R -> design searching sequences for mutations for cDNA
- Design_searching_sequences_for mutations_gDNA.R -> design searching sequences for mutations for gDNA
- Design_searching_sequences_for_synonymous_edit_cDNA.R -> design searching sequences for synonymous edit for cDNA
- Design_searching_sequences_for_synonymous_edit_gDNA.R -> design searching sequences for synonymous edit for gDNA
- Design_WT_searching_sequences_for_mutations_cDNA.R -> design WT sequences for mutations for cDNA
- Design_WT_searching_sequences_for_mutations_gDNA.R -> design WT sequences for mutations for gDNA
- Design_WT_searching_sequences_for_synonymous_edit_cDNA.R -> design WT sequences for synonymous edit for cDNA
- Design_WT_searching_sequences_for_synonymous_edit_gDNA.R -> design WT sequences for synonymous edit for gDNA

Mutations or synonymous edits that are missing in the final Overview_names.xlsx list are removed manually for example all InDels are removed.


### Odds_Ratio
Contains the R script and reference files to calculate the Odds Ratio.



## multi_exon_screen

### library_design
- create_gff_db.py -> creates a local database from a gff file which is used by the other scripts to get transcript information from the NF1 gene. Can be used from the command line like python3 create_gff_db.py input.gff output.db
- filter_vcf.py -> filters a vcf given a genomic range. Can be used from the command line like python3 filter_vcf.py -i /path/to/input.vcf -o /path/to/output.vcf -c chromosome -s start_position -e stop_position
- functions.py -> contains functions which are used in multiple scripts
- library_design.ipynb -> contains code that was used to compute and filter the variants and pegRNAs in the multi exon screen
- prioritize_exons.ipynb -> contains code that was used to get the most promising exons
- write_searching_sequences.ipynb -> contains code to compute searching sequences for the combined variants

### search_sequences
Contains the output of library_design/write_earching_sequences.ipynb. These sequences were input to count_edits_in_fastq/count_edits_fast.sh.

### count_edits_in_fastq
Contains the advanced pipeline to count edits. Simply run
count_edits_fast.sh -p /path/to/fastq_files_folder -s /path/to/searching_sequences.txt

### analysis
Contains code to reproduce figures 3 and 4 of the manuscript and other analyses related to the multi exon screen.

## ✍️ Authors
If you are interested in our work feel free to visit our website (https://www.boettcher-lab.net/) to learn more.
