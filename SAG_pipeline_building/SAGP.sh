#!/bin/bash

## pipleline for SAG file processes

FASTQFOLDER="./"
FASTAFOLDER="init_files"
CONTIGFOLDER="./"
BLASTDBNAME="NL10_vDNA_and_Archaeal_blastdb"
BLASTRESFOLDER="./"
BLASTRESADDFOLDER="./"
BLASTMINLENGTH=100
PARTITIONMEMBER="init_files/viral_partition_membership.csv"
MAINPARITITON="init_files/main_partitions.txt"
SUMMARYOUTPUTFOLDER="./"

BLASTEVAL=1e-20
BLASTIDEN=0.95
BLASTTHD=8


# echo Choose task to perform \
# 	1.\
# 	2.\
# 	3.\
# 	4.

# fastq raw reads to fasta

# fasta reads assembly

# fasta reads blastn

echo 1
./get_blastn_res.sh -f $FASTAFOLDER  -i $BLASTIDEN -e $BLASTEVAL -d $BLASTDBNAME -n $BLASTTHD -o out1

# blastn result viral partition identifier and summary

echo 2
./add_partition.sh -f out1 -p $PARTITIONMEMBER -l $BLASTMINLENGTH -o out2

echo 3
./sag_summary.sh -p $PARTITIONMEMBER -f out2 -s out4 -o test_sag_summary.csv

echo 4
./sag_donimated_partition.sh -f out2 -p $PARTITIONMEMBER -s out4 -o test_sag_dominated_vp_summary.csv

echo 5
./sag_perfile_summary.sh -f out2 -p $PARTITIONMEMBER -s out4 -o test_sag_perfile_summary.csv

echo 6
./per_sag_partition_number.sh -f out2 -m $MAINPARITITON -o out5

# viral hit reads extract

echo 7
./extract_vreads.sh -f $FASTAFOLDER -b out2 -o out3

echo 8
./check_both_ends_have_hits_and_separate.sh -e out3 -s out6 -o out7

# dominanted partition reads isolation

echo 9
./organize_reads_by_viral_partition.sh -e out3 -m out5_main -o out8






echo Whole pipeline finished
