#!/bin/bash

# check if both ends have hits or only one of the end have hits and get a summary file

echo
echo Checking if both ends have hits or single end have hits
echo

FASTAEXTRACTFOLDER="./"
OUTPUTSUMFOLDER="./"
OUTPUTFOLDER="./"

while [[ $# > 1 ]]; do
        key="$1"

        case $key in
                -e|--fasta_extracted_folder)
                        FASTAEXTRACTFOLDER="$2"
                        shift # past argument
                        ;;
                -s|--output_summary_folder)
                        OUTPUTSUMFOLDER="$2"
                        shift
                        ;;
		-o|--output_separate_fasta_folder)
			OUTPUTFOLDER="$2"
			shift
			;;
        esac
        shift
done

if [ ! -d $FASTAEXTRACTFOLDER ]; then
	echo Extracted fasta file folder not exist!
	exit 1
else
	mkdir -p $OUTPUTSUMFOLDER

        allsagsum="all_sag_pair_single_hits_summary"

	if [ -e $OUTPUTSUMFOLDER/$allsagsum ]; then
        	rm $OUTPUTSUMFOLDER/$allsagsum
	fi
	
	echo -e SAG_ID'\t'both_pair_have_hit'\t'only_single_end_have_hit'\t'total_reads_have_hits > $OUTPUTSUMFOLDER/$allsagsum

	for file in $(ls $FASTAEXTRACTFOLDER/*.fasta); do
	        persagsum=$(basename $file .fasta).read_hit_status_summary
	        awk 'NR%2==1' $file | sed -e 's/^>//;s/\/1$//g;s/\/2$//' | sort | uniq -c | sort -k1r  > $OUTPUTSUMFOLDER/$persagsum
		single=`cut -d' ' -f7 $OUTPUTSUMFOLDER/$persagsum | grep -c "1"`
       		paired=`cut -d' ' -f7 $OUTPUTSUMFOLDER/$persagsum | grep -c "2"`
		total=`grep -c "^>" $file`

		echo -e $(basename $file .fasta)'\t'$paired'\t'$single'\t'$total >> $OUTPUTSUMFOLDER/$allsagsum
	done
	echo Pair-end single-end summary generated
	echo

	mkdir -p $OUTPUTFOLDER
	Rscript separate_paired_single_hits_reads.R $FASTAEXTRACTFOLDER $OUTPUTFOLDER
	echo Fasta files separated
fi
echo
