#!/bin/bash

# separate blast hits based on viral partitions rather than SAG id

echo
echo Separating Fasta files based on viral paritions rather than SAG id
echo 

FASTAEXTRACTFOLDER="./"
MAINVPNFOLDER="./"
OUTPUTFOLDER="./"

while [[ $# > 1 ]]; do
        key="$1"

        case $key in
                -e|--fasta_extracted_folder)
                        FASTAEXTRACTFOLDER="$2"
                        shift # past argument
                        ;;
                -m|--main_viral_partition_folder)
                        MAINVPNFOLDER="$2"
                        shift
                        ;;
                -o|--output_folder)
                        OUTPUTFOLDER="$2"
                        shift
			;;
	esac
	shift
done

if [ ! -d $FASTAEXTRACTFOLDER ]; then
	echo Extracted fasta file folder not exist!
	exit 1
elif [ ! -d $MAINVPNFOLDER ]; then
	echo Main viral partition folder not exist!
	exit 1
else
	mkdir -p $OUTPUTFOLDER
	Rscript organize_rbvp.R $FASTAEXTRACTFOLDER $MAINVPNFOLDER $OUTPUTFOLDER
fi

echo Separation complete
echo
