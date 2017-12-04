#!/bin/bash

# to extract the reads which have hits in the database

echo
echo Extracting reads with hits in virus database
echo

BLASTRESADDFOLDER="./"
FASTAFOLDER="./"
OUTPUTFOLDER="./"

while [[ $# > 1 ]]; do
        key="$1"

        case $key in
                -f|--fasta_file_folder)
                        FASTAFOLDER="$2"
                        shift # past argument
                        ;;
                -b|--blast_res_add_folder)
                        BLASTRESADDFOLDER=$2
                        shift
                        ;;
		-o|--output_folder)
                        OUTPUTFOLDER="$2"
                        shift
                        ;;
                *)
                        # unknown option
                        ;;
        esac
        shift # past argument or value
done

if [ ! -d $BLASTRESADDFOLDER ]; then
	echo Blast hits folder not exist!
	exit 1
elif [ ! -d $FASTAFOLDER ]; then
	echo Original reads fasta folder not exist!
	exit 1
else
	mkdir -p $OUTPUTFOLDER
	python extract_v.py $BLASTRESADDFOLDER $FASTAFOLDER $OUTPUTFOLDER
	echo Reads extraction finished!
	echo
fi	
