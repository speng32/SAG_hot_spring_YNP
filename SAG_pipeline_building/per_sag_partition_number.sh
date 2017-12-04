#!/bin/bash

# get the per sag file reads-viral_partition_number relationship file, unknown partition will be removed

echo
echo Getting the reads-viral_partition_number relationship file
echo

BLASTRESADDFOLDER="./"
MAINVPNFILE="./"
OUTPUTFOLDER="./"

while [[ $# > 1 ]]; do
        key="$1"

        case $key in
                -f|--blast_res_add_folder)
                        BLASTRESADDFOLDER="$2"
                        shift # past argument
                        ;;
                -m|--main_viral_partition_file)
                        MAINVPNFILE="$2"
                        shift
                        ;;
                -o|--output_folder)
                        OUTPUTFOLDER="$2"
                        shift
                        ;;
	esac
	shift
done

OUTPUTFOLDERALL=$OUTPUTFOLDER"_all"
OUTPUTFOLDERMAIN=$OUTPUTFOLDER"_main"

if [ ! -d $BLASTRESADDFOLDER ]; then
	echo Blast result with viral info added folder not exist!
        exit 1
elif [ ! -e $MAINVPNFILE ]; then
	echo Main viral partition list file not exist!
	exit 1
else
	mkdir -p $OUTPUTFOLDERALL

	for file in $(ls $BLASTRESADDFOLDER/*.tabout); do
		id=`echo $file | cut -c8-10`
		outfile=$(basename $file .tabout).virpar.all
		cat $file | cut -f1,13 | sed '/unk/d' | sort -k1 -k2 | uniq > $OUTPUTFOLDERALL/$outfile
	done

	if [ "$(ls -A $OUTPUTFOLDERALL )" ]; then
		mkdir -p $OUTPUTFOLDERMAIN
		Rscript per_sag_partition_number_main.R $MAINVPNFILE $OUTPUTFOLDERALL $OUTPUTFOLDERMAIN
	fi

	echo Finished
	echo 
fi
