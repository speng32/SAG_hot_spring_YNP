#!/bin/bash

# add viral parition number at the end of the blast hits

echo
echo Adding viral partition information to the Blast result
echo

BLASTRESFOLDER="./"
PARTITIONMEMBER="./init_files/viral_partition_membership.csv"
OUTPUTFOLDER="./"
MINIMUMLENGTH=100

while [[ $# > 1 ]]; do
        key="$1"

        case $key in
                -f|--blast_res_folder)
                        BLASTRESFOLDER="$2"
                        shift # past argument
                        ;;
		-p|--partition_file)
			PARTITIONMEMBER="$2"
			shift
			;;
		-o|--output_folder)
			OUTPUTFOLDER="$2"
			shift
			;;
		-l|--min_len)
			MINIMUMLENGTH=$2
			shift
			;;
                *)
                        # unknown option
                        ;;
        esac
        shift # past argument or value
done

#echo $BLASTRESFOLDER $PARTITIONMEMBER $OUTPTUFOLDER $MINIMUMLENGTH

if [ ! -d $BLASTERESFOLDER ]; then
	echo Blast result folder not exist!
	exit 1
elif [ ! -e $PARTITIONMEMBER ]; then
	echo Viral parition member file not exist!
	exit 1
else
	python add_p.py $BLASTRESFOLDER $PARTITIONMEMBER $OUTPUTFOLDER $MINIMUMLENGTH
	echo Add viral partition number finished
	echo
fi
