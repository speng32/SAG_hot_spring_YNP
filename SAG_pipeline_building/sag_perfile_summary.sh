#!/bin/bash

# get sag per file summary after adding viral info into the blast result

echo
echo Generating per file summary of SAG blast hits
echo

PARTITIONMEMBER="./"
BLASTRESADDFOLDER="./"
OUTPUTFOLDER="./"
OUTPUTFILE=""

while [[ $# > 1 ]]; do
        key="$1"

        case $key in
                -f|--blast_res_add_folder)
                        BLASTRESADDFOLDER="$2"
                        shift # past argument
                        ;;
                -p|--partition_file)
                        PARTITIONMEMBER="$2"
                        shift
                        ;;
                -s|--output_folder)
                        OUTPUTFOLDER="$2"
                        shift
                        ;;
                -o|--output_file)
                        OUTPUTFILE=$2
                        shift
                        ;;
                *)
                        # unknown option
                        ;;
        esac
        shift # past argument or value
done

if [ ! -d $BLASTRESADDFOLDER ]; then
        echo Blast result with viral info added folder not exist!
        exit 1
elif [ ! -e $PARTITIONMEMBER ]; then
        echo Viral partition member file not exist!
        exit 1
else
        mkdir -p $OUTPUTFOLDER
        python sag_rs.py $PARTITIONMEMBER $BLASTRESADDFOLDER $OUTPUTFOLDER $OUTPUTFILE
        echo SAG per file summary generated
        echo
fi
