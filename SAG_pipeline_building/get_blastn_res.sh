#!/bin/bash

# run blastn on the fasta file

echo 

# Initialize own variables:
FASTAFOLDER="./"
IDENTITY="0.95"
EVALUE="1e-20"
show_help="test"
OUTPUTFOLDER="./"
DATABASE="test"
NUMTHREADS=8

while [[ $# > 1 ]]; do
	key="$1"

	case $key in
		-f|--input_folder)
			FASTAFOLDER="$2"
			shift # past argument
			;;
		-i|--identity)
			IDENTITY="$2"
			shift # past argument
			;;
		-e|--evalue)
			EVALUE="$2"
			shift # past argument
			;;
		-o|--output_folder)
			OUTPUTFOLDER="$2"
			shift
			;;
		-d|--database)
			DATABASE="$2"
			shift
			;;
		-n|--num_threads)
			NUMTHREADS="$2"
			shift
			;;	
		-h|--help)
			echo show_help
			exit 1
			;;
		*)
			# unknown option
			;;
	esac
	shift # past argument or value
done

if [[ -n $1 ]]; then
	echo "Last line of file specified as non-opt/last argument:"
	tail -1 $1
fi

#shift $((OPTIND-1))
#echo $FASTAFOLDER $IDENTITY $EVALUE $show_help $OUTPUTFOLDER $DATABASE $NUMTHREADS
#[ "$1" = "--" ] && shift

logfile="blastn".$DATABASE.$IDENTITY.$EVALUE."log"

if [ ! -d $FASTAFOLDER ]; then
	echo Input file path not exist!
	exit 1
else
	mkdir -p $OUTPUTFOLDER
	
	if [ -e $logfile ];then
		rm $logfile
	fi	

	echo Performing blastn
	for file in $(ls $FASTAFOLDER/*.fasta); do
		fout=$(basename $file .fasta).blastn.$DATABASE.$IDENTITY.$EVALUE.tabout
		echo Entering file $file >> $logfile
		blastn -db "$DATABASE" -query $file -out $OUTPUTFOLDER/$fout -outfmt 6 -num_threads $NUMTHREADS -evalue $EVALUE -perc_identity $IDENTITY >> $logfile 2>&1
		#echo Finished
		echo >> $logfile
	done
fi

echo 
echo Blastn finished
echo Output are under $OUTPUTFOLDER, check log at $logfile
echo
