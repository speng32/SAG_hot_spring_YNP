#!/bin/bash

check_sag=$1
check_list=$2
ref_spe_file=$3
tmpfolder=$4
outfolder=$5

sagid=`basename $check_sag .fasta | cut -c8-10`
echo $sagid enter

if [ ! -d $tmpfolder ]; then
	mkdir $tmpfolder
fi

#echo $sagid

CNT=1
#for file in `less $chekclist`; do
while read reff; do
	refname=`basename $reff .fasta | sed "s/ref_//"`
	#echo $refname
	ANI_mod1.pl -bl ~/tools/legacy_blast/bin/blastall -fd ~/tools/legacy_blast/bin/formatdb -qr $check_sag -sb $reff -od $tmpfolder/$CNT
	sed -i "s/$/\t$refname/" $tmpfolder/$CNT/raw.blast
	CNT=$((CNT+1))
	#echo $CNT
done < "$outfolder/$check_list"

#echo $tmpfolder/$sagid'_check_dc_blast_res.tsv'
find $tmpfolder -name raw.blast -exec cat {} + > $tmpfolder/$sagid'_check_dc_blast_res.tsv'

Rscript ~/data/sags/code/check_double_cell_class1_blast_overlap.R $tmpfolder/$sagid'_check_dc_blast_res.tsv' $ref_spe_file $outfolder/$sagid'_overlap_info.txt' $outfolder

echo $sagid done
#mv $tmpfolder 'temp_'$sagid
rm *.log
#rm -r $tmpfolder *.log
