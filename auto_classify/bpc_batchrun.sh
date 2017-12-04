#!/bin/bash

if [ ! -d BPC_batchrun ]; then
	mkdir BPC_batchrun
fi

Rscript ./code/get_ani_matrix_mod1_auto_classify.R /home/speng32/data/sags/ani_mod1/all_ani_output.tsv /home/speng32/data/sags/ani_mod1/Corrected_SAG_trimmed_contigs /home/speng32/data/sags/ani_mod1/Corrected_SAG_reference_genomes /home/speng32/data/sags/Reference_genomes_per_species.csv /home/speng32/data/sags/ref_ref_ani/ref_ref_background_overlap.csv df stringent 5 20 ./BPC_batchrun/original_bpc20-5_stringent
echo 20 done

Rscript ./code/get_ani_matrix_mod1_auto_classify.R /home/speng32/data/sags/ani_mod1/all_ani_output.tsv /home/speng32/data/sags/ani_mod1/Corrected_SAG_trimmed_contigs /home/speng32/data/sags/ani_mod1/Corrected_SAG_reference_genomes /home/speng32/data/sags/Reference_genomes_per_species.csv /home/speng32/data/sags/ref_ref_ani/ref_ref_background_overlap.csv df stringent 5 25 ./BPC_batchrun/original_bpc25-5_stringent
echo 25 done


Rscript ./code/get_ani_matrix_mod1_auto_classify.R /home/speng32/data/sags/ani_mod1/all_ani_output.tsv /home/speng32/data/sags/ani_mod1/Corrected_SAG_trimmed_contigs /home/speng32/data/sags/ani_mod1/Corrected_SAG_reference_genomes /home/speng32/data/sags/Reference_genomes_per_species.csv /home/speng32/data/sags/ref_ref_ani/ref_ref_background_overlap.csv df stringent 5 30 ./BPC_batchrun/original_bpc30-5_stringent
echo 30 done


Rscript ./code/get_ani_matrix_mod1_auto_classify.R /home/speng32/data/sags/ani_mod1/all_ani_output.tsv /home/speng32/data/sags/ani_mod1/Corrected_SAG_trimmed_contigs /home/speng32/data/sags/ani_mod1/Corrected_SAG_reference_genomes /home/speng32/data/sags/Reference_genomes_per_species.csv /home/speng32/data/sags/ref_ref_ani/ref_ref_background_overlap.csv df stringent 5 35 ./BPC_batchrun/original_bpc35-5_stringent
echo 35 done


Rscript ./code/get_ani_matrix_mod1_auto_classify.R /home/speng32/data/sags/ani_mod1/all_ani_output.tsv /home/speng32/data/sags/ani_mod1/Corrected_SAG_trimmed_contigs /home/speng32/data/sags/ani_mod1/Corrected_SAG_reference_genomes /home/speng32/data/sags/Reference_genomes_per_species.csv /home/speng32/data/sags/ref_ref_ani/ref_ref_background_overlap.csv df stringent 5 40 ./BPC_batchrun/original_bpc40-5_stringent
echo 40 done


Rscript ./code/get_ani_matrix_mod1_auto_classify.R /home/speng32/data/sags/ani_mod1/all_ani_output.tsv /home/speng32/data/sags/ani_mod1/Corrected_SAG_trimmed_contigs /home/speng32/data/sags/ani_mod1/Corrected_SAG_reference_genomes /home/speng32/data/sags/Reference_genomes_per_species.csv /home/speng32/data/sags/ref_ref_ani/ref_ref_background_overlap.csv df stringent 5 45 ./BPC_batchrun/original_bpc45-5_stringent
echo 45 done


Rscript ./code/get_ani_matrix_mod1_auto_classify.R /home/speng32/data/sags/ani_mod1/all_ani_output.tsv /home/speng32/data/sags/ani_mod1/Corrected_SAG_trimmed_contigs /home/speng32/data/sags/ani_mod1/Corrected_SAG_reference_genomes /home/speng32/data/sags/Reference_genomes_per_species.csv /home/speng32/data/sags/ref_ref_ani/ref_ref_background_overlap.csv df stringent 5 50 ./BPC_batchrun/original_bpc50-5_stringent
echo 50 done


Rscript ./code/get_ani_matrix_mod1_auto_classify.R /home/speng32/data/sags/ani_mod1/all_ani_output.tsv /home/speng32/data/sags/ani_mod1/Corrected_SAG_trimmed_contigs /home/speng32/data/sags/ani_mod1/Corrected_SAG_reference_genomes /home/speng32/data/sags/Reference_genomes_per_species.csv /home/speng32/data/sags/ref_ref_ani/ref_ref_background_overlap.csv df stringent 5 55 ./BPC_batchrun/original_bpc55-5_stringent
echo 55 done


Rscript ./code/get_ani_matrix_mod1_auto_classify.R /home/speng32/data/sags/ani_mod1/all_ani_output.tsv /home/speng32/data/sags/ani_mod1/Corrected_SAG_trimmed_contigs /home/speng32/data/sags/ani_mod1/Corrected_SAG_reference_genomes /home/speng32/data/sags/Reference_genomes_per_species.csv /home/speng32/data/sags/ref_ref_ani/ref_ref_background_overlap.csv df stringent 5 60 ./BPC_batchrun/original_bpc60-5_stringent
echo 60 done


Rscript ./code/get_ani_matrix_mod1_auto_classify.R /home/speng32/data/sags/ani_mod1/all_ani_output.tsv /home/speng32/data/sags/ani_mod1/Corrected_SAG_trimmed_contigs /home/speng32/data/sags/ani_mod1/Corrected_SAG_reference_genomes /home/speng32/data/sags/Reference_genomes_per_species.csv /home/speng32/data/sags/ref_ref_ani/ref_ref_background_overlap.csv df stringent 5 65 ./BPC_batchrun/original_bpc65-5_stringent
echo 65 done


Rscript ./code/get_ani_matrix_mod1_auto_classify.R /home/speng32/data/sags/ani_mod1/all_ani_output.tsv /home/speng32/data/sags/ani_mod1/Corrected_SAG_trimmed_contigs /home/speng32/data/sags/ani_mod1/Corrected_SAG_reference_genomes /home/speng32/data/sags/Reference_genomes_per_species.csv /home/speng32/data/sags/ref_ref_ani/ref_ref_background_overlap.csv df stringent 5 70 ./BPC_batchrun/original_bpc70-5_stringent
echo 70 done



echo 20
cut -d',' -f 2 ./BPC_batchrun/original_bpc20-5_stringent/auto_classification_results.txt | sort | uniq -c
echo

echo 25
cut -d',' -f 2 ./BPC_batchrun/original_bpc25-5_stringent/auto_classification_results.txt | sort | uniq -c
echo

echo 30
cut -d',' -f 2 ./BPC_batchrun/original_bpc30-5_stringent/auto_classification_results.txt | sort | uniq -c
echo

echo 35
cut -d',' -f 2 ./BPC_batchrun/original_bpc35-5_stringent/auto_classification_results.txt | sort | uniq -c
echo

echo 40
cut -d',' -f 2 ./BPC_batchrun/original_bpc40-5_stringent/auto_classification_results.txt | sort | uniq -c
echo

echo 45
cut -d',' -f 2 ./BPC_batchrun/original_bpc45-5_stringent/auto_classification_results.txt | sort | uniq -c
echo

echo 50
cut -d',' -f 2 ./BPC_batchrun/original_bpc50-5_stringent/auto_classification_results.txt | sort | uniq -c
echo

echo 55
cut -d',' -f 2 ./BPC_batchrun/original_bpc55-5_stringent/auto_classification_results.txt | sort | uniq -c
echo

echo 60
cut -d',' -f 2 ./BPC_batchrun/original_bpc60-5_stringent/auto_classification_results.txt | sort | uniq -c
echo

echo 65
cut -d',' -f 2 ./BPC_batchrun/original_bpc65-5_stringent/auto_classification_results.txt | sort | uniq -c
echo

echo 70
cut -d',' -f 2 ./BPC_batchrun/original_bpc70-5_stringent/auto_classification_results.txt | sort | uniq -c
echo
