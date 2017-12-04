options(stringsAsFactors = F)

# source("determine_sag_double_cell_sig_toolbox.R")
# args = commandArgs(trailingOnly=TRUE)
# 
# sag_ref_overlap_num = args[1]
# sag_ani_temp_folder = args[2]
# sag_ref_list_file = args[3]
# all_blast_record_file = args[4]
# ref_spe_header_ref_file = args[5]
# ref_spe_file = args[6]
# repeat_experiment_times = args[7]
#outfile = args[8]


#sag_ref_overlap_num = "/home/speng32/data/sags/test_stat/L13_overlap_info.txt"
#sag_ani_temp_folder = "/home/speng32/data/sags/test_stat/temp_L13"
#sag_ref_list_file = "/home/speng32/data/sags/test_stat/L13_list.txt"
all_blast_record_file = "/home/speng32/data/sags/ref_ref_ani/all_ref_ref_ani_blast_result_best_hit_only.tsv"
ref_spe_header_ref_file = "/home/speng32/data/sags/ref_ref_ani/ref_split_names_add_spe.txt"
ref_spe_file = "/home/speng32/data/sags/Reference_genomes_per_species.csv"
repeat_experiment_times = 200
#outfile = "/home/speng32/data/sags/ref_ref_ani/ref_background.txt"

all_blast_record = read.delim(all_blast_record_file, header = F)
colnames(all_blast_record) = c("q_header", "s_header", "identity", "alignment_length", "mismatches", "gap", "qstart", "qend", "sstart", "send", "evalue", "bit_score", "q_name", "s_name")
all_blast_record[, 13] = gsub("ref_", "", all_blast_record[, 13])
all_blast_record[, 14] = gsub("ref_", "", all_blast_record[, 14])

ref_spe_header = read.delim(ref_spe_header_ref_file, header = T)

ref_spe_lookup = read.csv(ref_spe_file, header = T)


all_temp_test = read.csv("/home/speng32/data/sags/test_stat/all_double_test.csv", header = F)
all_temp_test = as.data.frame(lapply(all_temp_test, function(x){paste0("/home/speng32/data/sags/test_stat/", x)}))

all_stat_res = apply(all_temp_test, 1, double_cell_stat_test)

# table(all_stat_res)
# all_temp_test[which(all_stat_res==T), 1]




#stat_test_result = double_cell_stat_test(sag_ref_overlap_num, sag_ani_temp_folder, sag_ref_list_file, 200)









# overlap_num = read.delim(sag_ref_overlap_num, header = F)
# colnames(overlap_num) = c("spe1", "spe2", "match_len1", "match_len2", "match_overlap")
# 
# sag_query_split_file = system(paste0("find ", sag_ani_temp_folder, " -name Query.split"), intern = T)[1]
# sag_1020_chunk_num = as.integer(system(paste0("grep -c '^>' ", sag_query_split_file), intern = T))
# per_experiment_draw_times = sag_1020_chunk_num
# 
# 
# sag_ref_list = read.table(sag_ref_list_file, header = F)
# sag_ref_list[, ncol(sag_ref_list) + 1] = gsub(".*ref_|\\.fasta", "", sag_ref_list[, 1])
# colnames(sag_ref_list) = c("file_loc", "ref_genome")
# ref_spe_lookup = read.csv(ref_spe_file, header = T)
# sag_ref_list_spe_table = merge(sag_ref_list, ref_spe_lookup, by.x = "ref_genome", by.y = "SAG1")
# ref_spe_names = unique(sag_ref_list_spe_table$Species)
# sag_ref_list_spe_table_l = split(sag_ref_list_spe_table, f = sag_ref_list_spe_table$Species)
# 
# 
# all_blast_record = read.delim(all_blast_record_file, header = F)
# colnames(all_blast_record) = c("q_header", "s_header", "identity", "alignment_length", "mismatches", "gap", "qstart", "qend", "sstart", "send", "evalue", "bit_score", "q_name", "s_name")
# all_blast_record[, 13] = gsub("ref_", "", all_blast_record[, 13])
# all_blast_record[, 14] = gsub("ref_", "", all_blast_record[, 14])
# 
# ref_spe_header = read.delim(ref_spe_header_ref_file, header = T)
# 
# if(length(ref_spe_names) <= 3){
#   sig_res = apply(overlap_num, 1, pairwise_sig_test)
# }else{
#   stop("Error! More than 3 ref species")
# }













