calc_ref1_vs_ref2_overlap_stat <- function(inref, per_experiment_draw_times, sag_id, sag_overlap){
  ref1 = inref[1]
  ref2 = inref[2]

  ref1_v_ref2_blast_record = all_blast_record[which(all_blast_record$q_name == ref1 & all_blast_record$s_name == ref2), ]
  
  ref1_spe_header = ref_spe_header[which(ref_spe_header$ref_genome == ref1), ]
  #ref2_spe_header = ref_spe_header[which(ref_spe_header$ref_species == ref2), ]
  
  ref1_rand_index = t(matrix(sample.int(nrow(ref1_spe_header), per_experiment_draw_times * repeat_experiment_times, replace = T), nrow = per_experiment_draw_times, ncol = repeat_experiment_times))
  
  ref1_v_ref2_overlap = integer()
  
  for(i in 1:repeat_experiment_times){
    cur_ref1_spe_header = ref1_spe_header[ref1_rand_index[i, ], ]
    cur_ref1_spe_header_merge_blast_res = merge(cur_ref1_spe_header, ref1_v_ref2_blast_record, by.x = "split_header", by.y = "q_header")
    if(nrow(cur_ref1_spe_header_merge_blast_res) > 0){
      cur_overlap = sum(cur_ref1_spe_header_merge_blast_res$qend - cur_ref1_spe_header_merge_blast_res$qstart + 1)
      ref1_v_ref2_overlap = c(ref1_v_ref2_overlap, cur_overlap)
    }else{
      ref1_v_ref2_overlap = c(ref1_v_ref2_overlap, 0)
    }
  }

  cur_lower_sig_5_percent = quantile(ref1_v_ref2_overlap, 0.05)
  hist(ref1_v_ref2_overlap, breaks = 1000, main = paste(sag_id, "", ref1, "vs", ref2, col = "gray"))
  abline(v = sag_overlap, lty = 2, col = "red")
  abline(v = cur_lower_sig_5_percent, lty = 2, col = "blue")

  #return(ref1_v_ref2_overlap)
  return(sag_overlap < cur_lower_sig_5_percent)
}


pairwise_sig_test <- function(x, sag_ref_list_spe_table_l, per_experiment_draw_times, sag_id){
  ref_s_A = as.character(x[1])
  ref_s_B = as.character(x[2])
  ml1 = as.numeric(x[3])
  ml2 = as.numeric(x[4])
  mo = as.numeric(x[5])
  
  #print(ref_s_A)
  #print(ref_s_B)
  
  ref_g_As = sag_ref_list_spe_table_l[[ref_s_A]]$ref_genome
  ref_g_Bs = sag_ref_list_spe_table_l[[ref_s_B]]$ref_genome
  
  #print(length(ref_g_As))
  #print(length(ref_g_Bs))
  
  ref_ref_pairwise_table = rbind(data.frame(ref1 = ref_g_As, ref2 = ref_g_Bs), data.frame(ref1 = ref_g_Bs, ref2 = ref_g_As))
  
  ref_s_A_B_lower_level_sig = apply(ref_ref_pairwise_table, 1, calc_ref1_vs_ref2_overlap_stat, per_experiment_draw_times, sag_id, mo)
  #ref_s_A_B_pairwise_overlap = apply(ref_ref_pairwise_table, 1, calc_ref1_vs_ref2_overlap_stat, per_experiment_draw_times, sag_id, mo)
  #ref_s_A_B_pairwise_overlap_flat = as.vector(ref_s_A_B_pairwise_overlap)
  #lower_sig_5_percent = quantile(ref_s_A_B_pairwise_overlap_flat, 0.05)
  
  # background_hist = paste0("background_hist_", sag_id, ".pdf")
  # pdf(background_hist)
  # hist(ref_s_A_B_pairwise_overlap_flat, breaks = 1000, main = paste(sag_id, "with", length(ref_g_As), ref_s_A, "and", length(ref_g_Bs), ref_s_B), col = "gray")
  # abline(v = mo, lty = 2, col = "red")
  # abline(v = lower_sig_5_percent, lty = 2, col = "blue")
  # dev.off()
  # return(mo < lower_sig_5_percent)

  return(sum(ref_s_A_B_lower_level_sig) > nrow(ref_ref_pairwise_table)/2)
}

#double_cell_stat_test <- function(sag_ref_overlap_num, sag_ani_temp_folder, sag_ref_list_file, repeat_experiment_times = 200){
double_cell_stat_test <- function(y){
  #print(class(y))
  sag_ref_overlap_num = as.character(y[1])
  sag_ani_temp_folder = as.character(y[2])
  sag_ref_list_file = as.character(y[3])
  
  sag_id = gsub(".*/|_overlap_info.txt", "", sag_ref_overlap_num)
  
  #print(sag_ref_overlap_num)
  #print(sag_ani_temp_folder)
  #print(sag_ref_list_file)
  
  print(paste0("Dealing ", sag_ref_overlap_num))
  
  overlap_num = read.delim(sag_ref_overlap_num, header = F)
  colnames(overlap_num) = c("spe1", "spe2", "match_len1", "match_len2", "match_overlap")
  
  sag_query_split_file = system(paste0("find ", sag_ani_temp_folder, " -name Query.split"), intern = T)[1]
  sag_1020_chunk_num = as.integer(system(paste0("grep -c '^>' ", sag_query_split_file), intern = T))
  per_experiment_draw_times = sag_1020_chunk_num
  
  
  sag_ref_list = read.table(sag_ref_list_file, header = F)
  sag_ref_list[, ncol(sag_ref_list) + 1] = gsub(".*ref_|\\.fasta", "", sag_ref_list[, 1])
  colnames(sag_ref_list) = c("file_loc", "ref_genome")

  sag_ref_list_spe_table = merge(sag_ref_list, ref_spe_lookup, by.x = "ref_genome", by.y = "SAG1")
  ref_spe_names = unique(sag_ref_list_spe_table$Species)
  sag_ref_list_spe_table_l = split(sag_ref_list_spe_table, f = sag_ref_list_spe_table$Species)
  
  
  if(length(ref_spe_names) <= 3){
    background_hist = paste0("background_hist_", sag_id, "_genome_level.pdf")
    pdf(background_hist)
    sig_res = apply(overlap_num, 1, pairwise_sig_test, sag_ref_list_spe_table_l, per_experiment_draw_times, sag_id)
    dev.off()
  }else{
    stop("Error! More than 3 ref species")
  }
  
  return(any(sig_res))
}