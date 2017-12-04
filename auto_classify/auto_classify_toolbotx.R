auto_host_cell_classification <- function(one_sag_ani_bpc, bpc_single_threshold){
  sag_id = unique(one_sag_ani_bpc$SAG)
  sag_id_code = str_extract(sag_id, "[A-P][0-9][0-9]")
  if(length(sag_id) != 1){
    return(0)
  }

  cur_tmp_folder_name = paste0(tmp_folder_name, "_", sag_id_code)
  
  sag_record_ani91 = one_sag_ani_bpc[which(one_sag_ani_bpc$ANI >= 91 & one_sag_ani_bpc$BPC > 0), ]
  sag_record_ani95 = one_sag_ani_bpc[which(one_sag_ani_bpc$ANI >= 95 & one_sag_ani_bpc$BPC > 0), ]
  sag_record_ani95_bpc30 = one_sag_ani_bpc[which(one_sag_ani_bpc$ANI >= 95 & one_sag_ani_bpc$BPC >= bpc_single_threshold), ]
  
  if(nrow(sag_record_ani91) > 0 & length(unique(sag_record_ani91$ref_species_name)) >= 2){
    related_file_name_list = paste0(ref_file_location, "/ref_", unique(sag_record_ani91$ref_genome_name), ".fasta")
    write.table(related_file_name_list, file = paste0(output_folder, "/", sag_id_code, "_list.txt"), row.names = F, col.names = F, quote = F)
    system(paste0("/home/speng32/data/sags/code/check_double_cell_class1.sh ", sag_file_lookup[which(sag_file_lookup[, 1] == sag_id_code), 2], " ", paste0(sag_id_code, "_list.txt "), ref_genome_species_file, " ", cur_tmp_folder_name, " ", output_folder))
    overlap_num = read.delim(paste0(output_folder, "/", sag_id_code, "_overlap_info.txt"), header = F)
    colnames(overlap_num) = c("spe1", "spe2", "match_len1", "match_len2", "match_overlap")

    overlap_num[, 6] = overlap_num[, 5]/overlap_num[, 3]
    overlap_num[, 7] = overlap_num[, 5]/overlap_num[, 4]
    overlap_num[, 8] = apply(overlap_num, 1, function(x){if(x[6] <= 0.3 & x[7] <= 0.3){return(T)}else{return(F)}})
    
    if(any(overlap_num[, 8] == T)){
     return(c(2, paste0(unique(sag_record_ani91$ref_species_name), collapse = " & ")))
    }
    else if(mode == "stringent"){
      return(c(2.1, paste0(unique(sag_record_ani91$ref_species_name), collapse = " & ")))
    }
  }
  
  if(nrow(sag_record_ani95_bpc30) > 0 & length(unique(sag_record_ani95_bpc30$ref_species_name)) == 1){
    return(c(1, unique(sag_record_ani95_bpc30$ref_species_name)))
  }
  else if(length(unique(sag_record_ani95$ref_species_name)) == 1){
    return(c(3, paste0("*", unique(sag_record_ani95$ref_species_name))))
  }
  else{
    return(c(4, "no_auto_class"))
  }
}

auto_host_cell_classification_single_first <- function(one_sag_ani_bpc, bpc_single_threshold){
  sag_id = unique(one_sag_ani_bpc$SAG)
  sag_id_code = str_extract(sag_id, "[A-P][0-9][0-9]")
  if(length(sag_id) != 1){
    return(0)
  }
  
  sag_record_ani91 = one_sag_ani_bpc[which(one_sag_ani_bpc$ANI >= 91 & one_sag_ani_bpc$BPC > 0), ]
  sag_record_ani95 = one_sag_ani_bpc[which(one_sag_ani_bpc$ANI >= 95 & one_sag_ani_bpc$BPC > 0), ]
  sag_record_ani95_bpc30 = one_sag_ani_bpc[which(one_sag_ani_bpc$ANI >= 95 & one_sag_ani_bpc$BPC >= bpc_single_threshold), ]
  
  if(nrow(sag_record_ani95_bpc30) > 0 & length(unique(sag_record_ani95_bpc30$ref_species_name)) == 1){
    return(c(1, unique(sag_record_ani95_bpc30$ref_species_name)))
  }
  else if(nrow(sag_record_ani91) > 0 & length(unique(sag_record_ani91$ref_species_name)) >= 2){
    related_file_name_list = paste0(ref_file_location, "/ref_", unique(sag_record_ani91$ref_genome_name), ".fasta")
    write.table(related_file_name_list, file = paste0(output_folder, "/", sag_id_code, "_list.txt"), row.names = F, col.names = F, quote = F)
    system(paste0("/home/speng32/data/sags/code/check_double_cell_class1.sh ", sag_file_lookup[which(sag_file_lookup[, 1] == sag_id_code), 2], " ", paste0(sag_id_code, "_list.txt "), ref_genome_species_file, " ", tmp_folder_name, " ", output_folder))
    overlap_num = read.csv(paste0(sag_id_code, "_overlap_info.txt"), header = F)
    colnames(overlap_num) = c("spe1", "spe2", "match_len1", "match_len2", "match_overlap")

    overlap_num[, 6] = overlap_num[, 5]/overlap_num[, 3]
    overlap_num[, 7] = overlap_num[, 5]/overlap_num[, 4]
    overlap_num[, 8] = apply(overlap_num, 1, function(x){if(x[6] <= 0.3 & x[7] <= 0.3){return(T)}else{return(F)}})
    
    if(any(overlap_num[, 8] == T)){
     return(c(2, paste0(unique(sag_record_ani91$ref_species_name), collapse = " & ")))
    }
    else if(mode == "stringent"){
      return(c(2.1, paste0(unique(sag_record_ani91$ref_species_name), collapse = " & ")))
    }
  }

  if(length(unique(sag_record_ani95$ref_species_name)) == 1){
    return(c(3, paste0("*", unique(sag_record_ani95$ref_species_name))))
  }
  else{
    return(c(4, "no_auto_class"))
  }
}

get_vs_plots_in_species_in_class <- function(x, ani_m, bpc_m, classn){
  sag_id = gsub(".*_", "", rownames(ani_m)[x])
  ref_names = colnames(ani_m)
  ref_names_spe_add = sapply(ref_names, function(x){return(ref_genome_species[which(ref_genome_species[, 1]==x), ])})
  plot(as.numeric(ani_m[x, ]), as.numeric(bpc_m[x, ]), pch=20, xlab = "ANI in percentage", ylab = "bp used in percentage", ylim = c(0,100), xlim = c(0, 100), main = sag_id, col = rgb(0,0,0,0.5))
  text(as.numeric(ani_m[x, ]), as.numeric(bpc_m[x, ]), labels=ref_names_spe_add[2, ], cex=0.4)
  if(classn == 1){
    abline(v = 95, h = 30, col = "red", lty = 2)
  }
  else if(classn == 2){
    abline(v = 91, col = "blue", lty = 2)
  }
  else if(classn == 3){
    abline(v = 95, h = 30, col = "green", lty = 2)
  }
  return(T)
}

determine_significant_chi_square_test <- function(sag_overlap_info){
  if(nrow(sag_overlap_info) > 1){
    multiple_res = apply(sag_overlap_info, 1, determine_significant_chi_square_test)
  }
  else{
    sag_two = c(sag_overlap_info[1, 3] + sag_overlap_info[1, 4], sag_overlap_info[1, 5])
    if(sag_overlap_info[1, 1] > sag_overlap_info[1, 2]){
      back_match_row = which(ref_ref_background[, 2] == sag_overlap_info[1, 1] & ref_ref_background[, 1] == sag_overlap_info[1, 2])
    }
    else{
      back_match_row = which(ref_ref_background[, 1] == sag_overlap_info[1, 1] & ref_ref_background[, 2] == sag_overlap_info[1, 2])
    }
    back_two = c(ref_ref_background[back_match_row, 3] + ref_ref_background[back_match_row, 4], ref_ref_background[back_match_row, 5])

    back_two_norm = back_two/sum(back_two)

    chisq_res = chisq_test(x = sag_two, p = back_two_norm)
  }
}