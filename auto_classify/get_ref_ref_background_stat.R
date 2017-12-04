options(stringsAsFactors = F)
library(GenomicRanges)
library(MASS)
library(stringr)
library(seqinr)
library(reshape2)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

ref_location = args[1]
ref_genome_species_file = args[2]
ref_ref_ani_blast_res_file = args[3]
output_file = args[4]

# ref_location = "/home/speng32/data/sags/ani_mod1/Corrected_SAG_reference_genomes"
# ref_genome_species_file = "/home/speng32/data/sags/Reference_genomes_per_species.csv"
# ref_ref_ani_blast_res_file = "/home/speng32/data/sags/ref_ref_ani/all_ref_ref_ani_blast_result.tsv"
# output_file = "/home/speng32/data/sags/ref_ref_ani/ref_ref_background_overlap.csv"

all_files = list.files(path = ref_location)
ref_files = all_files[startsWith(all_files, "ref_")]

ref_file_lookup = data.frame(ref_id = gsub("ref_|.fasta", "", ref_files), ref_loc = paste0(ref_location, "/", ref_files), ref_file_name = gsub(".fasta", "", ref_files))

ref_genome_species_table = read.csv(ref_genome_species_file, header = T)
colnames(ref_genome_species_table) = c("ref_id", "ref_species")

ref_genome_species_loc = merge(ref_genome_species_table, ref_file_lookup, by = "ref_id")

get_genome_size <- function(x){
  rid = as.character(x[1])
  rspe = as.character(x[2])
  rloc = as.character(x[3])
  rfile = as.character(x[4])
  
  curf = read.fasta(rloc)
  bp = sum(sapply(curf, length))
  return(c(rid, bp))
}


ref_genome_size = data.frame(t(apply(ref_genome_species_loc, 1, get_genome_size)))
colnames(ref_genome_size) = c("ref_id", "ref_size")
ref_genome_size[, 2] = as.numeric(ref_genome_size[, 2])

ref_genome_species_loc_size = merge(ref_genome_species_loc, ref_genome_size, by = "ref_id")

#ref_genome_species_loc_size1 = ref_genome_species_loc_size
#colnames(ref_genome_species_loc_size1) = paste0(colnames(ref_genome_species_loc_size1), "_1")
#ref_genome_species_loc_size2 = ref_genome_species_loc_size
#colnames(ref_genome_species_loc_size2) = paste0(colnames(ref_genome_species_loc_size2), "_2")



ani_blast_res = read.delim(file = ref_ref_ani_blast_res_file, header = F)
colnames(ani_blast_res) = c("qid", "sid", "identity", "alignment_length", "mismatches", "gap", "qstart", "qend", "sstart", "send", "evalue", "bit_score", "qname", "sname")

#ani_blast_res_by_qname = split(ani_blast_res, f = ani_blast_res$qname)

ani_blast_res_by_qname = split(ani_blast_res, list(ani_blast_res$qname, ani_blast_res$sname))
blast_res_names = names(ani_blast_res_by_qname)

calc_query_interval <- function(n, x){
  query_start_stop = x[[n]][, c("qid", "qstart", "qend", "qname", "sname")]
  query_start_stop_order = query_start_stop[order(query_start_stop$qid, query_start_stop$qstart, query_start_stop$qend), ]
  
  split_name = strsplit(n, ".ref_")
  
  qname = split_name[[1]][1]
  sname = paste0("ref_", split_name[[1]][2])

  if(nrow(query_start_stop) == 0){
    return(c(qname, sname, 0))
  }

  if(length(unique(query_start_stop_order$qname)) != 1 || length(unique(query_start_stop_order$sname)) != 1){
    stop("ERROR")
  }
  
  total_len = 0
  
  cur_id = query_start_stop_order[1, 1]
  cur_start = query_start_stop_order[1, 2]
  cur_stop = query_start_stop_order[1, 3]
  
  if(nrow(query_start_stop) == 1){
    len = cur_stop - cur_start + 1
    return(c(qname, sname, len))
  }
  
  for(i in 2:nrow(query_start_stop_order)){
    iid = query_start_stop_order[i, 1]
    istart = query_start_stop_order[i, 2]
    istop = query_start_stop_order[i, 3]
    if(iid == cur_id){
      if(istart > cur_stop){
        total_len = total_len + cur_stop - cur_start + 1
        cur_start = istart
        cur_stop = istop
      }
      else{
        cur_stop = istop
      }
    }
    else{
      total_len = total_len + cur_stop - cur_start + 1
      cur_id = iid
      cur_start = istart
      cur_stop = istop
    }
  }
  total_len = total_len + cur_stop - cur_start + 1
  
  return(c(qname, sname, total_len))
}

ref_ref_bp_covered = data.frame(t(sapply(blast_res_names, calc_query_interval, ani_blast_res_by_qname)))
ref_ref_bp_covered[, 3] = as.numeric(ref_ref_bp_covered[, 3])
colnames(ref_ref_bp_covered) = c("q_ref_name", "s_ref_name", "bp1v2")

#ref_ref_bp_covered[, 4] = rownames(ref_ref_bp_covered)
#colnames(ref_ref_bp_covered) = c("ref_file_name_1", "ref_file_name_2", "bp1v2", "vs_name")

ref_genome_species_loc_size_bp1v2 = merge(ref_genome_species_loc_size, ref_ref_bp_covered, by.x = "ref_file_name", by.y = "q_ref_name")

s_info_table = ref_genome_species_loc_size[, c("ref_id", "ref_file_name", "ref_species", "ref_size")]
colnames(s_info_table) = c("s_ref_id", "s_ref_file_name", "s_ref_species", "s_ref_size")

ref_genome_species_loc_size_bp1v2_add_sspe = merge(ref_genome_species_loc_size_bp1v2, s_info_table, by.x = "s_ref_name", by.y = "s_ref_file_name")

#### at species level, no ref_A vs ref_A
none_self_mapping_rows = which(ref_genome_species_loc_size_bp1v2_add_sspe$s_ref_species != ref_genome_species_loc_size_bp1v2_add_sspe$ref_species)
ref_genome_species_loc_size_bp1v2_add_sspe_none = ref_genome_species_loc_size_bp1v2_add_sspe[none_self_mapping_rows, ]

sub_full_table = ref_genome_species_loc_size_bp1v2_add_sspe_none[, c("ref_species", "s_ref_species", "ref_size", "s_ref_size", "bp1v2")]
#sub_full_table[, ncol(sub_full_table)+1] = sub_full_table$ref_size + sub_full_table$s_ref_size
#colnames(sub_full_table) = c("ref_species", "s_ref_species", "ref_size", "s_ref_size", "bp1v2", "total_size")


#sub_full_table_order = sub_full_table[order(sub_full_table$ref_species), ]

ref_ref_total_size_mean = dcast(sub_full_table, ref_species ~ s_ref_species, value.var = "ref_size", mean)
ref_ref_total_size_mean_l = melt(ref_ref_total_size_mean, value.name="ref_size_mean", na.rm = T)
ref_ref_total_overlap_mean = dcast(sub_full_table, ref_species ~ s_ref_species, value.var = "bp1v2", mean)
ref_ref_total_overlap_mean_l = melt(ref_ref_total_overlap_mean, value.name = "overlap_mean", na.rm = T)

ref_ref_size_overlap_ratio = cbind(ref_ref_total_size_mean_l, ref_ref_total_overlap_mean_l[, 3])
ref_ref_size_overlap_ratio[, ncol(ref_ref_size_overlap_ratio) + 1] = ref_ref_size_overlap_ratio[, 4]/ref_ref_size_overlap_ratio[, 3]
colnames(ref_ref_size_overlap_ratio) = c("ref_species", "s_ref_species", "ref_size_mean", "overlap_mean", "ratio")

s_ref_mean_size = unique(ref_ref_size_overlap_ratio[, c(1, 3)])
colnames(s_ref_mean_size) = c("s_ref_species", "s_ref_size_mean")

ref_ref_both_size_overlap_ratio = merge(ref_ref_size_overlap_ratio[, 1:4], s_ref_mean_size, by = "s_ref_species")

reorder_ref_ref_size_table <- function(x){
  x1 = x[1]
  x2 = x[2]
  x3 = x[3]
  x4 = x[4]
  x5 = x[5]
  
  if(x1 > x2){
    return(c(x2, x1, x3, x5, x4))
  }
  else{
    return(c(x1, x2, x5, x3, x4))
  }
}

ref_ref_both_size_overlap_ratio_reorder = data.frame(t(apply(ref_ref_both_size_overlap_ratio, 1, reorder_ref_ref_size_table)))
colnames(ref_ref_both_size_overlap_ratio_reorder) = c("ref_s1", "ref_s2", "ref_s1_mean_size", "ref_s2_mean_size", "overlap")
#any(ref_ref_both_size_overlap_ratio_reorder$overlap > ref_ref_both_size_overlap_ratio_reorder$ref_s1_mean_size)

ref_ref_both_size_overlap_ratio_reorder_tokened = data.table(id = paste(ref_ref_both_size_overlap_ratio_reorder[, 1], ref_ref_both_size_overlap_ratio_reorder[, 2], ref_ref_both_size_overlap_ratio_reorder[, 3], ref_ref_both_size_overlap_ratio_reorder[, 4], sep = "---"), overlap = as.numeric(ref_ref_both_size_overlap_ratio_reorder[, 5]))

ref_ref_both_size_overlap_ratio_reorder_tokened[, mean_overlap := mean(overlap), by = id]
ref_ref_both_size_overlap_ratio_reorder_tokened_uniq = as.data.frame(ref_ref_both_size_overlap_ratio_reorder_tokened[!duplicated(ref_ref_both_size_overlap_ratio_reorder_tokened$id), c(1, 3), with = F])

ref_ref_background_table = data.frame(str_split(ref_ref_both_size_overlap_ratio_reorder_tokened_uniq$id, "---", simplify = T), ref_ref_both_size_overlap_ratio_reorder_tokened_uniq$mean_overlap)
colnames(ref_ref_background_table) = c("ref_s1", "ref_s2", "ref_s1_mean_size", "ref_s2_mean_size", "mean_overlap")
ref_ref_background_table[, 3] = as.numeric(ref_ref_background_table[, 3])
ref_ref_background_table[, 4] = as.numeric(ref_ref_background_table[, 4])
#any(ref_ref_background_table$mean_overlap > ref_ref_background_table$ref_s1_mean_size)

write.csv(ref_ref_background_table, file = output_file, quote = F, row.names = F)

#head(ref_ref_background_table[order(ref_ref_background_table$ratio, decreasing = T), ])
