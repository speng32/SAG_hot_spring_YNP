options(stringsAsFactors = F)
library(reshape2)
library(dplyr)
library(plyr)

source("sag_toolbox.R")

data=read.table("SAG_spacers_to_NL10_vDNA_1E-5_80%id_25bp_partition_num.csv",sep=",",header = F)

data_id = unique(data.frame(SAG_id = sapply(data[,1], function(x){substr(x, 9,11)}), seq = data[,1], partition = as.numeric(trans_vp_from_original_to_ben(data[, ncol(data)]))))

data_id_no_na = data_id[which(!is.na(data_id[,3])),]

rown = sort(unique(data_id_no_na[,1]))
coln = sort(unique(data_id_no_na[,3]))
one_hit_recs_matrix = data.frame(matrix(0, ncol = length(coln), nrow = length(rown)))
rownames(one_hit_recs_matrix) = rown
colnames(one_hit_recs_matrix) = coln

uniq_recs_row_idx = !(duplicated(data_id_no_na[,2]) | duplicated(data_id_no_na[,2], fromLast = TRUE))

uniq_recs_list = data_id_no_na[uniq_recs_row_idx, ]
dup_rec_list = data_id_no_na[!uniq_recs_row_idx, ]

dup_hit_recs_matrix = one_hit_recs_matrix

for (j in 1:nrow(uniq_recs_list)) {
  one_hit_recs_matrix[uniq_recs_list[j, 1], as.character(uniq_recs_list[j, 3])] = one_hit_recs_matrix[uniq_recs_list[j, 1], as.character(uniq_recs_list[j, 3])] + 1
}

dup_rec_list_add = data.frame(dup_rec_list, known_counts = apply(dup_rec_list, 1,  function(x){return(one_hit_recs_matrix[x[1], as.character(as.integer(x[3]))])}))

dup_rec_list_add_max = ddply(dup_rec_list_add, ~SAG_id+seq, function(x){x[which.max(x$known_counts),]})

for (i in 1:nrow(dup_rec_list_add_max)) {
  one_hit_recs_matrix[dup_rec_list_add_max[i, 1], as.character(dup_rec_list_add_max[i, 3])] = one_hit_recs_matrix[dup_rec_list_add_max[i, 1], as.character(dup_rec_list_add_max[i, 3])] + 1
}

dup_cnt_cnt = dcast(dup_rec_list_add, seq~known_counts)
dup_cnt_cnt_num_g0 = apply(dup_cnt_cnt, 1, function(x){sum(x>0)})
dup_cnt_cnt_tie_sags = unique(substr(dup_cnt_cnt[which(dup_cnt_cnt_num_g0==1), 1], 9,11))

spacer_vp_mapping = rbind(uniq_recs_list, dup_rec_list_add_max[,1:3])
spacer_vp_mapping_order = spacer_vp_mapping[order(spacer_vp_mapping[,1]),]

#write.csv(one_hit_recs_matrix, file="matrix_crispr.csv", quote=F, row.names = T)
#write.csv(spacer_vp_mapping_order, file="unique_crispr_vp_list.csv", quote=F, row.names = F)


#### consider duplicate ####

uniq_recs_list_indi = cbind(uniq_recs_list, single_vp="True")
dup_rec_list_indi = cbind(dup_rec_list, single_vp="False")

spacer_vp_mapping_order_dup = rbind(uniq_recs_list_indi, dup_rec_list_indi)

for (j in 1:nrow(uniq_recs_list)) {
  dup_hit_recs_matrix[uniq_recs_list[j, 1], as.character(uniq_recs_list[j, 3])] = dup_hit_recs_matrix[uniq_recs_list[j, 1], as.character(uniq_recs_list[j, 3])] + 1
}

for (i in 1:nrow(dup_rec_list_indi)) {
  dup_hit_recs_matrix[dup_rec_list_indi[i, 1], as.character(dup_rec_list_indi[i, 3])] = dup_hit_recs_matrix[dup_rec_list_indi[i, 1], as.character(dup_rec_list_indi[i, 3])] + 1
}

write.csv(dup_hit_recs_matrix, file="matrix_crispr_dup.csv", quote=F, row.names = T)
write.csv(spacer_vp_mapping_order_dup, file="dup_crispr_vp_list.csv", quote=F, row.names = F)



