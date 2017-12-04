options(stringsAsFactors = F)
library(reshape2)
library(gplots)
source("../../sag_toolbox.R")

sag_vp_inf_m = read.table("../matrix_infection.csv", header = T, sep = ",")
sag_vp_crispr_m = read.table("../matrix_crispr_dup.csv", header = T, sep = ",")

rownames(sag_vp_inf_m) = substr(sag_vp_inf_m[,1], 8, 10)
sag_vp_inf_m = sag_vp_inf_m[,-1]
rownames(sag_vp_crispr_m) = sag_vp_crispr_m[,1]
sag_vp_crispr_m = sag_vp_crispr_m[,-1]

all_vp = sort(gsub("X", "", unique(c(colnames(sag_vp_inf_m), colnames(sag_vp_crispr_m)))))
all_sag = sort(unique(c(rownames(sag_vp_inf_m), rownames(sag_vp_crispr_m))))

sag_vp_inf_l = melt(as.matrix(sag_vp_inf_m))
sag_vp_inf_l[, 2] = as.numeric(gsub("X", "", sag_vp_inf_l[, 2]))
sag_vp_crispr_l = melt(as.matrix(sag_vp_crispr_m))
sag_vp_crispr_l[, 2] = as.numeric(gsub("X", "", sag_vp_crispr_l[, 2]))

sag_vp_inf_l_gt0 = sag_vp_inf_l[which(sag_vp_inf_l[,3] > 0), ]
sag_vp_inf_l_gt0[, ncol(sag_vp_inf_l_gt0) + 1] = paste0(sag_vp_inf_l_gt0[, 1], "+-*/", sag_vp_inf_l_gt0[, 2])
sag_vp_crispr_l_gt0 = sag_vp_crispr_l[which(sag_vp_crispr_l[,3] > 0), ]
sag_vp_crispr_l_gt0[, ncol(sag_vp_crispr_l_gt0) + 1] = paste0(sag_vp_crispr_l_gt0[, 1], "+-*/", sag_vp_crispr_l_gt0[,2])
pair_not_in_inf_idx = !(sag_vp_crispr_l_gt0[, 4] %in% sag_vp_inf_l_gt0[, 4])

sag_vp_crispr_l_gt0_add_pair = sag_vp_crispr_l_gt0[pair_not_in_inf_idx, ]

comb_l = rbind(sag_vp_inf_l_gt0[, 1:2], sag_vp_crispr_l_gt0_add_pair[, 1:2])

sag_spe_raw = read.csv("../auto_classification_results.txt", header = F)
sag_spe_raw = sag_spe_raw[which(sag_spe_raw[, 2] < 2), ]

sag_name_changer = read.csv("../species_name_change_auto_classify.csv", header = F)
sag_spe_raw[, 3] = sapply(sag_spe_raw[, 3], function(x){
  if(x %in% sag_name_changer[, 1]){
    y = which(as.character(sag_name_changer[, 1]) == as.character(x))
    return(sag_name_changer[y, 2])
  }
  else{
    return(x)
  }
})
colnames(sag_spe_raw) = c("SAG", "class", "Species.classification")
sag_spe = data.frame(sag_id = substr(sag_spe_raw$SAG, 8, 10), spe = sag_spe_raw$Species.classification)

comb_l_spe = merge(comb_l, sag_spe, by.x = "Var1", by.y = "sag_id")
comb_m_spe = dcast(comb_l_spe, spe ~ Var2, fill = 0)

unique_comb_l_spe = unique(comb_l_spe[, c(1,3)])


rownames(comb_m_spe) = comb_m_spe[,1]
comb_m_spe = comb_m_spe[, -1]

merged_sub_data = as.matrix(comb_m_spe)

sag_spe_count = data.frame(table(sag_spe$spe))
sag_spe_count_match = sag_spe_count[which(sag_spe_count[, 1] %in% rownames(merged_sub_data)), ]

merged_sub_data_norm = merged_sub_data/matrix(rep(sag_spe_count_match[, 2], ncol(merged_sub_data)), nrow = nrow(sag_spe_count_match))

#merged_sub_data_norm = merged_sub_data_norm[, -which(colSums(merged_sub_data_norm)==0)]

merged_sub_data_color = sapply(rownames(merged_sub_data), get_all_color_auto_class)
hmcol2 = colorRampPalette(c("lightblue","darkblue"))(100)



merged_sub_data_norm_0_converted = merged_sub_data_norm
merged_sub_data_norm_0_converted[merged_sub_data_norm_0_converted==0] = NA

merged_sub_data_norm_0_converted_colname_corrected = merged_sub_data_norm_0_converted


# column order by num sag
sag_vp = read.csv("../matrix_infection.csv", na.strings = 0)
sag_vp = melt(sag_vp, variable.name = "x", na.rm = T)
sag_vp[,2] = gsub("X", "", sag_vp[,2])
sag_vp[, ncol(sag_vp)+1] = substr(sag_vp[,1], 8, 10)
colnames(sag_vp) = c("sag_name", "vp", "num_reads", "sag_id")
vp_order = sort(table(sag_vp$vp), decreasing = T)
vp_order_name = names(vp_order)
vp_order_name_match = intersect(vp_order_name, colnames(merged_sub_data_norm_0_converted_colname_corrected))
vp_order_match = vp_order[vp_order_name_match]
vp_order_match_string = paste0(vp_order_name_match, " (", vp_order_match, ")")

#column order by num host
merged_sub_data_norm_0_converted_colname_corrected_orderby_hostcnt = merged_sub_data_norm_0_converted_colname_corrected[,order(apply(merged_sub_data_norm_0_converted_colname_corrected, 2, function(x){sum(!is.na(x))}))]
hostcnt_column_order = order(apply(merged_sub_data_norm_0_converted_colname_corrected, 2, function(x){sum(!is.na(x))}))
hostcnt_column = apply(merged_sub_data_norm_0_converted_colname_corrected, 2, function(x){sum(!is.na(x))})[hostcnt_column_order]
vp_order_name_match_orderby_hostcnt = intersect(colnames(merged_sub_data_norm_0_converted_colname_corrected_orderby_hostcnt), vp_order_name)
vp_order_match_orderby_hostcnt = vp_order[vp_order_name_match_orderby_hostcnt]
vp_order_match_string_orderby_hostcnt = paste0(vp_order_name_match_orderby_hostcnt, " (", hostcnt_column, ', ', vp_order_match_orderby_hostcnt, ")")

#column complex order
#remove multi-species
no_multi_row_idx = grep(" & ", rownames(merged_sub_data_norm_0_converted_colname_corrected), invert = T)
merged_sub_data_norm_0_converted_colname_corrected_no_multi = merged_sub_data_norm_0_converted_colname_corrected[no_multi_row_idx, ]
column_hostcnt_cnt = apply(merged_sub_data_norm_0_converted_colname_corrected_no_multi, 2, function(x){sum(!is.na(x))})

column_sagcnt_cnt = table(comb_l_spe$Var2)
#column_sagcnt_cnt = table(sag_vp$vp)
column_sagcnt_cnt = column_sagcnt_cnt[intersect(names(column_hostcnt_cnt), names(column_sagcnt_cnt))]
column_comb_cnt = rbind(column_hostcnt_cnt, column_sagcnt_cnt)
column_comb_cnt_order = column_comb_cnt[, order(column_comb_cnt[1,],column_comb_cnt[2,])]

column_comb_cnt_order_string = paste0(colnames(column_comb_cnt_order), " (", column_comb_cnt_order[1,], ", ", column_comb_cnt_order[2,], ")")


roworder = read.table("../fig3_row_order_based_on_phylo_class1_only.txt", header=F, sep=",")

rearrange_order = as.vector(sapply(t(as.vector(roworder)), function(x){which(x==rownames(merged_sub_data))}))

merged_sub_data_norm_0_converted_colname_corrected_rearranged = merged_sub_data_norm_0_converted_colname_corrected[rearrange_order,]
merged_sub_data_norm_0_converted_colname_corrected_rearranged_colreorder = merged_sub_data_norm_0_converted_colname_corrected_rearranged[, vp_order_name_match]
colnames(merged_sub_data_norm_0_converted_colname_corrected_rearranged_colreorder) = vp_order_match_string

merged_sub_data_norm_0_converted_colname_corrected_rearranged_byhostcnt_colreorder = merged_sub_data_norm_0_converted_colname_corrected_rearranged[, hostcnt_column_order]
colnames(merged_sub_data_norm_0_converted_colname_corrected_rearranged_byhostcnt_colreorder) = vp_order_match_string_orderby_hostcnt

merged_sub_data_norm_0_converted_colname_corrected_rearranged_complex_colreorder = merged_sub_data_norm_0_converted_colname_corrected_rearranged[, colnames(column_comb_cnt_order)]
colnames(merged_sub_data_norm_0_converted_colname_corrected_rearranged_complex_colreorder) = column_comb_cnt_order_string

merged_sub_data_color_rearranged = merged_sub_data_color[rearrange_order]



pdf("fig3_class1_only_raw.pdf", height = 6, width = 10)

heatmap.2(merged_sub_data_norm_0_converted_colname_corrected_rearranged_complex_colreorder,
          Rowv = F,
          Colv = F,
          col = hmcol2,
          trace = "none", 
          scale = "none",
          labCol = NULL,
          density.info = "none",
          main = "",
          dendrogram = "none",
          RowSideColors = merged_sub_data_color_rearranged,
          na.color = "white")

dev.off()