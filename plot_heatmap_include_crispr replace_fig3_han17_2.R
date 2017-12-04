# plot species level viral-host interaction heatmap
options(stringsAsFactors = F)
library(reshape2)
library(gplots)
source("sag_toolbox.R")

# read in matrix of viral signals and crispr signals
sag_vp_inf_m = read.table("matrix_infection.csv", header = T, sep = ",")
sag_vp_crispr_m = read.table("matrix_crispr_dup.csv", header = T, sep = ",")

rownames(sag_vp_inf_m) = substr(sag_vp_inf_m[,1], 8, 10)
sag_vp_inf_m = sag_vp_inf_m[,-1]
rownames(sag_vp_crispr_m) = sag_vp_crispr_m[,1]
sag_vp_crispr_m = sag_vp_crispr_m[,-1]

# collect all viral partitions number and sag id
all_vp = sort(gsub("X", "", unique(c(colnames(sag_vp_inf_m), colnames(sag_vp_crispr_m)))))
all_sag = sort(unique(c(rownames(sag_vp_inf_m), rownames(sag_vp_crispr_m))))

# change the matrices into list for further combination
sag_vp_inf_l = melt(as.matrix(sag_vp_inf_m))
sag_vp_inf_l[, 2] = as.numeric(gsub("X", "", sag_vp_inf_l[, 2]))
sag_vp_crispr_l = melt(as.matrix(sag_vp_crispr_m))
sag_vp_crispr_l[, 2] = as.numeric(gsub("X", "", sag_vp_crispr_l[, 2]))

# keep only the ones with either viral signals and crispr signals, exclude 0 in the matrices
sag_vp_inf_l_gt0 = sag_vp_inf_l[which(sag_vp_inf_l[,3] > 0), ]
sag_vp_inf_l_gt0[, ncol(sag_vp_inf_l_gt0) + 1] = paste0(sag_vp_inf_l_gt0[, 1], "+-*/", sag_vp_inf_l_gt0[, 2])
sag_vp_crispr_l_gt0 = sag_vp_crispr_l[which(sag_vp_crispr_l[,3] > 0), ]
sag_vp_crispr_l_gt0[, ncol(sag_vp_crispr_l_gt0) + 1] = paste0(sag_vp_crispr_l_gt0[, 1], "+-*/", sag_vp_crispr_l_gt0[,2])
pair_not_in_inf_idx = !(sag_vp_crispr_l_gt0[, 4] %in% sag_vp_inf_l_gt0[, 4])

sag_vp_crispr_l_gt0_add_pair = sag_vp_crispr_l_gt0[pair_not_in_inf_idx, ]

# combine viral signal and crispr signal into one list
comb_l = rbind(sag_vp_inf_l_gt0[, 1:2], sag_vp_crispr_l_gt0_add_pair[, 1:2])


# read in the sag host cell classification file and formalize the species names
sag_spe_raw = read.csv("SAG_species_classification_Jan17.csv", header = T)
sag_spe_raw[, 2] = gsub(" and ", " & ", sag_spe_raw[, 2])
sag_spe = data.frame(sag_id = substr(sag_spe_raw$SAG, 8, 10), spe = sag_spe_raw$Species.classification)
sag_name_changer = read.csv("species_Name_change_jan17.csv", header = F)
sag_spe[, 2] = sapply(sag_spe[, 2], function(x){
  if(x %in% sag_name_changer[, 1]){
    y = which(as.character(sag_name_changer[, 1]) == as.character(x))
    return(sag_name_changer[y, 2])
  }
  else{
    return(x)
  }
})

# merge the combined list with the sag classificaion to get species level signal list
comb_l_spe = merge(comb_l, sag_spe, by.x = "Var1", by.y = "sag_id")
comb_m_spe = dcast(comb_l_spe, spe ~ Var2, fill = 0)
unique_comb_l_spe = unique(comb_l_spe[, c(1,3)])

# remove duplicated records
rownames(comb_m_spe) = comb_m_spe[,1]
comb_m_spe = comb_m_spe[, -1]

merged_sub_data = as.matrix(comb_m_spe[-which(rownames(comb_m_spe)=="ND"),])

sag_spe_count = data.frame(table(sag_spe$spe))
sag_spe_count_match = sag_spe_count[which(sag_spe_count[, 1] %in% rownames(merged_sub_data)), ]
merged_sub_data_norm = merged_sub_data/matrix(rep(sag_spe_count_match[, 2], ncol(merged_sub_data)), nrow = nrow(sag_spe_count_match))


# get host species color for the heatmap
merged_sub_data_color = sapply(rownames(merged_sub_data),get_all_color_jan17)
hmcol2 = colorRampPalette(c("lightblue","darkblue"))(100)


# check to remove all zero columns or rows
merged_sub_data_norm_0_converted = merged_sub_data_norm
unknown_only_columns = which(colSums(merged_sub_data_norm_0_converted)==0)
merged_sub_data_norm_0_converted[merged_sub_data_norm_0_converted==0] = NA
merged_sub_data_norm_0_converted = merged_sub_data_norm_0_converted[,-unknown_only_columns]

merged_sub_data_norm_0_converted_colname_corrected = merged_sub_data_norm_0_converted


# change column order by num sag for the heatmap in a decresing order
sag_vp = read.csv("matrix_infection.csv", na.strings = 0)
sag_vp = melt(sag_vp, variable.name = "x", na.rm = T)
sag_vp[,2] = gsub("X", "", sag_vp[,2])
sag_vp[, ncol(sag_vp)+1] = substr(sag_vp[,1], 8, 10)
colnames(sag_vp) = c("sag_name", "vp", "num_reads", "sag_id")
vp_order = sort(table(sag_vp$vp), decreasing = T)
vp_order_name = names(vp_order)
vp_order_name_match = intersect(vp_order_name, colnames(merged_sub_data_norm_0_converted_colname_corrected))
vp_order_match = vp_order[vp_order_name_match]
vp_order_match_string = paste0(vp_order_name_match, " (", vp_order_match, ")")

# if there is a tie in the above order, then apply column order by number of host species that have the signal
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
column_sagcnt_cnt = column_sagcnt_cnt[intersect(names(column_hostcnt_cnt), names(column_sagcnt_cnt))]
column_comb_cnt = rbind(column_hostcnt_cnt, column_sagcnt_cnt)
column_comb_cnt_order = column_comb_cnt[, order(column_comb_cnt[1,],column_comb_cnt[2,])]

# the x-axis labels for the heatmap
column_comb_cnt_order_string = paste0(colnames(column_comb_cnt_order), " (", column_comb_cnt_order[1,], ", ", column_comb_cnt_order[2,], ")")


# change row orders based on phylogenetic information and formalized the heatmap row text labels
roworder = read.table("heatmap_row_order_jan17_2.txt", header=F, sep=",")
rearrange_order = as.vector(sapply(t(as.vector(roworder)), function(x){which(x==rownames(merged_sub_data))}))

merged_sub_data_norm_0_converted_colname_corrected_rearranged = merged_sub_data_norm_0_converted_colname_corrected[rearrange_order,]
merged_sub_data_norm_0_converted_colname_corrected_rearranged_colreorder = merged_sub_data_norm_0_converted_colname_corrected_rearranged[, vp_order_name_match]
colnames(merged_sub_data_norm_0_converted_colname_corrected_rearranged_colreorder) = vp_order_match_string

merged_sub_data_norm_0_converted_colname_corrected_rearranged_byhostcnt_colreorder = merged_sub_data_norm_0_converted_colname_corrected_rearranged[, hostcnt_column_order]
colnames(merged_sub_data_norm_0_converted_colname_corrected_rearranged_byhostcnt_colreorder) = vp_order_match_string_orderby_hostcnt

merged_sub_data_norm_0_converted_colname_corrected_rearranged_complex_colreorder = merged_sub_data_norm_0_converted_colname_corrected_rearranged[, colnames(column_comb_cnt_order)]
colnames(merged_sub_data_norm_0_converted_colname_corrected_rearranged_complex_colreorder) = column_comb_cnt_order_string

merged_sub_data_color_rearranged = merged_sub_data_color[rearrange_order]


# make the plot and save to file
pdf("./fig3/heatmap_species_level_no_unknown_combine_inf_crispr.pdf", height = 6, width = 10)

heatmap.2(merged_sub_data_norm_0_converted_colname_corrected_rearranged_complex_colreorder,
          Rowv = F,
          Colv = F,
          col = hmcol2,
          trace = "none", 
          scale = "none",
          labCol = NULL,
          density.info = "none",
          main = "",
          dendrogram = "column",
          RowSideColors = merged_sub_data_color_rearranged,
          na.color = "white")

dev.off()