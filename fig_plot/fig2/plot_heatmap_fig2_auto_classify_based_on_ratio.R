options(stringsAsFactors = F)
library(reshape2)
library(gplots)
library(RColorBrewer)
library(data.table)

source("../../sag_toolbox.R")

#sag_spe_raw = read.csv("../auto_classification_results_manual_curate_C07_L21.txt", header = F)
#sag_spe_raw = read.csv("../auto_classification_results_stringent_manual_curate_C07_L21.txt", header = F)
sag_spe_raw = read.csv("../auto_classification_results.txt", header = F)
sag_spe_raw = sag_spe_raw[which(sag_spe_raw[, 2] <= 2), ]

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
roworder = as.character(read.table("../fig2_row_order_based_on_phylo.txt", header=F, sep=",")[, 1])
spe_ordered = ordered(sag_spe[, 2], levels = roworder)
sag_spe_ordered = sag_spe[order(spe_ordered), ]

#sag_vp_inf_m = read.table("../matrix_infection.csv", header = T, sep = ",")
sag_vp_inf_m = read.table("../matrix_infection_by_bp.csv", header = T, sep = ",")
sag_vp_crispr_m = read.table("../matrix_crispr.csv", header = T, sep = ",")

#######
# sort columns

sortCols = function(data){
  data[, order(colSums(data), decreasing = T)]
}

sortRows <- function(data){
  order_string = paste("data[order(", paste(paste0("data[,", 1:ncol(data), "]"), collapse=","), ", decreasing = T), ]", collapse = "", sep = "")
  data = eval(parse(text = order_string))
}

sortAllRows <- function(data_m, row_spe){
  ordered_m = NULL
  for(i in unique(row_spe[, 2])){
    #print(i)
    i_idx = which(row_spe[, 2] == i)
    if(length(i_idx) > 1){
      ordered_m = rbind(ordered_m, sortRows(data_m[i_idx, ]))
    }else{
      rname = rownames(ordered_m)
      ordered_m = rbind(ordered_m, data_m[i_idx, ])
      rownames(ordered_m) = c(rname, rownames(data_m)[i_idx])
    }
  }
  return(ordered_m)
}


sortAllRows_order_by_heatmap_cnt <- function(data_m, row_spe){
  ordered_m = NULL
  for(i in names(sort(table(row_spe[,2]), decreasing=T))){
    #print(i)
    i_idx = which(row_spe[, 2] == i)
    if(length(i_idx) > 1){
      ordered_m = rbind(ordered_m, sortRows(data_m[i_idx, ]))
    }else{
      rname = rownames(ordered_m)
      ordered_m = rbind(ordered_m, data_m[i_idx, ])
      rownames(ordered_m) = c(rname, rownames(data_m)[i_idx])
    }
  }
  return(ordered_m)
}


####################################### bp lower 2 reads 300bp #######################################

sag_vp_reads_m = sag_vp_inf_m
sag_vp_reads_m[sag_vp_reads_m < 300] = 0
rownames(sag_vp_reads_m) = sag_vp_reads_m[, 1]
sag_vp_reads_m = sag_vp_reads_m[, -1]

sag_vp_reads_m = ifelse(sag_vp_reads_m>0, 1, 0)
sag_vp_reads_m[is.na(sag_vp_reads_m)] = 0

sag_spe_ordered_sub_reads = sag_spe_ordered[which(sag_spe_ordered$sag_id %in% row.names(sag_vp_reads_m)),]
sag_vp_reads_m = sag_vp_reads_m[sag_spe_ordered_sub_reads$sag_id, ]

keep_rows = rowSums(sag_vp_reads_m) > 0
keep_cols = colSums(sag_vp_reads_m) > 0
sag_vp_reads_m = sag_vp_reads_m[keep_rows, keep_cols]

#write.csv(sag_spe_ordered_sub_reads, file="fig2_row_order_auto_classify.csv", row.names = F, quote = F)

#sag_spe_ordered_new_order = sag_spe_ordered
#sag_spe_ordered_sub_reads_new_order = sag_spe_ordered[which(sag_spe_ordered$sag_id %in% row.names(sag_vp_reads_m)), ]
#sag_vp_reads_m_new_order = sag_vp_reads_m[sag_spe_ordered_sub_reads_new_order$sag_id, ]

np_reads_sag_reads_colsorted = sortCols(sag_vp_reads_m)
colnames(np_reads_sag_reads_colsorted) = gsub("X", "", colnames(np_reads_sag_reads_colsorted))
sub_sag_spe_ordered = sag_spe_ordered[which(sag_spe_ordered[, 1] %in% rownames(np_reads_sag_reads_colsorted)), ]
np_reads_sag_reads_colsorted = sortAllRows(np_reads_sag_reads_colsorted, sub_sag_spe_ordered)

np_reads_sag_reads_colsorted_no_unk = np_reads_sag_reads_colsorted
sag_spe_ordered_no_unk = sag_spe_ordered[which(sag_spe_ordered$spe != "Unknown"), ]

matrix_sag_order = rownames(np_reads_sag_reads_colsorted_no_unk)
nature_order = sag_spe_ordered_no_unk[, 1]

common_sags = intersect(matrix_sag_order, nature_order)
leftout_sags = setdiff(nature_order, common_sags)

np_reads_sag_spe_reordered_sub_no_unk = sag_spe_ordered_no_unk[which(sag_spe_ordered_no_unk$sag_id %in% common_sags), ]
np_reads_sag_spe_reordered_sub_no_unk = np_reads_sag_spe_reordered_sub_no_unk[match(matrix_sag_order, np_reads_sag_spe_reordered_sub_no_unk$sag_id), ]

np_row_color_no_unk = sapply(np_reads_sag_spe_reordered_sub_no_unk[, 2], get_all_color_auto_class)

np_row_color_legend_c_no_unk = c(unique(np_row_color_no_unk), "#ffffff")
np_row_color_legend_n_no_unk = unique(np_reads_sag_spe_reordered_sub_no_unk[,2])

sag_spe_cnt = table(sub_sag_spe_ordered[,2])
final_species_with_hits = sag_spe_cnt[unique(sub_sag_spe_ordered[,2])]
legend_labels = c(paste(names(final_species_with_hits), " (", final_species_with_hits, ")", sep = ""), paste0("Total # of SAGs: ", nrow(sag_vp_reads_m)))

np_reads_sag_reads_colsorted_no_unk = np_reads_sag_reads_colsorted_no_unk[, -ncol(np_reads_sag_reads_colsorted_no_unk)]

pdf("fig2_lower_1.pdf", height = 12, width = 8)
par(las=1, cex = 1.2, mar=c(2,5,2,20))
image(t(as.matrix(rev(1:length(np_row_color_legend_c_no_unk)))), col = np_row_color_legend_c_no_unk, xlab = "", xaxt = 'n', yaxt = 'n', bty = 'n')
axis(4, at=seq(0, 1, by = 1/(length(legend_labels)-1)), labels=rev(legend_labels), tick = F)
dev.off()

col = brewer.pal(9, "Set1")[1:3]
heatmap_col = colorRampPalette(c("#F0F0F0", col[2]))(100)

pdf("fig2_lower_2.pdf", height = 10, width = 12)
heatmap.2(as.matrix(np_reads_sag_reads_colsorted_no_unk),Rowv = F,Colv = F,col = heatmap_col,trace = "none", scale ="none",labCol = NULL,density.info = "none",main = "",dendrogram ="none",na.color = "white",key = F,cexRow = 0.25,cexCol= 1, RowSideColors = np_row_color_no_unk, labRow = F)
dev.off()

write.csv(np_reads_sag_spe_reordered_sub_no_unk, file="fig2_row_order_auto_classify.csv", row.names = F, quote = F)

pdf("fig2_lower_3.pdf", height = 4, width = 9.4)
par(cex=1, las = 2)
plot(colSums(np_reads_sag_reads_colsorted_no_unk), pch=20, type="l", col="black", bty = 'n', xlab = "", ylab = "", lwd = 2, xaxt = 'n', ylim = c(0, 120))
#axis(side = 2, at = c(0, 20, 40, 60, 80, 100, 120), labels = c(0, 20, 40, 60, 80, 100, 120))
title(xlab = "Viral Partition Number", line=0)
title(ylab = "# of SAGs", line = 2.2)
points(x = 1:ncol(np_reads_sag_reads_colsorted_no_unk), y=colSums(np_reads_sag_reads_colsorted_no_unk), type = "p", pch = 20)
dev.off()
 



####################################### bp upper 5 reads 750bp #######################################

sag_vp_reads_m = sag_vp_inf_m
sag_vp_reads_m[sag_vp_reads_m < 750] = 0
rownames(sag_vp_reads_m) = sag_vp_reads_m[, 1]
sag_vp_reads_m = sag_vp_reads_m[, -1]

sag_vp_reads_m = ifelse(sag_vp_reads_m>0, 1, 0)
sag_vp_reads_m[is.na(sag_vp_reads_m)] = 0

sag_spe_ordered_sub_reads = sag_spe_ordered[which(sag_spe_ordered$sag_id %in% row.names(sag_vp_reads_m)),]
sag_vp_reads_m = sag_vp_reads_m[sag_spe_ordered_sub_reads$sag_id, ]

keep_rows = rowSums(sag_vp_reads_m) > 0
keep_cols = colSums(sag_vp_reads_m) > 0
sag_vp_reads_m = sag_vp_reads_m[keep_rows, keep_cols]
#write.csv(sag_spe_ordered_sub_reads, file="fig2_row_order_auto_classify.csv", row.names = F, quote = F)

#sag_spe_ordered_new_order = sag_spe_ordered
#sag_spe_ordered_sub_reads_new_order = sag_spe_ordered[which(sag_spe_ordered$sag_id %in% row.names(sag_vp_reads_m)), ]
#sag_vp_reads_m_new_order = sag_vp_reads_m[sag_spe_ordered_sub_reads_new_order$sag_id, ]

np_reads_sag_reads_colsorted = sortCols(sag_vp_reads_m)
colnames(np_reads_sag_reads_colsorted) = gsub("X", "", colnames(np_reads_sag_reads_colsorted))
sub_sag_spe_ordered = sag_spe_ordered[which(sag_spe_ordered[, 1] %in% rownames(np_reads_sag_reads_colsorted)), ]
np_reads_sag_reads_colsorted = sortAllRows(np_reads_sag_reads_colsorted, sub_sag_spe_ordered)

np_reads_sag_reads_colsorted_no_unk = np_reads_sag_reads_colsorted
sag_spe_ordered_no_unk = sag_spe_ordered[which(sag_spe_ordered$spe != "Unknown"), ]

matrix_sag_order = rownames(np_reads_sag_reads_colsorted_no_unk)
nature_order = sag_spe_ordered_no_unk[, 1]

common_sags = intersect(matrix_sag_order, nature_order)
leftout_sags = setdiff(nature_order, common_sags)

np_reads_sag_spe_reordered_sub_no_unk = sag_spe_ordered_no_unk[which(sag_spe_ordered_no_unk$sag_id %in% common_sags), ]
np_reads_sag_spe_reordered_sub_no_unk = np_reads_sag_spe_reordered_sub_no_unk[match(matrix_sag_order, np_reads_sag_spe_reordered_sub_no_unk$sag_id), ]

np_row_color_no_unk = sapply(np_reads_sag_spe_reordered_sub_no_unk[, 2], get_all_color_auto_class)

np_row_color_legend_c_no_unk = c(unique(np_row_color_no_unk), "#ffffff")
np_row_color_legend_n_no_unk = unique(np_reads_sag_spe_reordered_sub_no_unk[,2])

sag_spe_cnt = table(sub_sag_spe_ordered[,2])
final_species_with_hits = sag_spe_cnt[unique(sub_sag_spe_ordered[,2])]
legend_labels = c(paste(names(final_species_with_hits), " (", final_species_with_hits, ")", sep = ""), paste0("Total # of SAGs: ", nrow(sag_vp_reads_m)))

np_reads_sag_reads_colsorted_no_unk = np_reads_sag_reads_colsorted_no_unk[, -ncol(np_reads_sag_reads_colsorted_no_unk)]

pdf("fig2_upper_1.pdf", height = 12, width = 8)
par(las=1, cex = 1.2, mar=c(2,5,2,20))
image(t(as.matrix(rev(1:length(np_row_color_legend_c_no_unk)))), col = np_row_color_legend_c_no_unk, xlab = "", xaxt = 'n', yaxt = 'n', bty = 'n')
axis(4, at=seq(0, 1, by = 1/(length(legend_labels)-1)), labels=rev(legend_labels), tick = F)
dev.off()

col = brewer.pal(9, "Set1")[1:3]
heatmap_col = colorRampPalette(c("#F0F0F0", col[2]))(100)

pdf("fig2_upper_2.pdf", height = 10, width = 12)
heatmap.2(as.matrix(np_reads_sag_reads_colsorted_no_unk),Rowv = F,Colv = F,col = heatmap_col,trace = "none", scale ="none",labCol = NULL,density.info = "none",main = "",dendrogram ="none",na.color = "white",key = F,cexRow = 0.25,cexCol= 1, RowSideColors = np_row_color_no_unk, labRow = F)
dev.off()

write.csv(np_reads_sag_spe_reordered_sub_no_unk, file="fig2_row_order_auto_classify.csv", row.names = F, quote = F)

pdf("fig2_upper_3.pdf", height = 4, width = 9.3)
par(cex=1, las = 2)
plot(colSums(np_reads_sag_reads_colsorted_no_unk), pch=20, type="l", col="black", bty = 'n', xlab = "", ylab = "", lwd = 2, xaxt = 'n', ylim = c(0, 80))
#axis(side = 2, at = c(0, 30, 60, 90, 120, 150, 180), labels = c(0, 30, 60, 90, 120, 150, 180))
title(xlab = "Viral Partition Number", line=0)
title(ylab = "# of SAGs", line = 2.2)
points(x = 1:ncol(np_reads_sag_reads_colsorted_no_unk), y=colSums(np_reads_sag_reads_colsorted_no_unk), type = "p", pch = 20)
dev.off()


