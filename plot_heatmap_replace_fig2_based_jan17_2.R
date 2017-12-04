# plot cell level viral-host interaction heatmap

options(stringsAsFactors = F)
library(reshape2)
library(gplots)
library(RColorBrewer)
library(data.table)

source("sag_toolbox.R")

#read-in host cell species classification and format species names
sag_spe_raw = read.csv("SAG_species_classification_Jan17.csv", header = T)
sag_spe_raw[, 2] = gsub(" and ", " & ", sag_spe_raw[, 2])
sag_spe_raw = sag_spe_raw[grep("ND", sag_spe_raw[, 2], invert = T), ]

sag_name_changer = read.csv("species_Name_change_jan17.csv", header = F)
sag_spe_raw[, 2] = sapply(sag_spe_raw[, 2], function(x){
  if(x %in% sag_name_changer[, 1]){
    y = which(as.character(sag_name_changer[, 1]) == as.character(x))
    return(sag_name_changer[y, 2])
  }
  else{
    return(x)
  }
})

sag_spe = data.frame(sag_id = substr(sag_spe_raw$SAG, 8, 10), spe = sag_spe_raw$Species.classification)
roworder = as.character(read.table("heatmap_row_order_jan17_2.txt", header=F, sep=",")[, 1])
spe_ordered = ordered(sag_spe[,2], levels = roworder)
sag_spe_ordered = sag_spe[order(spe_ordered), ]

sag_vp_inf_m = read.table("matrix_infection.csv", header = T, sep = ",")
sag_vp_crispr_m = read.table("matrix_crispr.csv", header = T, sep = ",")

#######
# sort columns

# sortCols = function(data){
#   data[, order(colSums(data), decreasing = T)]
# }

# sortRows <- function(data){
#   order_string = paste("data[order(", paste(paste0("data[,", 1:ncol(data), "]"), collapse=","), ", decreasing = T), ]", collapse = "", sep = "")
#   data = eval(parse(text = order_string))
# }

# sortAllRows <- function(data_m, row_spe){
#   ordered_m = NULL
#   for(i in unique(row_spe[, 2])){
#     #print(i)
#     i_idx = which(row_spe[, 2] == i)
#     if(length(i_idx) > 1){
#       ordered_m = rbind(ordered_m, sortRows(data_m[i_idx, ]))
#     }else{
#       rname = rownames(ordered_m)
#       ordered_m = rbind(ordered_m, data_m[i_idx, ])
#       rownames(ordered_m) = c(rname, rownames(data_m)[i_idx])
#     }
#   }
#   return(ordered_m)
# }


# sortAllRows_order_by_heatmap_cnt <- function(data_m, row_spe){
#   ordered_m = NULL
#   for(i in names(sort(table(row_spe[,2]), decreasing=T))){
#     #print(i)
#     i_idx = which(row_spe[, 2] == i)
#     if(length(i_idx) > 1){
#       ordered_m = rbind(ordered_m, sortRows(data_m[i_idx, ]))
#     }else{
#       rname = rownames(ordered_m)
#       ordered_m = rbind(ordered_m, data_m[i_idx, ])
#       rownames(ordered_m) = c(rname, rownames(data_m)[i_idx])
#     }
#   }
#   return(ordered_m)
# }


####################################### reads only #######################################

# create heatmap for viral signal only
sag_vp_reads_m = sag_vp_inf_m
rownames(sag_vp_reads_m) = substr(sag_vp_reads_m[, 1], 8, 10)
sag_vp_reads_m = sag_vp_reads_m[, -1]
sag_vp_reads_m = ifelse(sag_vp_reads_m>0, 1, 0)
sag_vp_reads_m[is.na(sag_vp_reads_m)] = 0

sag_spe_ordered_sub_reads = sag_spe_ordered[which(sag_spe_ordered$sag_id %in% row.names(sag_vp_reads_m)),]
sag_vp_reads_m = sag_vp_reads_m[sag_spe_ordered_sub_reads$sag_id, ]

# save heatmap row order to file
#write.csv(sag_spe_ordered_sub_reads, file="heatmap_row_order_read_only_jan17.csv", row.names = F, quote = F)


sag_spe_ordered_new_order = sag_spe_ordered

sag_spe_ordered_sub_reads_new_order = sag_spe_ordered_new_order[which(sag_spe_ordered_new_order$sag_id %in% row.names(sag_vp_reads_m)), ]
sag_vp_reads_m_new_order = sag_vp_reads_m
sag_vp_reads_m_new_order = sag_vp_reads_m_new_order[sag_spe_ordered_sub_reads_new_order$sag_id, ]

# sort columns for heatmap
np_reads_sag_reads_colsorted = sortCols(sag_vp_reads_m_new_order)
colnames(np_reads_sag_reads_colsorted) = gsub("X", "", colnames(np_reads_sag_reads_colsorted))
sub_sag_spe_ordered = sag_spe_ordered[which(sag_spe_ordered[, 1] %in% rownames(np_reads_sag_reads_colsorted)), ]
np_reads_sag_reads_colsorted = sortAllRows(np_reads_sag_reads_colsorted, sub_sag_spe_ordered)

# remove undecided samples for the heatmap
np_reads_sag_reads_colsorted_no_unk = np_reads_sag_reads_colsorted
sag_spe_ordered_no_unk = sag_spe_ordered[which(sag_spe_ordered$spe != "Unknown"), ]


np_reads_sag_spe_reordered_sub_no_unk = sag_spe_ordered_no_unk[which(sag_spe_ordered_no_unk$sag_id %in% row.names(np_reads_sag_reads_colsorted_no_unk)),]
np_reads_sag_spe_reordered_sub_no_unk = np_reads_sag_spe_reordered_sub_no_unk[match(row.names(np_reads_sag_reads_colsorted_no_unk), np_reads_sag_spe_reordered_sub_no_unk$sag_id), ]

# get host species color for the heatmap
np_row_color_no_unk = sapply(np_reads_sag_spe_reordered_sub_no_unk[, 2], get_all_color_jan17)

# get legends for heatmap host species
np_row_color_legend_c_no_unk = unique(np_row_color_no_unk)
np_row_color_legend_n_no_unk = unique(np_reads_sag_spe_reordered_sub_no_unk[,2])

sag_spe_cnt = table(sub_sag_spe_ordered[,2])
final_species_with_hits = sag_spe_cnt[unique(sub_sag_spe_ordered[,2])]
legend_labels = paste(names(final_species_with_hits), " (", final_species_with_hits, ")", sep = "")

np_reads_sag_reads_colsorted_no_unk = np_reads_sag_reads_colsorted_no_unk[, -ncol(np_reads_sag_reads_colsorted_no_unk)]

# plot host species color key for the heatmap
pdf("./fig2/np_reads_only_color_key_no_unk_jan17.pdf", height = 12, width = 8)
par(las=1, cex = 1.2, mar=c(2,5,2,20))
image(t(as.matrix(rev(1:length(np_row_color_legend_c_no_unk)))), col = np_row_color_legend_c_no_unk, xlab = "", xaxt = 'n', yaxt = 'n', bty = 'n')
axis(4, at=seq(0, 1, by = 1/11), labels=rev(legend_labels), tick = F)
dev.off()

# plot interaction network heatmap
col = brewer.pal(9, "Set1")[1:3]
heatmap_col = colorRampPalette(c("#F0F0F0", col[2]))(100)
pdf("./fig2/np_reads_only_vs_vp_no_unk_jan17.pdf", height = 10, width = 12)
heatmap.2(as.matrix(np_reads_sag_reads_colsorted_no_unk),Rowv = F,Colv = F,col = heatmap_col,trace = "none", scale ="none",labCol = NULL,density.info = "none",main = "",dendrogram ="none",na.color = "white",key = F,cexRow = 0.25,cexCol= 1, RowSideColors = np_row_color_no_unk)
dev.off()

# save the order of rows for the heatmap
write.csv(np_reads_sag_spe_reordered_sub_no_unk, file="./fig2/np_heatmap_row_order_read_only_no_unk_jan17.csv", row.names = F, quote = F)

# plot line plot for number of host cells the viral signal found in a decreasing order
pdf("./fig2/np_hist_reads_only_no_unk_jan17.pdf", height = 4, width = 9.4)
par(cex=1)
plot(colSums(np_reads_sag_reads_colsorted_no_unk), pch=20, type="l", col="black", bty = 'n', xlab = "", ylab = "", lwd = 2, xaxt = 'n', ylim = c(0, 200))
#axis(side = 2, at = c(0, 30, 60, 90, 120, 150, 180), labels = c(0, 30, 60, 90, 120, 150, 180))
title(xlab = "Viral Partition Number", line=0)
title(ylab = "# of SAGs", line = 2.2)
points(x = 1:ncol(np_reads_sag_reads_colsorted_no_unk), y=colSums(np_reads_sag_reads_colsorted_no_unk), type = "p", pch = 20)
dev.off()


