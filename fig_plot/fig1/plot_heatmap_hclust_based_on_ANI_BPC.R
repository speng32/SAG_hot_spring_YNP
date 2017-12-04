options(stringsAsFactors = F)
library(gplots)
library(grDevices)
library(devtools)

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
source("../../sag_toolbox.R")

mat_ani = read.csv("../ani_matrix_by_perl_mod1.csv", check.names = F)
mat_bpc = read.csv("../bpc_matrix_by_perl_mod1.csv", check.names = F)
#mat_ani = read.csv("../ani_matrix_by_perl_mod1_stringent.csv", check.names = F)
#mat_bpc = read.csv("../bpc_matrix_by_perl_mod1_stringent.csv", check.names = F)
new_spe_name = read.csv("../Reference_genomes_per_species.csv", header = T)
#classification_file = read.csv("../auto_classification_results_manual_curate_C07_L21.txt", header = F)
classification_file = read.csv("../auto_classification_results.txt", header = F)
name_change_table = read.csv("../species_name_change_auto_classify.csv", header = F)

classification_file_filter12 = classification_file[which(classification_file[, 2] <= 2), ]
classification_file_filter12[, ncol(classification_file_filter12)+1] = substr(classification_file_filter12[, 1], 8, 10)
colnames(classification_file_filter12) = c("sag_full", "class", "spe", "sagid")

# exclude bpc 0's ani
mat_ani[mat_bpc == 0] = 0

mat_ani_bin = mat_ani[, -1]
row.names(mat_ani_bin) = gsub("AD_903_", "", mat_ani[, 1])
mat_ani_bin[is.na(mat_ani_bin)] = 0
mat_ani_bin = mat_ani_bin/100

mat_ani_bin_idx = ifelse(mat_ani_bin >= 0.7, 1, 0)
mat_ani_bin_70_only = mat_ani_bin * mat_ani_bin_idx
mat_ani_bin_filtered = mat_ani_bin_70_only[classification_file_filter12$sagid, ]


col_order = c("Escherichia_coli_str._K-12_substr._MDS42_DNA","Hydrogenobaculum_sp._3684_complete_genome","Metallosphaera_yellowstonensis_MK1","AC-742_N10","Acidianus_hospitalis_W1_complete_genome","Sulfolobus_islandicus_HVE10_4_chromosome_complete_genome","Sulfolobus_solfataricus_P2_complete_genome","AB-777_J04","A_nanophilium","AB-777_K09","AB-777_K20","AB-777_J03","Sulfolobus_tokodaii_str._7_chromosome_complete_genome","Sulfolobus_acidocaldarius_DSM_639_chromosome_complete_genome","AB-777_G06","AB-777_G05","AB-777_L09","Ignicoccus_hospitalis_KIN4_I_complete_genome","AC-742_M05","AC-742_E15","Acidilobus_saccharovorans_345-15","Acidilobus_sp._7A_complete_genome","Acidilobus_sulfurireducans","Thermoproteus_tenax_Kra_1","AB-777_J10","Vulcanisaeta_distributa_DSM_14429_complete_genome","Vulcanisaeta_moutnovskia_768-28_complete_genome","Nanoarchaeum_equitans","Nanoarchaeota_archaeon_7A_complete_genome","N._stetteri","AB-777_F03_Nanoarchaea_sequences","AB-777_O03")



mat_ani_bin_filtered_col_ordered = mat_ani_bin_filtered[, col_order]
sag_class = classification_file_filter12[, c(4, 2)]
#sag_class_color = t(as.matrix(ifelse(sag_class[, 2] == 1, "#33cc33", "#006600")))
sag_class_color = t(as.matrix(ifelse(sag_class[, 2] == 1, "#cccccc", "#000000")))
sag_class_spe = classification_file_filter12[, c(3, 2)]
sag_class_spe[, 1] = sapply(sag_class_spe[, 1], function(x){
  if(x %in% name_change_table[, 1]){
    y = which(as.character(name_change_table[, 1]) == as.character(x))
    return(name_change_table[y, 2])
  }
  else{
    return(x)
  }
})
sag_class_spe_color = sapply(sag_class_spe[, 1], get_all_color_auto_class)

rowsidebar = rbind(sag_class_color, sag_class_spe_color)
rownames(rowsidebar) = c("Class", "Species")

col_order_rename = new_spe_name[match(col_order, new_spe_name[, 1]), 1]
col_order_rename_spe = new_spe_name[match(col_order, new_spe_name[, 1]), 2]

colnames(mat_ani_bin_filtered_col_ordered) = col_order_rename
mat_ani_bin_filtered_col_ordered_spe = mat_ani_bin_filtered_col_ordered
colnames(mat_ani_bin_filtered_col_ordered_spe) = col_order_rename_spe


#myc = colorRampPalette(c(rep("white", 4), "orange", rep("red",1)) )(100)
#myc = c(rep("#FFFFFF", 70), colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(30))
#myc = c(rep("#FFFFFF", 70), colorRampPalette(c("#ffff00", "#ff9900", "#cc6600", "#ff3300", "#ff0000", "#cc0000"))(30))
myc = c(rep("#FFFFFF", 70), colorRampPalette(c("#ffff00", "#ff9900", "#ff6600", "#ff0000", "#cc0000"))(30))


mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

pdf("fig1.pdf", width = 8, height = 10)
#par(cex.main=1)
ani_hm = heatmap.3(as.matrix(mat_ani_bin_filtered_col_ordered), 
          hclustfun = myclust, 
          distfun = mydist, 
          scale = "none", 
          trace = "none",
          na.color = "white",
          col = myc,
          dendrogram = "row", 
          margins = c(6,12),
          Rowv = TRUE, 
          Colv = FALSE, 
          RowSideColors = rowsidebar, 
          symbreaks = FALSE, 
          key = TRUE, 
          symkey = FALSE,
          density.info = "none", 
          main = "", 
          labRow = FALSE, 
          #ColSideColorsSize = 7, 
          RowSideColorsSize = nrow(rowsidebar), 
          KeyValueName="ANI"
          # key.xtickfun = function() {
          #   breaks = seq(0, 1, length = 11)
          #   list(at = parent.frame()$scale01(breaks),
          #        labels = breaks)}
          )

legend("topright",legend=c("Class 1", "Class 2", "Acidocryptum nanophilium", "Acidilobus sp", "Vulcanisaeta sp", "Sulfolobus sp 2", "Sulfolobus sp 1", "Acidianus hospitalis", "Hydrogenobaculum sp", "Nanoarchaea", "Nanoarchaea & Acidocryptum nanophilium", "Nanoarchaea & Sulfolobus sp 1", "Nanoarchaea & Vulcanisaeta sp", "Nanoarchaea & Sulfolobus sp 2", "A. nanophilium & Sulfolobus sp 1"),
       fill=c("#cccccc", "#000000", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#cc00cc"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()


pdf("heatmap_based_on_ani_with_class_spe.pdf", width = 8, height = 10)
#par(cex.main=1)
heatmap.3(as.matrix(mat_ani_bin_filtered_col_ordered_spe), 
          hclustfun = myclust, 
          distfun = mydist, 
          scale = "none", 
          trace = "none",
          na.color = "white",
          col = myc,
          dendrogram = "row", 
          margins = c(6,12),
          Rowv = TRUE, 
          Colv = FALSE, 
          RowSideColors = rowsidebar, 
          symbreaks = FALSE, 
          key = TRUE, 
          symkey = FALSE,
          density.info = "none", 
          main = "", 
          labRow = FALSE, 
          #ColSideColorsSize = 7, 
          RowSideColorsSize = nrow(rowsidebar), 
          KeyValueName="ANI"
          # key.xtickfun = function() {
          #   breaks = seq(0, 1, length = 11)
          #   list(at = parent.frame()$scale01(breaks),
          #        labels = breaks)}
)

legend("topright",legend=c("Class 1", "Class 2", "Acidocryptum nanophilium", "Acidilobus sp", "Vulcanisaeta sp", "Sulfolobus sp 2", "Sulfolobus sp 1", "Acidianus hospitalis", "Hydrogenobaculum sp", "Nanoarchaea", "Nanoarchaea & Acidocryptum nanophilium", "Nanoarchaea & Sulfolobus sp 1", "Nanoarchaea & Vulcanisaeta sp", "Nanoarchaea & Sulfolobus sp 2", "A. nanophilium & Sulfolobus sp 1"),
       fill=c("#cccccc", "#000000", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#cc00cc"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()


#heatmap.2(as.matrix(mat_ani_bin_filtered_col_ordered), trace = "none", scale = "none", density.info = "none", cexRow = 0.25, na.color = "white", col = myc, Colv = F, dendrogram = "row")

# temp = heatmap.2(as.matrix(mat_ani_bin_filtered), trace = "none", scale = "none", density.info = "none", cexRow = 0.25, na.color = "white", col = myc)
# mat_ani_bin_filtered_reordered = mat_ani_bin_filtered[rev(temp$rowInd), temp$colInd]
# 
# write.csv(mat_ani_bin_filtered_reordered, file = "reordered_perl_ani_2way_dendrogram_matrix.csv", quote = F)

#### get ANI order to use for BPC heatmap

ani_hm_order = ani_hm$rowInd

mat_bpc_rowname_change = mat_bpc
mat_bpc_rowname_change[is.na(mat_bpc_rowname_change)] = 0
row.names(mat_bpc_rowname_change) = gsub("AD_903_", "", mat_bpc_rowname_change[, 1])
mat_bpc_rowname_change = mat_bpc_rowname_change[, -1]
mat_bpc_rowname_change_ani_gt70_only = mat_bpc_rowname_change*mat_ani_bin_idx/100
mat_bpc_filtered_ani_gt70_only = mat_bpc_rowname_change_ani_gt70_only[classification_file_filter12$sagid, ]
mat_bpc_filtered_ani_gt70_ani_order = mat_bpc_filtered_ani_gt70_only[rev(ani_hm_order), col_order]

myc_bpc = c(colorRampPalette(c("#ffffff", "#99ccff", "#6699ff", "#3366ff", "#0000ff"))(70), rep("#0000ff", 30))
#myc_bpc = c(colorRampPalette(c("#ffffff", "#99ccff", "#6699ff", "#3366ff", "#0000ff"))(30), rep("#0000ff", 70))


pdf("bpc_heatmap_based_on_ani_order.pdf", width = 8, height = 10)
bpc_hm = heatmap.3(as.matrix(mat_bpc_filtered_ani_gt70_ani_order), 
                   scale = "none", 
                   trace = "none",
                   na.color = "white",
                   col = myc_bpc,
                   dendrogram = "none", 
                   margins = c(6,12),
                   Rowv = FALSE, 
                   Colv = FALSE, 
                   symbreaks = FALSE, 
                   key = TRUE, 
                   symkey = FALSE,
                   density.info = "none", 
                   main = "", 
                   labRow = FALSE, 
                   KeyValueName="BPC")

#legend("topright",legend=c("Class 1", "Class 2", "Acidocryptum nanophilium", "Acidilobus sp", "Vulcanisaeta sp", "Sulfolobus sp 2", "Sulfolobus sp 1", "Acidianus hospitalis", "Hydrogenobaculum sp", "Nanoarchaea", "Nanoarchaea & Acidocryptum nanophilium", "Nanoarchaea & Sulfolobus sp 1", "Nanoarchaea & Vulcanisaeta sp", "Nanoarchaea & Sulfolobus sp 2"), fill=c("#cccccc", "#000000", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()
