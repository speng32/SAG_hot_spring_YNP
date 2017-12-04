# on cluster
options(stringsAsFactors = F)
options(warn=-1)
suppressMessages(library(GenomicRanges))
suppressMessages(library(Vennerable))
suppressMessages(library(stringr))
suppressMessages(library(reshape2))
source("/home/speng32/data/sags/code/auto_classify_toolbotx.R")

args = commandArgs(trailingOnly=TRUE)
all_ani_list_file = args[1]
sag_file_location = args[2]
ref_file_location = args[3]
ref_genome_species_file = args[4]
ref_ref_background_file = args[5]
classify_rules = args[6]	#could be "sf" or "df"
mode = args[7]	#could be "stringent" or "normal"
bpc_cutoff = as.numeric(args[8]) # set bpc cutoff
bpc_single_threshold =as.numeric(args[9])
output_folder = args[10]
tmp_folder_name = paste0(output_folder, ifelse(grepl("/$", output_folder), "temp", "/temp"))

# all_ani_list_file = "/home/speng32/data/sags/ani_mod1/all_ani_output.tsv"
# sag_file_location = "/home/speng32/data/sags/ani_mod1/Corrected_SAG_trimmed_contigs"
# ref_file_location = "/home/speng32/data/sags/ani_mod1/Corrected_SAG_reference_genomes"
# ref_genome_species_file = "/home/speng32/data/sags/Reference_genomes_per_species.csv"
# ref_ref_background_file = "/home/speng32/data/sags/ref_ref_ani/ref_ref_background_overlap.csv"
# classify_rules = "df"
# mode = "stringent"
# bpc_cutoff = 5
# bpc_single_threshold = 30
# output_folder = "./"
# tmp_folder_name = paste0(output_folder, "temp")

# print(all_ani_list_file)
# print(sag_file_location)
# print(ref_file_location)
# print(ref_genome_species_file)
# print(ref_ref_background_file)
# print(mode)
# print(bpc_cutoff)
# print(output_folder)

if(classify_rules != "sf" & classify_rules != "df"){
	stop("Classification Rule error!")
}

if(mode != "stringent" & mode != "normal"){
	stop("Mode error!")
}

if(output_folder != "./"){
	dir.create(output_folder)
}


all_ani_list = read.delim(all_ani_list_file, header = F)
all_ani_list[, ncol(all_ani_list) + 1] = all_ani_list[7]/all_ani_list[6]*100
all_ani_list[, 1] = gsub(".*/|\\.fasta", "", all_ani_list[,1])
all_ani_list[, 2] = gsub(".*/|\\.fasta", "", all_ani_list[,2])
colnames(all_ani_list) = c("f1", "f2", "ANI", "Count_Contig_After_Cut", "Num_Hit", "Total_Length", "Query_Hit_length", "bp_coverage")

if(bpc_cutoff > 0){
	print(paste0("Will cut BPC at ", bpc_cutoff, "%"))
	all_ani_list[all_ani_list[, 8] < bpc_cutoff, 8] = 0
}else{
	print("No BPC cut")
}

print(paste0("Will use ", mode, " mode"))

all_ani_m = dcast(all_ani_list, f1 ~ f2, value.var = "ANI")
all_bpc_m = dcast(all_ani_list, f1 ~ f2, value.var = "bp_coverage")

all_file_name = unique(all_ani_list[, 1])
ref_file = all_file_name[grepl("^ref", all_file_name)]
sag_file = setdiff(all_file_name, ref_file)

sag_file_lookup = cbind(str_extract(sag_file, "[A-P][0-9][0-9]"), paste0(sag_file_location, "/", sag_file, ".fasta"))

row.names(all_ani_m) = all_ani_m[, 1]
all_ani_m = all_ani_m[, -1]
row.names(all_bpc_m) = all_bpc_m[, 1]
all_bpc_m = all_bpc_m[, -1]


#### create ani matrix #####

ani_m = all_ani_m[sag_file, ref_file]
cname = gsub("ref_", "", colnames(ani_m))
colnames(ani_m) = cname
rownames(ani_m) = gsub("_contigs", "", rownames(ani_m))

bpc_m = all_bpc_m[sag_file, ref_file]
cname = gsub("ref_", "", colnames(bpc_m))
colnames(bpc_m) = cname
rownames(bpc_m) = gsub("_contigs", "", rownames(bpc_m))

write.csv(ani_m, file = paste0(output_folder, "/ani_matrix_by_perl_mod1.csv"), quote = F)
write.csv(bpc_m, file = paste0(output_folder, "/bpc_matrix_by_perl_mod1.csv"), quote = F)

# read in mapping ref genome species files
ref_genome_species = read.csv(ref_genome_species_file, header = T)
colnames(ref_genome_species) = c("ref_genome_name", "ref_species_name")

all_ani_l = melt(as.matrix(ani_m), na.rm = F, value.name = "ANI")
all_bpc_l = melt(as.matrix(bpc_m), na.rm = F, value.name = "BPC")

all_ani_bpc_l = cbind(all_ani_l, all_bpc_l[, 3])
colnames(all_ani_bpc_l) = c("SAG", "ref_genome_name", "ANI", "BPC")
all_ani_bpc_l_wspe = merge(ref_genome_species, all_ani_bpc_l, by = "ref_genome_name")
all_ani_bpc_l_wspe_nonull = all_ani_bpc_l_wspe[which(!is.na(all_ani_bpc_l_wspe$BPC)), ]

# read in ref ref background info 

ref_ref_background = read.csv(ref_ref_background_file, header = T)


####### cell classification #######

all_ani_bpc_l_wspe_nonull_lists = split(all_ani_bpc_l_wspe_nonull, f = all_ani_bpc_l_wspe_nonull$SAG)

#sublist = list(all_ani_bpc_l_wspe_nonull_lists[[215]], all_ani_bpc_l_wspe_nonull_lists[[216]], all_ani_bpc_l_wspe_nonull_lists[[217]])
#auto_classify_results = t(sapply(sublist, auto_host_cell_classification))

if(classify_rules == "sf"){
	print("Will perfom single first")
	auto_classify_results = t(sapply(all_ani_bpc_l_wspe_nonull_lists, auto_host_cell_classification_single_first, bpc_single_threshold))
}else{
	print("Will perfom double first")
	auto_classify_results = t(sapply(all_ani_bpc_l_wspe_nonull_lists, auto_host_cell_classification, bpc_single_threshold))
}

write.table(auto_classify_results, paste0(output_folder, "/", "auto_classification_results.txt"), sep = ",", col.names = F, quote = F)

### plot 1
c1_l = all_ani_bpc_l_wspe[which(all_ani_bpc_l_wspe$SAG %in% names(which(auto_classify_results[, 1] == "1"))), ]
c1_ani_m = dcast(c1_l, SAG ~ ref_genome_name, value.var = "ANI")
rownames(c1_ani_m) = c1_ani_m[, 1]
c1_ani_m = c1_ani_m[, -1]
c1_bpc_m = dcast(c1_l, SAG ~ ref_genome_name, value.var = "BPC")
rownames(c1_bpc_m) = c1_bpc_m[, 1]
c1_bpc_m = c1_bpc_m[, -1]

pdf(file = paste0(output_folder, "/ani_vs_bpc_plots_in_species_class1.pdf"))
hstat_species = sapply(1:nrow(c1_ani_m), get_vs_plots_in_species_in_class, c1_ani_m, c1_bpc_m, 1)
dev.off()

### plot 2
c2_l = all_ani_bpc_l_wspe[which(all_ani_bpc_l_wspe$SAG %in% names(which(auto_classify_results[, 1] == "2"))), ]
c2_ani_m = dcast(c2_l, SAG ~ ref_genome_name, value.var = "ANI")
rownames(c2_ani_m) = c2_ani_m[, 1]
c2_ani_m = c2_ani_m[, -1]
c2_bpc_m = dcast(c2_l, SAG ~ ref_genome_name, value.var = "BPC")
rownames(c2_bpc_m) = c2_bpc_m[, 1]
c2_bpc_m = c2_bpc_m[, -1]

pdf(file = paste0(output_folder, "/ani_vs_bpc_plots_in_species_class2.pdf"))
hstat_species = sapply(1:nrow(c2_ani_m), get_vs_plots_in_species_in_class, c2_ani_m, c2_bpc_m, 2)
dev.off()

### plot 2.1
if(mode == "stringent"){
	c21_l = all_ani_bpc_l_wspe[which(all_ani_bpc_l_wspe$SAG %in% names(which(auto_classify_results[, 1] == "2.1"))), ]
	c21_ani_m = dcast(c21_l, SAG ~ ref_genome_name, value.var = "ANI")
	rownames(c21_ani_m) = c21_ani_m[, 1]
	c21_ani_m = c21_ani_m[, -1]
	c21_bpc_m = dcast(c21_l, SAG ~ ref_genome_name, value.var = "BPC")
	rownames(c21_bpc_m) = c21_bpc_m[, 1]
	c21_bpc_m = c21_bpc_m[, -1]

	pdf(file = paste0(output_folder, "/ani_vs_bpc_plots_in_species_class2.1.pdf"))
	hstat_species = sapply(1:nrow(c21_ani_m), get_vs_plots_in_species_in_class, c21_ani_m, c21_bpc_m, 2)
	dev.off()
}

### plot 3
c3_l = all_ani_bpc_l_wspe[which(all_ani_bpc_l_wspe$SAG %in% names(which(auto_classify_results[, 1] == "3"))), ]
c3_ani_m = dcast(c3_l, SAG ~ ref_genome_name, value.var = "ANI")
rownames(c3_ani_m) = c3_ani_m[, 1]
c3_ani_m = c3_ani_m[, -1]
c3_bpc_m = dcast(c3_l, SAG ~ ref_genome_name, value.var = "BPC")
rownames(c3_bpc_m) = c3_bpc_m[, 1]
c3_bpc_m = c3_bpc_m[, -1]

pdf(file = paste0(output_folder, "/ani_vs_bpc_plots_in_species_class3.pdf")
hstat_species = sapply(1:nrow(c3_ani_m), get_vs_plots_in_species_in_class, c3_ani_m, c3_bpc_m, 3)
dev.off()

### plot 4
c4_l = all_ani_bpc_l_wspe[which(all_ani_bpc_l_wspe$SAG %in% names(which(auto_classify_results[, 1] == "4"))), ]
c4_ani_m = dcast(c4_l, SAG ~ ref_genome_name, value.var = "ANI")
rownames(c4_ani_m) = c4_ani_m[, 1]
c4_ani_m = c4_ani_m[, -1]
c4_bpc_m = dcast(c4_l, SAG ~ ref_genome_name, value.var = "BPC")
rownames(c4_bpc_m) = c4_bpc_m[, 1]
c4_bpc_m = c4_bpc_m[, -1]

pdf(file = paste0(output_folder, "/ani_vs_bpc_plots_in_species_class4.pdf")
hstat_species = sapply(1:nrow(c4_ani_m), get_vs_plots_in_species_in_class, c4_ani_m, c4_bpc_m, 4)
dev.off()









