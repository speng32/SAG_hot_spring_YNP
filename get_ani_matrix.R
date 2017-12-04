# create ANI result matrix based on the perl code calculation

options(stringsAsFactors = F)
library(reshape2)

# read-in the file
all_ani_list = read.delim("all_ani_results_combined.tsv", header = F)
colnames(all_ani_list) = c("f1", "f2", "ANI")

# change list to matrix
all_ani_m = dcast(all_ani_list, f1 ~ f2, value.var = "ANI")

all_file_name = unique(all_ani_list[, 1])
ref_file = all_file_name[grepl("^ref", all_file_name)]
sag_file = setdiff(all_file_name, sag_file)

row.names(all_ani_m) = all_ani_m[, 1]
all_ani_m = all_ani_m[, -1]

#### missing one sag and check #####
classification = read.csv("new_classification_mod.csv", heade = T)
sag_file_mod = gsub("_", "-", substr(sag_file, 1, 10))
missing = setdiff(classification$SAG, sag_file_mod)

#### create ani matrix #####

ani_m = all_ani_m[sag_file, ref_file]

cname = gsub("ref_", "", colnames(ani_m))
colnames(ani_m) = cname
rownames(ani_m) = gsub("_contigs", "", rownames(ani_m))

write.csv(ani_m, file = "ani_matrix_by_perl.csv", quote = F)




