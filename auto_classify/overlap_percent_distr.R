options(stringsAsFactors = F)

overlap_info = read.csv("all_overlap_numbers.txt", header = F)

overlap_info[, ncol(overlap_info) + 1] = overlap_info[, 3]/overlap_info[, 1]
overlap_info[, ncol(overlap_info) + 1] = overlap_info[, 3]/overlap_info[, 2]

higher_percent = apply(overlap_info, 1, function(x){if(x[4] < x[5]){return(x[5])}else{return(x[4])}})
png("overlap_percent_hist.png")
hist(higher_percent, breaks = 20)
dev.off()

all_percent = c(overlap_info[, 4], overlap_info[, 5])
hist(all_percent, breaks = 40)


total_per = overlap_info[, 3]/(overlap_info[, 1] + overlap_info[, 2] - overlap_info[, 3])