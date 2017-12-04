# on cluster
options(stringsAsFactors = F)
options(warn=-1)
suppressMessages(library(GenomicRanges))
suppressMessages(library(Vennerable))
suppressMessages(library(stringr))

args = commandArgs(trailingOnly=TRUE)

checklistfile = args[1]
refspeciesfile = args[2]
outfile = args[3]
outfolder = args[4]

sag_id = str_extract(checklistfile, "[A-Z][0-9][0-9]")

#blast_res = read.table(file.choose(), sep = "\t", header = F)
blast_res = read.table(checklistfile, sep = "\t", header = F)
colnames(blast_res) = c("qid", "sid", "identity", "alignment_length", "mismatches", "gap", "qstart", "qend", "sstart", "send", "evalue", "bit_score", "ref_genome_name")

#ref_genome_species = read.csv("Reference_genomes_per_species.csv", header = T)
ref_genome_species = read.csv(refspeciesfile, header = T)
ref_genome_species = ref_genome_species[, 1:2]
colnames(ref_genome_species) = c("ref_genome_name", "ref_species_name")

blast_res_wspe = merge(blast_res, ref_genome_species, by = "ref_genome_name")

blast_res_wspe_by_refspe = split(blast_res_wspe, f = blast_res_wspe$ref_species_name)

calc_query_interval <- function(x){
  query_start_stop = x[, c("qid", "qstart", "qend")]
  query_start_stop_order = query_start_stop[order(query_start_stop$qid, query_start_stop$qstart, query_start_stop$qend), ]
  
  interval = data.frame(character(), integer(), integer())
  cur_id = query_start_stop_order[1, 1]
  cur_start = query_start_stop_order[1, 2]
  cur_stop = query_start_stop_order[1, 3]
  
  for(i in 2:nrow(query_start_stop_order)){
    iid = query_start_stop_order[i, 1]
    istart = query_start_stop_order[i, 2]
    iend = query_start_stop_order[i, 3]
    if(iid != cur_id){
      interval = rbind(interval, c(cur_id, cur_start, cur_stop))
      cur_id = iid
      cur_start = istart
      cur_stop = iend
    }
    else if(istart > cur_stop){
      interval = rbind(interval, c(cur_id, cur_start, cur_stop))
      cur_start = istart
      cur_stop = iend
    }
    else{
      cur_stop = iend
    }
  }
  interval = rbind(interval, c(cur_id, cur_start, cur_stop))
  colnames(interval) = c("qid","qstart", "qstop")
  return(interval)
}

blast_reduce_interval = lapply(blast_res_wspe_by_refspe, calc_query_interval)

if(length(blast_reduce_interval) == 2){
  r1 = IRanges(start = as.numeric(blast_reduce_interval[[1]]$qstart), end = as.numeric(blast_reduce_interval[[1]]$qstop))
  r2 = IRanges(start = as.numeric(blast_reduce_interval[[2]]$qstart), end = as.numeric(blast_reduce_interval[[2]]$qstop))
  refspe1 = GRanges(seqnames = blast_reduce_interval[[1]]$qid, ranges = r1, strand = NULL)
  refspe2 = GRanges(seqnames = blast_reduce_interval[[2]]$qid, ranges = r2, strand = NULL)
  
  overlapres = findOverlaps(refspe1, refspe2, ignore.strand = T)
  overlapres_t = ranges(overlapres, r1, r2)
  
  spe1_match_length = sum(width(r1))
  spe2_match_length = sum(width(r2))
  overlap_length = sum(width(overlapres_t))
  
  plotfilename = paste0(outfolder, "/", sag_id, "_overlap_venn.pdf")
  pdf(file = plotfilename)
  venn_obj = Venn(SetNames = names(blast_reduce_interval), Weight = c(0, spe1_match_length - overlap_length, spe2_match_length - overlap_length, overlap_length))
  plot(venn_obj)
  grid.text(sag_id, y=0.9, gp=gpar(col="black", cex=2))
  dev.off()
 
  number_info = c(names(blast_reduce_interval)[1], names(blast_reduce_interval)[2], spe1_match_length, spe2_match_length, overlap_length)
}else if(length(blast_reduce_interval) == 3){
  r1 = IRanges(start = as.numeric(blast_reduce_interval[[1]]$qstart), end = as.numeric(blast_reduce_interval[[1]]$qstop))
  r2 = IRanges(start = as.numeric(blast_reduce_interval[[2]]$qstart), end = as.numeric(blast_reduce_interval[[2]]$qstop))
  r3 = IRanges(start = as.numeric(blast_reduce_interval[[3]]$qstart), end = as.numeric(blast_reduce_interval[[3]]$qstop))
  refspe1 = GRanges(seqnames = blast_reduce_interval[[1]]$qid, ranges = r1, strand = NULL)
  refspe2 = GRanges(seqnames = blast_reduce_interval[[2]]$qid, ranges = r2, strand = NULL)
  refspe3 = GRanges(seqnames = blast_reduce_interval[[3]]$qid, ranges = r3, strand = NULL)
  
  overlapres1 = findOverlaps(refspe1, refspe2, ignore.strand = T)
  overlapres_t1 = ranges(overlapres1, r1, r2)
  overlapres2 = findOverlaps(refspe1, refspe3, ignore.strand = T)
  overlapres_t2 = ranges(overlapres2, r1, r3)
  overlapres3 = findOverlaps(refspe2, refspe3, ignore.strand = T)
  overlapres_t3 = ranges(overlapres3, r2, r3)
  
  spe1_match_length = sum(width(r1))
  spe2_match_length = sum(width(r2))
  spe3_match_length = sum(width(r3))
  
  overlap_length1 = sum(width(overlapres_t1))
  overlap_length2 = sum(width(overlapres_t2))
  overlap_length3 = sum(width(overlapres_t3))
  
  plotfilename = paste0(outfolder, "/", sag_id, "_overlap_venn.pdf")
  pdf(file = plotfilename)
  venn_obj1 = Venn(SetNames = names(blast_reduce_interval)[1:2], Weight = c(0, spe1_match_length - overlap_length1, spe2_match_length - overlap_length1, overlap_length1))
  plot(venn_obj1)
  grid.text(sag_id, y=0.9, gp=gpar(col="black", cex=2))
  venn_obj2 = Venn(SetNames = names(blast_reduce_interval)[c(1, 3)], Weight = c(0, spe1_match_length - overlap_length2, spe3_match_length - overlap_length2, overlap_length2))
  plot(venn_obj2)
  grid.text(sag_id, y=0.9, gp=gpar(col="black", cex=2))
  venn_obj3 = Venn(SetNames = names(blast_reduce_interval)[2:3], Weight = c(0, spe2_match_length - overlap_length3, spe3_match_length - overlap_length3, overlap_length3))
  plot(venn_obj3)
  grid.text(sag_id, y=0.9, gp=gpar(col="black", cex=2))
  dev.off()

  number_info = cbind(c(names(blast_reduce_interval)[1], names(blast_reduce_interval)[2], spe1_match_length, spe2_match_length, overlap_length1), c(names(blast_reduce_interval)[1], names(blast_reduce_interval)[3], spe1_match_length, spe3_match_length, overlap_length2), c(names(blast_reduce_interval)[2], names(blast_reduce_interval)[3], spe2_match_length, spe3_match_length, overlap_length3))
}else{
  stop("ERROR! Species should be 2 or 3!")
}

#colnames(number_info) = c("spe1", "spe2", "match_len1", "match_len2", "match_overlap")

write.table(t(number_info), outfile, sep = "\t", row.names = F, col.names = F, quote = F)
