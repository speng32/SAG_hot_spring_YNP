# Rscript to merge all reads into one file if they hit to the same viral partition regardless of SAG ID

options(stringsAsFactors = F)
suppressMessages(library(dplyr))

args = commandArgs(trailingOnly=T)

FASTAEXTRACTFOLDER=args[1]
MAINVPNFOLDER=args[2]
OUTPUTFOLDER=args[3]

allfile = list.files(path=FASTAEXTRACTFOLDER)

fasta_file = allfile[grep(".fasta$",allfile)]

partition_base <- function(inf){
  sag_id = gsub(".vhits.fasta","",inf)
  f1 = paste0(FASTAEXTRACTFOLDER,'/',inf)
  f2 = paste0(MAINVPNFOLDER,'/',sag_id,".virpar.main")

  if(file.exists(f1) && file.exists(f2) && file.info(f1)$size != 0 && file.info(f2)$size !=0){
    seqs = read.delim(f1, header=F)
    headers = data.frame(head=seqs[seq(1,nrow(seqs),2),],rows=seq(1,nrow(seqs),2))
    #seqs[seq(1,nrow(seqs),2),1] = paste0(seqs[seq(1,nrow(seqs),2),1]," SAG_", sag_id)
  
    main_hits = read.delim(paste0(MAINVPNFOLDER,'/',sag_id,".virpar.main"), header=F)
    main_hits[,1] = paste0("> ",main_hits[,1])
    colnames(main_hits) = c("head","viral_partition")
  
    combined = left_join(x = main_hits, y = headers, by = "head")
  
    uniq_vrp_ids = unique(combined[,2])

    res_in = lapply(uniq_vrp_ids, function(x){
      sub_combined = combined[combined[,2]==x,]
      outfile_name = paste0("viral_partition_",x,".fasta")
    
      lines_to_output = sort(c(sub_combined[,3],sub_combined[,3]+1))
    
      ret = write.table(seqs[lines_to_output,], file=paste0(OUTPUTFOLDER,'/',outfile_name), row.names = F, col.names = F, quote = F, append = T)
    })

    #print(paste0("file ",inf, " finished"))
  }
  else{
    print("file not exist")
  }
}

res_out = lapply(fasta_file,partition_base)
