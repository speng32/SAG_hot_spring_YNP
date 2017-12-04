# Rscript to separate reads with both pair end have hits and only a single end have hits

options(stringsAsFactors = F)
suppressMessages(library(dplyr))

args = commandArgs(trailingOnly=T)

FASTAEXTRACTFOLDER=args[1]
OUTPUTFOLDER=args[2]

allfile = list.files(path=FASTAEXTRACTFOLDER)

fasta_file = allfile[grep(".fasta$",allfile)]

split_file <- function(inf){
  if(file.info(paste0(FASTAEXTRACTFOLDER,'/',inf))$size > 0){
    sag_id = gsub(".fasta","",inf)
    file_sname = paste0(sag_id,".single.fasta")
    file_pname = paste0(sag_id,".paired.fasta")
  
    seqfile = read.delim(paste0(FASTAEXTRACTFOLDER,"/",inf), header=F)
    headers = data.frame(head=seqfile[seq(1,nrow(seqfile),2),],rows=seq(1,nrow(seqfile),2))
    headers[,1] = substr(headers[,1], 1, nchar(headers[,1])-2)
    headers_sorted = headers[order(headers[,1]),]
  
    single_line = !(duplicated(headers_sorted[,1]) | duplicated(headers_sorted[,1], fromLast = TRUE))
  
    single_header = headers_sorted[single_line,]
    paired_header = headers_sorted[!single_line,]
  
    if(nrow(single_header)>0){
      single_idx = sort(c(single_header[,2],single_header[,2]+1))
      write.table(seqfile[single_idx,], file=paste0(OUTPUTFOLDER,'/',file_sname),row.names = F, col.names = F, quote = F)
    }
    if(nrow(paired_header)>0){
      paired_idx = sort(c(paired_header[,2],paired_header[,2]+1))
      write.table(seqfile[paired_idx,], file=paste0(OUTPUTFOLDER,'/',file_pname),row.names = F, col.names = F, quote = F)
    }
  }
  else
    print(paste0("File is empty!"))
}

res = lapply(fasta_file, split_file)
