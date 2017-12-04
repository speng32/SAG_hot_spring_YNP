# Rscript to get main viral partition number from all blast hits

args = commandArgs(trailingOnly=T)

MAINVPNFILE = args[1]
OUTPUTFOLDERALL = args[2]
OUTPUTFOLDERMAIN = args[3]

options(stringsAsFactors = F)
files = list.files(path=OUTPUTFOLDERALL)
virpar_files = files[grep(".virpar.all$",files)]

main_virpar = as.matrix(read.table(MAINVPNFILE, header=F))

only_main_partition <- function(infile){
  inf=paste0(OUTPUTFOLDERALL,"/",infile)
  if(file.info(inf)$size==0){
    print(paste(infile,"is empty, excluded"))
  }
  else{
    outf = gsub(".all",".main",infile)
    content = read.table(inf, header=F, sep="\t")
    content_main = content[content[,2] %in% main_virpar,]
    write.table(content_main, paste0(OUTPUTFOLDERMAIN,"/",outf), row.names=F, col.names=F, quote=F, sep="\t")
    #print(paste(inf,"finished"))
  }
}

res = lapply(virpar_files, only_main_partition)
