trans_vp_from_original_to_ben <- function(x){
  trans_table = read.csv("Partition_num_change.csv", header=T)
  row_idx = match(x, trans_table$Original.name)
  res = as.character(trans_table$Ben.Name[row_idx])
  return(res)
}

trans_vp_from_ben_to_original <- function(y){
  trans_table = read.csv("Partition_num_change.csv", header=T)
  row_idx = match(y, trans_table$Ben.Name)
  res = as.character(trans_table$Original.name[row_idx])
  return(res)
}

get_color = function(x){
  y = NULL
  if(grepl("Acidianus hospitalis", x, ignore.case = T))
    y = "#ff7f00"
  else if(grepl("Acidilobus sp", x, ignore.case = T))
    y = "#fb9a99"
  else if(grepl("Acidocryptum nanopilium", x, ignore.case = T))
    y = "#a6cee3"
  else if(grepl("Metallospheara sp", x, ignore.case = T))
    y = "#e31a1c"
  else if(grepl("Nanobsidianus stetteri", x, ignore.case = T))
    y = "#fdbf6f"
  else if(grepl("Sulfolobus sp 1", x, ignore.case = T))
    y = "#1f78b4"
  else if(grepl("Sulfolobus Acidocaldarius DSM 639", x, ignore.case = T))
    y = "#cab2d6"
  else if(grepl("Sulfolobus islandicus", x, ignore.case = T))
    y = "#6a3d9a"
  else if(grepl("Sulfolobus sp 2", x, ignore.case = T))
    y = "#33a02c"
  else if(grepl("Sulfolobus_tokodaii_str._7", x, ignore.case = T))
    y = "#ffff99"
  else if(grepl("Vulcanisaeta sp", x, ignore.case = T))
    y = "#b2df8a"
  else
    y = "-1"
  
  return(y)
}

get_color_jan17 <- function(x){
  y = NULL
  if(grepl("A. nanophilium", x, ignore.case = T))
    y = "#a6cee3"
  else if(grepl("Acidiolobus sp", x, ignore.case = T))
    y = "#1f78b4"
  else if(grepl("Vulcanisaeta sp", x, ignore.case = T))
    y = "#b2df8a"
  else if(grepl("Sulfolobus sp 2", x, ignore.case = T))
    y = "#33a02c"
  else if(grepl("Sulfolobus sp 1", x, ignore.case = T))
    y = "#fb9a99"
  else if(grepl("N. stetteri", x, ignore.case = T))
    y = "#e31a1c"
  else if(grepl("A. hospitalis", x, ignore.case = T))
    y = "#fdbf6f"
  else if(grepl("Hydrogenobaculum sp", x, ignore.case = T))
    y = "#ff7f00"
  else
    y = "-1"
  
  return(y)
}

multi_get_color_jan17 <- function(x){
  y = NULL
  if(grepl("A. nanophilium & N. stetteri", x, ignore.case = T))
    y = "#e41a1c"
  else if(grepl("Acidilobus sp & N. stetteri", x, ignore.case = T))
    y = "#377eb8"
  else if(grepl("N. stetteri & Vulcanisaeta", x, ignore.case = T))
    y = "#4daf4a"
  else if(grepl("Sulfolobus sp 1 & N. stetteri", x, ignore.case = T))
    y = "#ff7f00"
  else if(grepl("Sulfolobus sp 2 & N. stetteri", x, ignore.case = T))
    y = "#ffff33"
  else
    y = "-1"
  
  return(y)
}

get_all_color_jan17 <- function(x){
  y = NULL
  if(x == "Acidocryptum nanophilium")
    y = "#a6cee3"
  else if(x == "Acidiolobus sp")
    y = "#1f78b4"
  else if(x == "Vulcanisaeta sp")
    y = "#b2df8a"
  else if(x == "Sulfolobus sp 2")
    y = "#33a02c"
  else if(x == "Sulfolobus sp 1")
    y = "#fb9a99"
  else if(x == "Nanobsidianus stetteri")
    y = "#e31a1c"
  else if(x == "Acidianus hospitalis")
    y = "#fdbf6f"
  else if(x == "Hydrogenobaculum sp")
    y = "#fdbf6f"
  else if(x == "Acidocryptum nanophilium & Nanobsidianus stetteri")
    y = "#ff7f00"
  else if(x == "Acidilobus sp & Nanobsidianus stetteri")
    y = "#cab2d6"
  else if(x == "Vulcanisaeta sp & Nanobsidianus stetteri")
    y = "#6a3d9a"
  else if(x == "Sulfolobus sp 1 & Nanobsidianus stetteri")
    y = "#ffff99"
  else if(x == "Sulfolobus sp 2 & Nanobsidianus stetteri")
    y = "#b15928"
  else
    y = "-1"
  return(y)
}


get_all_color_auto_class <- function(x){
  y = NULL
  if(x == "Acidocryptum nanophilium")
    y = "#a6cee3"
  else if(x == "Acidilobus sp")
    y = "#1f78b4"
  else if(x == "Vulcanisaeta sp")
    y = "#b2df8a"
  else if(x == "Sulfolobus sp 2")
    y = "#33a02c"
  else if(x == "Sulfolobus sp 1")
    y = "#fb9a99"
  else if(x == "Acidianus hospitalis")
    y = "#e31a1c"
  else if(x == "Hydrogenobaculum sp")
    y = "#fdbf6f"
  else if(x == "Nanoarchaea")
    y = "#ff7f00"
  else if(x == "Nanoarchaea & Acidocryptum nanophilium")
    y = "#cab2d6"
  else if(x == "Nanoarchaea & Sulfolobus sp 1")
    y = "#6a3d9a"
  else if(x == "Nanoarchaea & Sulfolobus sp 2")
    y = "#b15928"
  else if(x == "Nanoarchaea & Vulcanisaeta sp")
    y = "#ffff99"
  else if(x == "A. nanophilium & Sulfolobus sp 1")
    y = "#cc00cc"
  else
    y = "-1"
  return(y)
}



# get_color = function(x){
#   y = NULL
#   if(grepl("Acidocryptum nanopilium", x, ignore.case = T))
#     y = "#a6cee3"
#   else if(grepl("Stygiolobus sp", x, ignore.case = T))
#     y = "#1f78b4"
#   else if(grepl("Vulcanisaeta sp", x, ignore.case = T))
#     y = "#b2df8a"
#   else if(grepl("Sulfolobus sp", x, ignore.case = T))
#     y = "#33a02c"
#   else if(grepl("Acidilobus sp", x, ignore.case = T))
#     y = "#fb9a99"
#   else if(grepl("Metallospheara sp", x, ignore.case = T))
#     y = "#e31a1c"
#   else if(grepl("Nanobsidianus stetteri", x, ignore.case = T))
#     y = "#fdbf6f"
#   else if(grepl("Acidianus hospitalis", x, ignore.case = T))
#     y = "#ff7f00"
#   else if(grepl("Sulfolobus Acidocaldarius DSM 639", x, ignore.case = T))
#     y = "#cab2d6"
#   else if(grepl("Sulfolobus islandicus", x, ignore.case = T))
#     y = "#6a3d9a"
#   else if(grepl("Sulfolobus_tokodaii_str._7", x, ignore.case = T))
#     y = "#ffff99"
#   else
#     y = "-1"
#   
#   return(y)
# }

get_matrix_from_edge_list <- function(x){
  x[,2] = as.numeric(x[,2])
  colnames(x) = c("SAG_id","vp","value")
  y = dcast(x, SAG_id ~ vp)
  y[is.na(y)] = 0
  return(y)
}

transform_list_to_matrix_view_crispr <- function(l){
  no_unk_l = cbind(l[complete.cases(l),], Counts = 1)
  m = dcast(data = no_unk_l, SAG_id ~ vp, value.var = "Counts", fun.aggregate = sum)
  return(m)
}