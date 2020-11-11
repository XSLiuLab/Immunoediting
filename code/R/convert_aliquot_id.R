###convert id
nameconvert <- function(x,all){
  path <- paste("~/hrdect/HRDDATA/TCGA/segtxt/",x,"txt/",sep = "")
  setwd(path)
  vcf <- list.files()
  vcfs <- purrr::map(vcf, ~data.table::fread(.))
  names(vcfs) <- sub(".ascat2.allelic_specific.seg.txt", "", vcf)
  vcfs <- data.table::rbindlist(vcfs, use.names = FALSE, idcol = "sample")
  names(vcfs)[2] <- c("aliquot_id")
  
  vcfs <- inner_join(all,vcfs,by="aliquot_id")
  filepath <- paste("~/hrdect/HRDDATA/TCGA/segtxt/rdsdata/",x,".rds",sep = "")
  saveRDS(vcfs, file = filepath)
}
c <- list.files("~/hrdect/HRDDATA/TCGA/segtxt/")
c <- sub("txt", "", c)

for (i in seq_along(c)){
  nameconvert(c[i],all = allquote.rds)
}