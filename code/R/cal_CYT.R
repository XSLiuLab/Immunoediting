##cytolytic activity
##Cytolytic activity (CYT) was calculated as the geometric mean of GZMA and PRF1 (as expressed in TPM, 0.01 offset)
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}
sample <- data.frame(sample=NA,cyt=NA)

file <- dir("~/neo_dep/RNA-seq/TCGA_pancancer/exp_matrix_tpm/")
for(i in 1:33){
  load(paste("~/neo_dep/RNA-seq/TCGA_pancancer/exp_matrix_tpm/",file[i],sep = ""))
  s <- data.frame(sample=colnames(cancer_exp)[-(1:2)])
  cancer_exp <- cancer_exp %>% filter(gene=="GZMA"|gene=="PRF1")
  s$cyt <- apply(cancer_exp[,3:ncol(cancer_exp)],2,gm_mean)
  sample <- rbind(sample,s)
  rm(cancer_exp)
}
sample <- na.omit(sample)
saveRDS(sample,file = "~/neo_dep/ssNEA/data/pancancer_cyt.rds")
