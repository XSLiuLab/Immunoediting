fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

fpkm2tpm <- function(cancertype){
  ##read fpkm data
  f2 = paste("~/neo_dep/RNA-seq/TCGA_pancancer/",cancertype,"/",cancertype,".htseq_fpkm.tsv",sep = "")
  cancer_exp <- data.table::fread(
    file = f2,
    stringsAsFactors = FALSE,
    verbose = FALSE,
    data.table = TRUE,
    showProgress = TRUE,
    header = T,
    quote = ""
  ) 
  
  ##tranform fpkm to tpm
  #the data in the matrix is log2(fpkm+1)
  cancer_exp[,2:ncol(cancer_exp)] <- 2^(cancer_exp[,2:ncol(cancer_exp)])-1
  cancer_exp <- as.data.frame(cancer_exp)
  cancer_exp[,2:ncol(cancer_exp)] <- apply(cancer_exp[,2:ncol(cancer_exp)],2,fpkmToTpm)
  
  ##remove normal samples
  i = which(substr(colnames(cancer_exp),14,15)%in%c("11"))
  cancer_exp = cancer_exp %>% select(c(1,i))
  ## add gene symbol
  map = data.table::fread(
    file = "~/neo_dep/RNA-seq/TCGA_pancancer/acc/gencode.v22.annotation.gene.probeMap",
    stringsAsFactors = FALSE,
    verbose = FALSE,
    data.table = TRUE,
    showProgress = TRUE,
    header = T,
    quote = ""
  )
  cancer_exp <- left_join(cancer_exp,map %>% select(id,gene),by=c("Ensembl_ID"="id")) %>% 
    select(Ensembl_ID,gene,everything())
  if(ncol(cancer_exp)==2){
    return("this cancer doesn't have ‘solid tissue normal’ samples")
  }else{
    saveRDS(cancer_exp,file = paste("~/neo_dep/RNA-seq/TCGA_pancancer/exp_tpm_mat_normal/",cancertype,"_exp.rds",sep = "")) 
  }
}

cancer <- read.table("~/neo_dep/RNA-seq/TCGA_pancancer/cancer")
for (i in cancer$V1){
  fpkm2tpm(i)
}
