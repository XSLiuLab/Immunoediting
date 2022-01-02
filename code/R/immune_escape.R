library(tidyr)
abs_seg <- data.table::fread("~/data/TCGA_mastercalls.abs_segtabs.fixed.txt",sep = "\t",header = T,data.table = F)
chr6_loh <- abs_seg %>% filter(Chromosome==6) %>% 
  filter(LOH==1)
###get HLA gene
anno <- data.table::fread("~/data/gencode.v32lift37.annotation.gtf.gz",data.table = F)
anno <- anno %>% 
  separate(col = V9,into = c("id","type","name"),sep = ";") %>% 
  filter(V3=="gene") %>% 
  mutate(name=gsub("gene_name","",name)) %>% 
  mutate(name=gsub('"',"",name)) %>% 
  mutate(name=gsub(" ","",name))
hla <- anno %>% 
  filter(name %in% c("HLA-A","HLA-B","HLA-C"))

check_loh <- function(start,end,hla_start,hla_end){
  if((hla_end>start & hla_end<end) | (hla_start>start & hla_start<end) | (hla_start<start & hla_end>end)){
    return("yes")
  }else{return("no")}
}
chr6_loh$hla_a_in <- mapply(check_loh,chr6_loh$Start,chr6_loh$End,29909037,29913661)
chr6_loh$hla_b_in <- mapply(check_loh,chr6_loh$Start,chr6_loh$End,31237268,31324965)
chr6_loh$hla_c_in <- mapply(check_loh,chr6_loh$Start,chr6_loh$End,31236526,31239907)

chr6_loh %>% group_by(Sample) %>% 
  summarise(HLA_A=any(hla_a_in=="yes"),
            HLA_B=any(hla_b_in=="yes"),
            HLA_C=any(hla_c_in=="yes")) -> hla_loh
saveRDS(hla_loh,file = "data/HLA_LOH.rds")
####apm
apm_gene <- readRDS("data/ap_pathway.rds")
# apm_gene <- apm_gene[c(grep("HLA-A",apm_gene),grep("B2M",apm_gene),grep("HLA-B",apm_gene),
#                        grep("HLA-C",apm_gene))]
non_slient_mutation <- data.table::fread("~/data/mc3.v0.2.8.PUBLIC.nonsilentGene.xena",data.table = F)

which(!(apm_gene %in% non_slient_mutation$sample))
apm_gene[31]
##基因CD8B2改名了，之前是叫CD8BP
"CD8BP" %in% non_slient_mutation$sample

apm_gene[31] <- "CD8BP"
non_slient_mutation <- non_slient_mutation %>%
  dplyr::filter(sample %in% apm_gene) 
check <- apply(non_slient_mutation[,2:ncol(non_slient_mutation)],2,function(x){any(x==1)}) %>% as.data.frame()
colnames(check) <- "apm_mut"
check$sample <- rownames(check)
saveRDS(check,file = "data/apm_mut_sample.rds")
###checkpoint over expression CD274 CTLA-4
tpm_gene <- readRDS("~/data/tpm_trans.rds")
CD274 <- tpm_gene %>% filter(gene == "CD274") %>%
  select(c(3:ncol(tpm_gene))) %>% t() %>% as.data.frame()
CD274$sample <- rownames(CD274)
colnames(CD274)[1] <- "PDL1"
CTLA4 <- tpm_gene %>% filter(gene == "CTLA4") %>%
  select(c(3:ncol(tpm_gene))) %>% t() %>% as.data.frame()
CTLA4$sample <- rownames(CTLA4)
colnames(CTLA4)[1] <- "CTLA4"

check_point <- left_join(CD274,CTLA4,by="sample")
saveRDS(check_point,file = "data/checkpoint_exp.rds")
table(substr(check_point$sample,14,15))

normal <- data.frame(sample=colnames(tpm_gene)[which(substr(colnames(tpm_gene),14,15)==11)])
normal <- left_join(normal,check_point,by="sample")
normal$cancer <- EasyBioinfo::get_cancer_type(normal$sample)
normal %>%
  group_by(cancer) %>%
  summarise(mean_pdl1 = mean(PDL1),
            sd_pdl1 = sd(PDL1),
            mean_ctla4 = mean(CTLA4),
            sd_ctla4 = sd(CTLA4)) -> mean_sd

cancer <- check_point %>% 
  filter(as.numeric(substr(sample,14,15)) < 11)
cancer$cancer <- EasyBioinfo::get_cancer_type(cancer$sample)
cancer <- left_join(cancer,mean_sd,by="cancer")
cancer <- na.omit(cancer)

cancer <- cancer %>%
  mutate(PDL1_over = ifelse(PDL1 > (mean_pdl1+2*sd_pdl1),"yes","no"),
         CTLA4_over = ifelse(CTLA4 > (mean_ctla4+2*sd_ctla4),"yes","no")) 
saveRDS(cancer,file = "data/checkpoint_over.rds")
# HLA_LOH <- readRDS("HLA_LOH.rds")
# cancer <- left_join(
#   cancer,
#   HLA_LOH %>% rename(sample=Sample),
#   by="sample"
# )
# 
# cancer$HLA_LOH <- mapply(function(x,y,z){
#   if(all(is.na(c(x,y,z)))){
#     return(NA)
#   }else if (any(c(x,y,z))){
#     return("yes")
#   }else{
#     return("no")
#   }
# },cancer$HLA_A,cancer$HLA_B,cancer$HLA_C)
# 
# saveRDS(cancer,file = "escape_3index.rds")