library(tidyr)
####apm
apm_gene <- readRDS("data/ap_pathway.rds")
non_slient_mutation <- data.table::fread("~/data/mc3.v0.2.8.PUBLIC.nonsilentGene.xena",data.table = F)
apm_gene <- apm_gene$V1
which(!(apm_gene %in% non_slient_mutation$sample))
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