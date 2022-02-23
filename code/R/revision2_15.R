es_exp <- readRDS("data/es_exp_filter_driver.rds")
es_ccf <- readRDS("data/neo_nes_ccf06_1_remove_driver_samples_addp.rds")
escape_all <- readRDS("data/escape_all.rds")
pancancer_subtcells <- readRDS("data/pancancer_subtcells.rds")

es_ccf$cancer <- get_cancer_type(es_ccf$sample)
ccf_immune <- left_join(
  es_ccf,
  pancancer_subtcells
) %>%
  mutate(es_type=ifelse(es<0 & p <0.05,"yes","no")) %>% na.omit(.)
exp_immune <- left_join(
  escape_all,
  pancancer_subtcells %>% mutate(sample=substr(sample,1,12))
) %>% distinct(sample,.keep_all = T) %>% 
  mutate(escape=ifelse(apm_mut==TRUE | over_exp=="yes" | LOH_status=="yes" | es_down=="yes","yes","no")) 
anno <- data.table::fread("data/annotation-tcga.tsv",data.table = F)

ccf_immune_anno <- left_join(
  ccf_immune %>% mutate(sample=substr(sample,1,12)),
  anno %>% rename(sample=V1)
)
MSI <- ccf_immune_anno %>% filter(!is.na(MSI)) %>% 
  filter(nchar(MSI)>=1)

exp_immune <- left_join(
  es_exp,
  pancancer_subtcells
) %>%
  mutate(escape=ifelse(es<0 & p_value <0.05,"yes","no")) %>% na.omit(.)
exp_immune_anno <- left_join(
  exp_immune %>% mutate(sample=substr(sample,1,12)),
  anno %>% rename(sample=V1)
)
exp_MSI <- exp_immune_anno %>% filter(!is.na(MSI)) %>% 
  filter(nchar(MSI)>=1)

p0 <- ggplot(data=MSI,aes(x=es_type,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8 T + NK")
p01 <- ggplot(data=MSI,aes(x=es_type,y=iTreg+nTreg))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="iTreg + nTreg")
p02 <- ggplot(data=MSI,aes(x=es_type,y=log((CD8_T+NK+0.00001)/(iTreg+nTreg+0.00001))))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8_T + NK / iTreg + nTreg")
p02
p1 <- ggplot(data=MSI,aes(x=es_type,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8 T + NK")+
  facet_wrap(~MSI)
p2 <- ggplot(data=MSI,aes(x=es_type,y=iTreg+nTreg))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="iTreg + nTreg")+
  facet_wrap(~MSI)
p3 <- ggplot(data=MSI,aes(x=es_type,y=log((CD8_T+NK+0.00001)/(iTreg+nTreg+0.00001))))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8_T + NK / iTreg + nTreg")+
  facet_wrap(~MSI)
p0 + p01 + p02 + p1 + p2 +p3

p0 <- ggplot(data=ccf_immune,aes(x=es_type,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8 T + NK")+
  facet_wrap(~cancer)
p01 <- ggplot(data=ccf_immune,aes(x=es_type,y=iTreg+nTreg))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="iTreg + nTreg")+
  facet_wrap(~cancer)
p02 <- ggplot(data=ccf_immune,aes(x=es_type,y=log((CD8_T+NK+0.00001)/(iTreg+nTreg+0.00001))))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8_T + NK / iTreg + nTreg")+
  facet_wrap(~cancer)
p02
exp_immune$cancer <- get_cancer_type(exp_immune$sample)
p0 <- ggplot(data=exp_immune,aes(x=escape,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="CD8 T + NK")+
  facet_wrap(~cancer)
p01 <- ggplot(data=exp_immune,aes(x=escape,y=iTreg+nTreg))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="iTreg + nTreg")+
  facet_wrap(~cancer)
p02 <- ggplot(data=exp_immune,aes(x=escape,y=log((CD8_T+NK+0.00001)/(iTreg+nTreg+0.00001))))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="CD8_T + NK / iTreg + nTreg")+
  facet_wrap(~cancer)
p02


p0 <- ggplot(data=ccf_immune,aes(x=es,y=CD8_T+NK))+
  geom_point()+
  theme_prism()+
  labs(x="ES CCF",y="CD8 T + NK")
p01 <- ggplot(data=exp_immune,aes(x=es,y=CD8_T+NK))+
  geom_point()+
  theme_prism()+
  labs(x="ES RNA",y="CD8_T+NK")
p02 <- ggplot(data=ccf_immune,aes(x=es,y=iTreg + nTreg))+
  geom_point()+
  theme_prism()+
  labs(x="ES CCF",y="iTreg + nTreg")
p03 <- ggplot(data=exp_immune,aes(x=es,y=iTreg + nTreg))+
  geom_point()+
  theme_prism()+
  labs(x="ES RNA",y="iTreg + nTreg")
p0 + p01 + p02 +p03

p4 <- ggplot(data=exp_MSI,aes(x=es,y=iTreg+nTreg))+
  geom_point()+
  geom_smooth(method = lm)+
  stat_cor()+
  theme_prism()+
  labs(x="ES EXP",y="iTreg+nTreg")
p5 <- ggplot(data=exp_MSI,aes(x=es,y=CD8_T+NK))+
  geom_point()+
  geom_smooth(method = lm)+
  stat_cor()+
  theme_prism()+
  labs(x="ES EXP",y="CD8_T+NK")
p6 <- ggplot(data=exp_MSI,aes(x=es,y=iTreg+nTreg))+
  geom_point()+
  geom_smooth(method = lm)+
  stat_cor()+
  theme_prism()+
  labs(x="ES EXP",y="iTreg+nTreg")+
  facet_wrap(~MSI)
p7 <- ggplot(data=exp_MSI,aes(x=es,y=CD8_T+NK))+
  geom_point()+
  geom_smooth(method = lm)+
  stat_cor()+
  theme_prism()+
  labs(x="ES EXP",y="CD8_T+NK")+
  facet_wrap(~MSI)
(p4 + p5) / (p6 + p7)

ggplot(data = es_ccf,aes(x=p,y=log(all_counts)))+
  geom_point()+
  theme_prism()

##########免疫治疗
#clonal TMB
all_mut_ccf_ici <- readRDS("data/Immunotherapy/all_mut_ccf_ici.rds")
all_clinical <- readRDS("data/Immunotherapy/all_clinical.rds")

neo_nes <- readRDS("data/Immunotherapy/nes_immunetherapy.rds")
neo_nes <- left_join(neo_nes,all_clinical,by="sample")
neo_nes <- neo_nes %>% filter(!is.na(response2))
neo_nes <- neo_nes %>% 
  filter(!(sample %in% c("nadeem_Pt85","liu_Patient132")))
neo_missense %>% group_by(sample) %>% 
  summarise(neo_counts=sum(neo=="yes"),all_counts=n()) -> summ

all_mut_ccf_ici %>%
  group_by(sample) %>%
  summarise(all_tmb=n()/38,clonal_tmb=sum(cancer_cell_frac>=0.9)/38) -> tmb

neo_nes <- left_join(neo_nes,tmb)
library(ezcox)
show_forest(neo_nes,covariates = "es",time = "OS.time",status = "OS")

neo_nes <- neo_nes %>%
  mutate(obj=ifelse(response2=="response",1,0))
model2 <- glm(obj ~ es, data = neo_nes, family = "binomial")
summary(model2)
performance::performance_hosmer(model2)
ggplot(data=neo_nes,aes(x=response2,y=es))+
  geom_boxplot()+
  stat_compare_means()

####Indel mutation
liu_indel <- data.table::fread("data/Immunotherapy/liu_indel_counts",data.table = F)
liu_indel <- liu_indel %>% 
  mutate(sample=gsub("_pass.vcf","",V1)) %>% 
  mutate(sample=paste0("liu_",sample)) %>% 
  mutate(indel_counts=V2+V3) %>% 
  select(sample,indel_counts)
nadeem_indel <- data.table::fread("data/Immunotherapy/nadeem_indel_counts",data.table = F)
nadeem_indel <- nadeem_indel %>% 
  mutate(sample=gsub("_pass.vcf","",V1)) %>% 
  mutate(sample=paste0("nadeem_",sample)) %>% 
  mutate(indel_counts=V2+V3) %>% 
  select(sample,indel_counts)
willy_indel <- data.table::fread("data/Immunotherapy/willy_hugo_indel_counts",data.table = F)
willy_indel <- willy_indel %>% 
  mutate(sample=gsub("_pass.vcf","",V1)) %>% 
  mutate(sample=paste0("willy_",sample)) %>% 
  mutate(indel_counts=V2+V3) %>% 
  select(sample,indel_counts)

all_indel <- bind_rows(liu_indel,nadeem_indel,willy_indel)
saveRDS(all_indel,file = "data/Immunotherapy/all_indel_counts.rds")

neo_nes <- left_join(neo_nes,all_indel)
show_forest(neo_nes,covariates = "indel_counts",time = "OS.time",status = "OS")
model2 <- glm(obj ~ es, data = neo_nes, family = "binomial")
summary(model2)
performance::performance_hosmer(model2)
###注释一下免疫治疗的数据，看看有没有新抗原突变和driver在同一个基因上的样本
##ref alt
library(readr)
files <- list.files("../../evolution/immunetherapy/")
files <- files[grepl("neoantigens.unfiltered.txt",files)]
re <- vector("list",25)
for (i in 1:25){
  mut <- read_table2(paste("../../evolution/immunetherapy/",files[i],sep = ""),col_names = paste0(rep("V",27),c(1:27)))
  
  mut <- mut %>%
    select(c(1,3:10,20:27))
  colnames(mut) <- c("sample","chr","position",
                     "ref","alt","gene","exp","res_pos",
                     "hla","score_el","rank_el","score_ba",
                     "rank_ba","IC50","candidate","bindlevel",
                     "novelty")
  mut <- mut %>%
    mutate(novelty=ifelse(is.na(novelty),0,novelty)) %>%
    mutate(index=paste(sample,chr,position,ref,alt,sep = ":"))
  
  neo <- mut %>%
    filter(IC50<50 & bindlevel=="SB" & novelty==1 & exp>1) %>%
    distinct(index,.keep_all = T)
  
  mut <- mut %>%
    distinct(index,.keep_all = T) %>%
    mutate(neo = ifelse(index %in% neo$index , "neo","not_neo"))
  mut <- mut %>% 
    select(sample,chr,position,gene,exp,neo,ref,alt) %>% 
    mutate(gene=gsub("\\:.+","",gene))
  mut$sample <- paste(gsub("_.+","",files[i]),mut$sample,sep = "_")
  re[[i]] <- mut
}

all_mut <- bind_rows(re)
saveRDS(all_mut,file = "../tmp/all_mut_ici_withrefalt.rds")

all_mut_ici <- readRDS("../tmp/all_mut_ici_withrefalt.rds")
all_ccf <- readRDS("data/Immunotherapy/all_ccf.rds")
all_mut_ccf <- left_join(
  all_mut_ici %>% 
    mutate(index=paste(sample,chr,position,sep = ":")),
  all_ccf %>% 
    mutate(Chromosome=paste0("chr",Chromosome)) %>% 
    mutate(index=paste(sample,Chromosome,Start_position,sep = ":")) %>% 
    select(index,cancer_cell_frac,purity),
  by="index"
)
all_mut_ccf <- all_mut_ccf %>% filter(!is.na(cancer_cell_frac))
saveRDS(all_mut_ccf,file = "../tmp/all_mut_ccf_ici_withrefalt.rds")

all_mut_ccf_ici <- readRDS("../tmp/all_mut_ccf_ici_withrefalt.rds")
anno_input <- all_mut_ccf_ici %>%
  mutate(end=position) %>% 
  select(chr,position,end,ref,alt,index)
write.table(anno_input,file = "../tmp/anno_input_ici.txt",sep = "\t",row.names = F,col.names = F,quote = F)

anno_out <- data.table::fread("../tmp/ici_anno.exonic_variant_function",data.table = F)
anno_out <- anno_out %>% 
  rowwise() %>% 
  mutate(gene= strsplit(V3,split = ":")[[1]][1])
anno_out$is_driver <- NA
for (i in 1:nrow(anno_out)){
  dt <- driver_mutations %>% 
    filter(gene == anno_out$gene[i]) %>% 
    filter(CODE == "SKCM")
  a <- sapply(dt$protein_change,
              function(x,y){grep(x,y)},
              y=strsplit(x = anno_out$V3[i],split = ":")[[1]])
  anno_out$is_driver[i] <- ifelse(any(lengths(a)>=1),"yes","no")
}
anno_out <- anno_out %>% select(V9,is_driver) %>% rename(index=V9)
all_mut_ccf_ici <- left_join(
  all_mut_ccf_ici,anno_out
)
sample_remove <- all_mut_ccf_ici %>% 
  group_by(sample) %>%
  summarise(inter_gene=intersect(gene[neo=="neo"],
                                 gene[is_driver=="yes"]))
##两个样本需要去掉
# A tibble: 2 x 2
# Groups:   sample [2]
# sample         inter_gene
# <chr>          <chr>     
# 1 liu_Patient132 MAP2K1    
# 2 nadeem_Pt85    FAT1 

##signature
all_mut_ccf_ici <- readRDS("../tmp/all_mut_ccf_ici_withrefalt.rds")
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)
all_mut_ccf_ici <- as.data.frame(all_mut_ccf_ici)
sigs.input <- mut.to.sigs.input(mut.ref = all_mut_ccf_ici, 
                                sample.id = "sample", 
                                chr = "chr", 
                                pos = "position", 
                                ref = "ref", 
                                alt = "alt",bsg = BSgenome.Hsapiens.UCSC.hg38)

mm <- matrix(data=rep(1,30*207),nrow = 207,ncol = 30)
colnames(mm) <- rownames(signatures.cosmic)
sig_ici <- as.data.frame(mm) %>% mutate(sample=rownames(sigs.input))
for (i in 1:nrow(sig_ici)){
  sample_1 <-  whichSignatures(tumor.ref = sigs.input, 
                               signatures.ref = signatures.cosmic, 
                               sample.id = sig_ici$sample[i], 
                               contexts.needed = TRUE,
                               tri.counts.method = 'default')
  sig_ici[i,1:30] <- sample_1$weights
}
saveRDS(sig_ici,file = "data/sig_ici.rds")
##筛选需要的sig：4，7，2，13
sig_ici <- sig_ici %>% 
  select(sample,Signature.2,Signature.4,Signature.13,Signature.7)
saveRDS(sig_ici,file = "data/sig_ici_need.rds")

###
escape_all <- escape_all %>% 
  mutate(escape=ifelse(apm_mut==TRUE | over_exp=="yes" | LOH_status=="yes" | es_down=="yes","yes","no")) 
sim_all <- readRDS("../tmp/sim_2000_not_filter_all_ccf.rds")
ccf_escape <- left_join(
  es_ccf %>% mutate(sample=substr(sample,1,12)),
  escape_all %>% select(sample,escape)
) %>% filter(!is.na(escape))
escape <- ccf_escape %>% filter(escape=="yes")
no_escape <- ccf_escape %>% filter(escape=="no")
sim_escape <- sim_all %>% 
  filter(sample %in% escape$sample)
sim_no_escape <- sim_all %>% 
  filter(sample %in% no_escape$sample)
sim_escape %>%
  group_by(sim_num) %>%
  summarise(median_es=median(es)) -> summ1
sim_no_escape %>%
  group_by(sim_num) %>%
  summarise(median_es=median(es)) -> summ2
p1 <- WVPlots::ShadedDensity(frame = summ1, 
                            xvar = "median_es",
                            threshold = median(escape$es),
                            title = "",
                            tail = "left")

p1$layers[[1]]$aes_params$colour <- "red"
p1$layers[[1]]$aes_params$size <- 1
p1$layers[[2]]$aes_params$fill <- "blue"     #geom_ribbon
p1$layers[[3]]$aes_params$colour <- "black" 
p1$layers[[3]]$aes_params$size <- 1
p1 <- p1 + labs(x="Simulation median es")+
  theme_prism()
p1

p2 <- WVPlots::ShadedDensity(frame = summ2, 
                             xvar = "median_es",
                             threshold = median(no_escape$es),
                             title = "",
                             tail = "left")

p2$layers[[1]]$aes_params$colour <- "red"
p2$layers[[1]]$aes_params$size <- 1
p2$layers[[2]]$aes_params$fill <- "blue"     #geom_ribbon
p2$layers[[3]]$aes_params$colour <- "black" 
p2$layers[[3]]$aes_params$size <- 1
p2 <- p2 + labs(x="Simulation median es")+
  theme_prism()
p2

ggplot(data=es_ccf,aes(x=escape,y=es))+
  geom_boxplot()+
  stat_compare_means()
##
checkpoint_over <- readRDS("data/checkpoint_over.rds") %>% 
  mutate(type=ifelse(PDL1_over=="yes" | CTLA4_over=="yes","yes","no"))
checkpoint_over_escape <- checkpoint_over %>% filter(type=="yes")
apm_mut_sample <- readRDS("data/apm_mut_sample.rds")
apm_mut_escape <- apm_mut_sample %>% filter(apm_mut)
LOH <- readRDS("data/HLA_loh.rds")
LOH_escape <- LOH %>% filter(LOH_status=="yes")
nes_exp <- readRDS("data/es_exp_filter_driver.rds")
es_escape <- nes_exp %>% filter(es<0 & p_value<0.05)

escape_samples <- Reduce(union,list(substr(checkpoint_over_escape$sample,1,12),
                            substr(apm_mut_escape$sample,1,12),
                            substr(LOH_escape$Sample_Barcode,1,12),
                            substr(es_escape$sample,1,12)))
es_ccf <- es_ccf %>% 
  mutate(escape=ifelse(substr(sample,1,12) %in% escape_samples,"yes","no"))
table(es_ccf$escape)
escape <- es_ccf %>% filter(escape=="yes")
no_escape <- es_ccf %>% filter(escape=="no")
