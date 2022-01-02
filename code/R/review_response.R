library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggprism)
library(patchwork)
library(parallel)
##review1

# 使用cd8/treg --------------------------------------------------------------


pancancer_subtcells <- readRDS("data/pancancer_subtcells.rds")
nes_ccf <- readRDS("data/neo_nes_ccf06_1_remove_driver_samples_addp.rds")
nes_exp <- readRDS("data/nes_exp_addp.rds")

ccf_immune <- left_join(
  nes_ccf,
  pancancer_subtcells
) %>%
  mutate(es_type=ifelse(es<0 & p <0.05,"yes","no")) %>% 
  mutate(`CD8/Treg`=exp(CD8_T)/exp(iTreg+nTreg)) %>% 
  mutate(`NK/Treg`=exp(NK)/exp(iTreg+nTreg))

exp_immune <- left_join(
  nes_exp ,
  pancancer_subtcells
) %>%
  mutate(es_type=ifelse(es<0 & p_value <0.05,"yes","no")) %>% 
  mutate(`CD8/Treg`=exp(CD8_T)/exp(iTreg+nTreg)) %>% 
  mutate(`NK/Treg`=exp(NK)/exp(iTreg+nTreg))

p1 <- ggplot(data=ccf_immune,aes(x=es_type,y=`CD8/Treg`))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8/Treg")
p2 <- ggplot(data=ccf_immune,aes(x=es_type,y=`NK/Treg`))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="NK/Treg")
p1 + p2

p3 <- ggplot(data=exp_immune,aes(x=es_type,y=`CD8/Treg`))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="CD8/Treg")
p4 <- ggplot(data=exp_immune,aes(x=es_type,y=`NK/Treg`))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="NK/Treg")
p3 + p4



####Reviewer2 ---------------
# 免疫逃逸的其他机制 ---------------------------------------------------------------
##重新算一下逃逸的样本--immune_escape.R

checkpoint_exp <- readRDS("~/Immunoediting/data/checkpoint_exp.rds")
checkpoint_ccf <- left_join(
  nes_ccf %>% mutate(sample=substr(sample,1,15)),
  checkpoint_exp
)
checkpoint_ccf$cancer <- EasyBioinfo::get_cancer_type(checkpoint_ccf$sample)
checkpoint_ccf %>% 
  group_by(cancer) %>% 
  summarise(median_es=median(es),
            median_pdl1=median(PDL1),
            median_ctla4=median(CTLA4)) -> summ1
cor.test(summ1$median_es,summ1$median_ctla4)

checkpoint_over <- readRDS("~/Immunoediting/data/checkpoint_over.rds")
over_ccf <- left_join(
  nes_ccf %>% mutate(sample=substr(sample,1,15)),
  checkpoint_over
) %>% na.omit() %>% 
  mutate(es_type=ifelse(es<0 & p <0.05,"yes","no")) %>% 
  mutate(over_type=ifelse(PDL1_over=="yes" | CTLA4_over=="yes","over","no_over"))
ggplot(data=over_ccf,aes(x=CTLA4_over,y=es))+
  geom_boxplot()+
  stat_compare_means()
table(over_ccf$es_type,over_ccf$over_type) %>% fisher.test()

apm_mut_sample <- readRDS("~/Immunoediting/data/apm_mut_sample.rds")
apm_mut_ccf <- left_join(
  nes_ccf %>% mutate(sample=substr(sample,1,15)),
  apm_mut_sample
) %>% na.omit() %>% 
  mutate(es_type=ifelse(es<0 & p <0.05,"yes","no"))
ggplot(data=apm_mut_ccf,aes(x=apm_mut,y=es))+
  geom_boxplot()+
  stat_compare_means()
table(apm_mut_ccf$es_type,apm_mut_ccf$apm_mut) %>% chisq.test()

HLA_LOH <- readRDS("~/Immunoediting/data/HLA_LOH.rds")
loh_ccf <- left_join(
  nes_ccf %>% mutate(sample=substr(sample,1,15)),
  HLA_LOH %>% rename(sample=Sample)
) %>% na.omit() %>% 
  mutate(es_type=ifelse(es<0 & p <0.05,"yes","no")) %>% 
  mutate(loh_type=ifelse(HLA_A | HLA_B | HLA_C ,"LOH","no_LOH"))
ggplot(data=loh_ccf,aes(x=loh_type,y=es))+
  geom_boxplot()+
  stat_compare_means()
table(loh_ccf$es_type,loh_ccf$loh_type) %>% fisher.test()

##rna
checkpoint_rna <- left_join(
  nes_exp %>% mutate(sample=substr(sample,1,15)),
  checkpoint_exp
)
checkpoint_ccf$cancer <- EasyBioinfo::get_cancer_type(checkpoint_ccf$sample)
checkpoint_ccf %>% 
  group_by(cancer) %>% 
  summarise(median_es=median(es),
            median_pdl1=median(PDL1),
            median_ctla4=median(CTLA4)) -> summ1
cor.test(summ1$median_es,summ1$median_ctla4)

over_rna <- left_join(
  nes_exp%>% mutate(sample=substr(sample,1,15)),
  checkpoint_over
) %>% na.omit() %>% 
  mutate(es_type=ifelse(es<0 & p_value <0.05,"yes","no")) %>% 
  mutate(over_type=ifelse(PDL1_over=="yes" | CTLA4_over=="yes","over","no_over"))
ggplot(data=over_rna,aes(x=CTLA4_over,y=es))+
  geom_boxplot()+
  stat_compare_means()
table(over_rna$es_type,over_rna$over_type) %>% fisher.test()

apm_mut_exp <- left_join(
  nes_exp %>% mutate(sample=substr(sample,1,15)),
  apm_mut_sample
) %>% na.omit() %>% 
  mutate(es_type=ifelse(es<0 & p_value <0.05,"yes","no"))
ggplot(data=apm_mut_exp,aes(x=apm_mut,y=es))+
  geom_boxplot()+
  stat_compare_means()
table(apm_mut_exp$es_type,apm_mut_exp$apm_mut) %>% fisher.test()

loh_exp <- left_join(
  nes_exp %>% mutate(sample=substr(sample,1,15)),
  HLA_LOH %>% rename(sample=Sample)
) %>% na.omit() %>% 
  mutate(es_type=ifelse(es<0 & p_value <0.05,"yes","no")) %>% 
  mutate(loh_type=ifelse(HLA_A & HLA_B & HLA_C ,"LOH","no_LOH"))
ggplot(data=loh_exp,aes(x=loh_type,y=es))+
  geom_boxplot()+
  stat_compare_means()
table(loh_exp$es_type,loh_exp$loh_type) %>% fisher.test()







#--------------------------

# 不同突变类型的CCF分布 ------------------------------------------------------------
all_mut_mis_ccf <- readRDS("~/Immunoediting/data/all_mut_mis_ccf.rds")
pancancer_mutation <- readRDS("~/Immunoediting/data/pancancer_mutation.rds")
pancancer_mutation <- pancancer_mutation %>% 
  mutate(index=paste(gene,Sample_ID,chrom,start,ref,alt,sep = ":")) %>% 
  select(Sample_ID,Amino_Acid_Change,index)
pancancer_mutation_ccf <- left_join(
  pancancer_mutation %>% select(-Sample_ID),
  all_mut_mis_ccf %>% mutate(index=paste(Hugo_Symbol,index,sep = ":")) %>% select(index,ccf_hat),by="index"
) %>% distinct(index,.keep_all = T) %>% na.omit() %>% 
  tidyr::separate(col = index,into = c("gene","sample","chromosome","position","ref","alt"),sep = ":")

pancancer_mutation_ccf_mt <- get_mutation_type(pancancer_mutation_ccf)
pancancer_mutation_ccf_mt <- as.data.frame(pancancer_mutation_ccf_mt)

##验证是不是两两分布都没有差异
dt <- matrix(c(1:96*96),nrow = 96,ncol=96)
dt <- as.data.frame(dt)
colnames(dt) <- unique(pancancer_mutation_ccf_mt$mutation_type)
rownames(dt) <- colnames(dt)

for (i in rownames(dt)){
  for (j in colnames(dt)){
    a <- pancancer_mutation_ccf_mt %>% 
      filter(mutation_type==i)
    b <- pancancer_mutation_ccf_mt %>% 
      filter(mutation_type==j)
    dt[i,j] <- ks.test(a$ccf_hat,b$ccf_hat)$p.value
  }
}

ggplot(data=a,aes(x=ccf_hat,..scaled..))+
  geom_density()
ggplot(data=b,aes(x=ccf_hat,..scaled..))+
  geom_density()

##检测审稿人提到的原癌基因和抑癌基因中的突变的CCF分布差别
oncogene <- data.table::fread("data/ongene_human.txt",data.table = F)
tsg <- data.table::fread("data/Human_TSGs.txt",data.table = F)

oncogene_mt <- pancancer_mutation_ccf_mt %>% filter(gene %in% oncogene$OncogeneName)
tsg_mt <- pancancer_mutation_ccf_mt %>% filter(gene %in% tsg$GeneSymbol)
ks.test(oncogene_mt$ccf_hat,tsg_mt$ccf_hat)

library(ggprism)
dt <- data.frame(ccf=c(oncogene_mt$ccf_hat,tsg_mt$ccf_hat),
                 type=c(rep("Oncogenic gene",nrow(oncogene_mt)),rep("Suppressor gene",nrow(tsg_mt))))
ggplot(data=dt,aes(x=ccf,..scaled..,color=type))+
  geom_density()+
  theme_prism()+
  labs(x="CCF",y="Density")+
  annotate("text", x=0.50, y=0.75, label="Kolmogorov-Smirnov test, p-value = 0.98",size=7)





#---------------------------

# 把强的 Escape 信号的样本排除掉，然后再看 ESCCF 的信号 --------------------------------------
nes_ccf <- readRDS("data/neo_nes_ccf06_1_remove_driver_samples_addp.rds")
nes_exp <- readRDS("data/nes_exp_addp.rds")

both <- left_join(
  nes_ccf,
  nes_exp %>% select(es,sample,p_value) %>% rename(es_exp=es,p_exp=p_value)
) %>% na.omit() %>% 
  mutate(exp_type=ifelse(es_exp<0 & p_exp<0.05,"escape","no_escape"))

median(both[both$exp_type=="no_escape","es"])
median(both[both$exp_type=="escape","es"])

ggplot(data=both,aes(x=exp_type,y=es))+
  geom_boxplot()+
  stat_compare_means()

filter_escape <- both %>% filter(exp_type=="no_escape")
saveRDS(filter_escape,file = "data/filter_escape.rds")
##使用这些样本进行模拟sim_filter_escape_sample_and_driver.R
neo_nes <- readRDS("data/filter_escape.rds")
sim_all <- readRDS("~/data/sim_2000_not_filter_all_ccf.rds")
sim_no_escape <- sim_all %>% filter(sample %in% neo_nes$sample)

sim_no_escape %>%
  group_by(cancer,sim_num) %>%
  summarise(median_es=median(es)) -> summ2
neo_nes$cancer <- EasyBioinfo::get_cancer_type(neo_nes$sample)
neo_nes_summ <- neo_nes %>% 
  group_by(cancer) %>% summarise(median_es=median(es))
sim_all %>%
  group_by(sim_num) %>%
  summarise(median_es=median(es)) -> summ
p <- WVPlots::ShadedDensity(frame = summ, 
                            xvar = "median_es",
                            threshold = median(neo_nes$es),
                            title = "",
                            tail = "left")

p$layers[[1]]$aes_params$colour <- "red"
p$layers[[1]]$aes_params$size <- 1
p$layers[[2]]$aes_params$fill <- "blue"     #geom_ribbon
p$layers[[3]]$aes_params$colour <- "black" 
p$layers[[3]]$aes_params$size <- 1
#p$layers[[4]]$aes_params$label <- "Actual median ES" #geom_text
p1 <- p + labs(x="Simulation median es")+
  theme_prism()+
  labs(y="Density")
p1

neo_nes_summ <- neo_nes_summ %>%
  rowwise() %>%
  mutate(p=mean(summ2$median_es[summ2$cancer==cancer] <= median_es))
saveRDS(neo_nes_summ,file = "data/neo_sim_filter_excapeAnddriver.rds")

get_f1 <- function(dt,pancancer_p,dt2,median_es){
  p1 <- ggplot(data=dt,aes(x=1,y=es))+
    geom_violin(alpha=0.7,width=0.5)+
    geom_boxplot(width=0.2)+
    theme_prism()+
    labs(x=paste0("median es = ",round(median_es,digits = 3),"\n n = ",nrow(dt)),size=4,y="ESccf")+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank())+
    theme(text = element_text(size = 7))+
    annotate(geom="text", x=1, y=1.1, label=paste0("p = ",pancancer_p),
             color="red",size=4)
  dt$cancer <- get_cancer_type(dt$sample)
  dt %>% group_by(cancer) %>%
    summarise(median_est=median(es),c=n()) %>%
    arrange(median_est) %>%
    mutate(label=paste0(cancer,"\n(n=",c,")"))-> summ1
  summ1 <- left_join(summ1,dt2 %>% select(cancer,p))
  summ1 <- summ1 %>%
    mutate(sig=case_when(
      p < 0.05 & p > 0.01 ~ "*",
      p < 0.01 ~ "**",
      TRUE ~ "ns"
    ))
  dt <- left_join(dt,summ1)
  dt$label <- factor(dt$label,levels = summ1$label)
  df2 <- data.frame(x = 1:nrow(summ1), y = 1.1, family = summ1$sig)
  p2 <- ggplot(data=dt,aes(x=label,y=es))+
    geom_boxplot()+
    theme_prism()+
    labs(y="ESccf")+
    theme(axis.title.x = element_blank())+
    theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1))+
    theme(text = element_text(size = 6))+
    geom_text(data=df2,aes(x=x,y=y,label = family))+
    geom_hline(yintercept=0,
               color = "red", size=1)
  p3 <- p1 + p2 + plot_layout(widths = c(1, 8))
}


neo_nes <- readRDS("data/filter_escape.rds") %>% select(-p)
neo_nes_summ <- readRDS("data/neo_sim_filter_excapeAnddriver.rds")
p3 <- get_f1(neo_nes,pancancer_p = 0,dt2 = neo_nes_summ,median_es = median(neo_nes$es))
p3
##考虑其他的逃逸因素,分别做逃逸和不逃逸的模拟
checkpoint_over <- readRDS("~/Immunoediting/data/checkpoint_over.rds")
apm_mut_sample <- readRDS("~/Immunoediting/data/apm_mut_sample.rds")
#HLA_LOH <- readRDS("~/Immunoediting/data/HLA_LOH.rds")

inter_sample <- Reduce(intersect, list(checkpoint_over$sample,apm_mut_sample$sample,
                                       substr(both$sample,1,15)))
both_other <- Reduce(left_join,list(both %>% mutate(sample=substr(sample,1,15)) %>% 
                                      dplyr::filter(sample %in% inter_sample),
                                    checkpoint_over %>% select(sample,PDL1_over,CTLA4_over) %>% 
                                      dplyr::filter(sample %in% inter_sample),
                                    apm_mut_sample %>% dplyr::filter(sample %in% inter_sample)))
both_other <- both_other %>% 
  mutate(escape=ifelse(exp_type=="escape" | PDL1_over == "yes" | CTLA4_over == "yes" | apm_mut,"yes","no"))

ggplot(data=both_other,aes(x=escape,y=es))+
  geom_boxplot()+
  stat_compare_means()+
  labs(x="Escape",y="ESccf")+
  theme_prism()
median(both_other[both_other$escape=="no","es"])
median(both_other[both_other$escape=="yes","es"])
table(both_other$escape)
both_other$cancer <- get_cancer_type(both_other$sample)
table(both_other$escape,both_other$cancer)
saveRDS(both_other,file = "data/all_escape_mech.rds")

neo_nes <- readRDS("~/Immunoediting/data/all_escape_mech.rds")
neo_nes$cancer <- EasyBioinfo::get_cancer_type(neo_nes$sample)
escape <- neo_nes %>% filter(escape=="yes")
no_escape <- neo_nes %>% filter(escape!="yes")
sim_escape <- sim_all %>% 
  filter(substr(sample,1,15) %in% escape$sample)
sim_no_escape <- sim_all %>% 
  filter(substr(sample,1,15) %in% no_escape$sample)

sim_escape %>%
  group_by(cancer,sim_num) %>%
  summarise(median_es=median(es)) -> summ1
escape_summ <- escape %>% 
  group_by(cancer) %>% summarise(median_es=median(es))
sim_escape %>%
  group_by(sim_num) %>%
  summarise(median_es=median(es)) -> summ2
p <- WVPlots::ShadedDensity(frame = summ2, 
                            xvar = "median_es",
                            threshold = median(escape$es),
                            title = "",
                            tail = "left")

p$layers[[1]]$aes_params$colour <- "red"
p$layers[[1]]$aes_params$size <- 1
p$layers[[2]]$aes_params$fill <- "blue"     #geom_ribbon
p$layers[[3]]$aes_params$colour <- "black" 
p$layers[[3]]$aes_params$size <- 1
p1 <- p + labs(x="Simulation median es")+
  theme_prism()+
  labs(y="Density")
p1

noescape_summ <- no_escape %>% 
  group_by(cancer) %>% summarise(median_es=median(es))
sim_no_escape %>%
  group_by(sim_num) %>%
  summarise(median_es=median(es)) -> summ3
p <- WVPlots::ShadedDensity(frame = summ3, 
                            xvar = "median_es",
                            threshold = median(no_escape$es),
                            title = "",
                            tail = "left")

p$layers[[1]]$aes_params$colour <- "red"
p$layers[[1]]$aes_params$size <- 1
p$layers[[2]]$aes_params$fill <- "blue"     #geom_ribbon
p$layers[[3]]$aes_params$colour <- "black" 
p$layers[[3]]$aes_params$size <- 1
p2 <- p + labs(x="Simulation median es")+
  theme_prism()+
  labs(y="Density")
p2

p1+p2

sim_no_escape %>%
  group_by(cancer,sim_num) %>%
  summarise(median_es=median(es)) -> summ4

noescape_summ <- noescape_summ %>%
  rowwise() %>%
  mutate(p=mean(summ4$median_es[summ4$cancer==cancer] <= median_es))
escape_summ <- escape_summ %>%
  rowwise() %>%
  mutate(p=mean(summ1$median_es[summ1$cancer==cancer] <= median_es))
saveRDS(escape_summ,file = "data/neo_sim_filter_driver_escape.rds")
saveRDS(noescape_summ,file = "data/neo_sim_filter_driver_noescape.rds")


escape <- escape %>% select(-p)
no_escape <- no_escape %>% select(-p)
p3 <- get_f1(escape,pancancer_p = 0.74,dt2 = escape_summ,median_es = median(escape$es))
p3

p4 <- get_f1(no_escape,pancancer_p = 0.014,dt2 = noescape_summ,median_es = median(no_escape$es))
p4

p3/p4
#--------------------------

# es和突变数量的关系 --------------------------------------------------------------
##CCF
all_mut_ccf <- readRDS("data/all_mut_ccf_tpm.rds")
all_mut_ccf <- all_mut_ccf %>%
  rename(ccf=ccf_hat) %>%
  mutate(neo=ifelse(neo=="neo","yes","no"))
samples_has_subclonal <- all_mut_ccf %>% filter(ccf<0.6) %>% select(sample) %>%
  distinct(sample)
all_mut_ccf <- all_mut_ccf %>%
  mutate(gene_protein_change=paste(Hugo_Symbol,Protein_Change,sep = "-"))
driver_mutations <- readRDS("~/Immunoediting/data/driver_mutations.rds")
driver_mutations <- driver_mutations %>%
  mutate(gene_protein_change=paste(gene,protein_change,sep = "-"))
all_mut_ccf <- all_mut_ccf %>%
  mutate(is_driver=ifelse(gene_protein_change %in% driver_mutations$gene_protein_change,"yes","no"))
sum(all_mut_ccf$is_driver=="yes" & all_mut_ccf$neo=="yes") / sum(all_mut_ccf$neo=="yes")
all_mut_ccf %>% group_by(sample) %>%
  summarise(inter_gene=intersect(Hugo_Symbol[neo=="yes"],
                                 Hugo_Symbol[is_driver=="yes"])) -> aaa
all_mut_ccf <- all_mut_ccf %>%
  mutate(sample_neo_index=paste(sample,neo,Hugo_Symbol,sep = ","))
aaa <- aaa %>% mutate(sample_neo_index=paste(sample,"yes",inter_gene,sep = ","))
all_mut_ccf %>%
  mutate(in_aaa = ifelse(sample_neo_index %in% aaa$sample_neo_index,"yes","no")) %>%
  group_by(sample) %>%
  summarise(need_sample=ifelse(any(in_aaa=="yes"),"no","yes")) %>%
  filter(need_sample=="yes") -> summ2
need_samples <- intersect(samples_has_subclonal$sample,summ2$sample)
all_mut_ccf  %>%
  filter(sample %in% need_samples) %>%
  group_by(sample) %>%
  summarise(c_n=sum(neo=="yes"),c_m=sum(neo=="no")) %>% filter(c_n>=1 & c_m >=1) -> summ
neo_missense <- all_mut_ccf %>% filter(sample %in% summ$sample)
neo_missense <- neo_missense %>% select(sample,neo,ccf) %>% filter(!is.na(ccf))
saveRDS(neo_missense,file = "data/mut_cal_es_filter_driver_samples.rds")

neo_missense <- readRDS("~/Immunoediting/data/mut_cal_es_filter_driver_samples.rds")
neo_missense %>% group_by(sample) %>% 
  summarise(neo_counts=sum(neo=="yes"),all_counts=n()) -> summ
summ <- left_join(summ,nes_ccf)

library(ggplot2)
library(ggprism)
summ <- summ %>% 
  mutate(type=ifelse(es<0 & p<0.05,"yes","no"))
ggplot(data=summ,aes(x=log(neo_counts),y=es,color=type))+
  geom_point()+
  theme_prism()+
  labs(x="log(Neoantigen counts)",y="ESccf")

summ$inter <- cut(summ$neo_counts,c(0,10,20,30,40,50,1064))
summ %>% 
  group_by(inter) %>% 
  summarise(sig_n=sum(type=="yes")/n(),
            not_sig_n=sum(type=="no")/n()) -> a
ggplot(data=a,aes(x=inter,y=sig_n))+
  geom_bar(stat = "identity")+
  theme_prism()+
  labs(x="Interval of Neoantigen Counts",y="Precent of Immune elimination sample")

##RNA
all_mut_exp <- readRDS("data/all_mut_tpm_not_filter.rds")
all_mut_exp <- all_mut_exp %>% select(sample,tpm_exp,neo) %>% rename(exp=tpm_exp)

all_mut_exp %>% group_by(sample) %>% 
  summarise(neo_counts=sum(neo=="neo"),all_counts=n()) -> summ
summ <- left_join(nes_exp,summ)
summ <- summ %>% 
  mutate(type=ifelse(es<0 & p_value<0.05,"yes","no"))
ggplot(data=summ,aes(x=log(neo_counts),y=es,color=type))+
  geom_point()+
  theme_prism()+
  labs(x="log(Neoantigen counts)",y="ESexp")

summ$inter <- cut(summ$neo_counts,c(0,10,20,30,40,50,1628))
summ %>% 
  group_by(inter) %>% 
  summarise(sig_n=sum(type=="yes")/n(),
            not_sig_n=sum(type=="no")/n()) -> a
ggplot(data=a,aes(x=inter,y=sig_n))+
  geom_bar(stat = "identity")+
  theme_prism()+
  labs(x="Interval of Neoantigen Counts",y="Precent of Immune escape sample")

###是不是需要标准化





#---------------------------

# 测序噪声影响CCF分布 -------------------------------------------------------------
all_mut_mis_ccf <- readRDS("~/Immunoediting/data/all_mut_mis_ccf.rds")
pancancer_mutation <- readRDS("~/Immunoediting/data/pancancer_mutation.rds") %>% as.data.frame()
pancancer_mutation <- pancancer_mutation %>% 
  mutate(index=paste(gene,Sample_ID,chrom,start,ref,alt,sep = ":")) %>% 
  select(Sample_ID,Amino_Acid_Change,index)
pancancer_mutation_ccf <- left_join(
  pancancer_mutation %>% select(-Sample_ID),
  all_mut_mis_ccf %>% mutate(index=paste(Hugo_Symbol,index,sep = ":")) %>% select(index,ccf_hat),by="index"
) %>% distinct(index,.keep_all = T) %>% na.omit() %>% 
  tidyr::separate(col = index,into = c("gene","sample","chromosome","position","ref","alt"),sep = ":")

TCGA_sample_depth <- readRDS("~/Immunoediting/data/TCGA_sample_depth.rds")
pancancer_mutation_ccf <- left_join(
  pancancer_mutation_ccf,TCGA_sample_depth %>% rename(sample=samples)
)
pancancer_mutation_ccf <- pancancer_mutation_ccf %>% 
  filter(!is.na(depth))
quantile(pancancer_mutation_ccf$depth)
unique(pancancer_mutation_ccf$sample) %>% length()

pancancer_mutation_ccf <- pancancer_mutation_ccf %>% 
  filter(is_wgs =="no")
pancancer_mutation_ccf <- pancancer_mutation_ccf %>% 
  filter(!is.na(depth)) %>% 
  mutate(type=case_when(
    depth < quantile(pancancer_mutation_ccf$depth)[2] ~ "low",
    depth > quantile(pancancer_mutation_ccf$depth)[4] ~ "high",
    TRUE ~ "medium"
  ))

pancancer_mutation_ccf <- pancancer_mutation_ccf %>% 
  filter(type!= "medium")
ggplot(data = pancancer_mutation_ccf,aes(x=ccf_hat,..scaled..,color=type))+
  geom_density(size=1)+
  theme_prism()+
  labs(y="Density",x="CCF")

low <- pancancer_mutation_ccf %>% filter(type=="low")
high <- pancancer_mutation_ccf %>% filter(type=="high")
ks.test(high$ccf_hat,low$ccf_hat)
#---------------------------
# 其他免疫浸润定量的方法 -------------------------------------------------------------
immune <- data.table::fread("data/infiltration_estimation_for_tcga.csv",check.names = F,data.table = F) %>%
  rename(sample=cell_type)
nes_ccf <- readRDS("data/neo_nes_ccf06_1_remove_driver_samples_addp.rds")
nes_exp <- readRDS("data/nes_exp_addp.rds")
colnames(immune)
cibersort <- immune %>% 
  select(sample,ends_with("CIBERSORT-ABS"))
colnames(cibersort)
get_immune_plot <- function(ccf_es_dt,exp_es_dt,immune_dt,sample_len,cd8_nk_names,treg_names){
  ccf_immune <- left_join(
    ccf_es_dt %>% mutate(sample=substr(sample,1,sample_len)),
    immune_dt
  ) %>%
    mutate(es_type=ifelse(es<0 & p <0.05,"yes","no"))
  
  exp_immune <- left_join(
    exp_es_dt %>% mutate(sample=substr(sample,1,sample_len)),
    immune_dt
  ) %>%
    mutate(es_type=ifelse(es<0 & p_value <0.05,"yes","no"))
  
  p1 <- ggplot(data=ccf_immune,aes(x=es_type,y=CD8_NK))+
    geom_boxplot()+
    stat_compare_means()+
    theme_prism()+
    labs(x="Elimination",y=cd8_nk_names)
  p2 <- ggplot(data=ccf_immune,aes(x=es_type,y=Treg))+
    geom_boxplot()+
    stat_compare_means()+
    theme_prism()+
    labs(x="Elimination",y=treg_names)
  
  p3 <- ggplot(data=exp_immune,aes(x=es_type,y=CD8_NK))+
    geom_boxplot()+
    stat_compare_means()+
    theme_prism()+
    labs(x="Escape",y=cd8_nk_names)
  p4 <- ggplot(data=exp_immune,aes(x=es_type,y=Treg))+
    geom_boxplot()+
    stat_compare_means()+
    theme_prism()+
    labs(x="Escape",y=treg_names)
  
  p5 <- ggplot(data=ccf_immune,aes(x=es_type,y=CD8_treg))+
    geom_boxplot()+
    stat_compare_means()+
    theme_prism()+
    labs(x="Elimination",y="CD8/Treg")
  p6 <- ggplot(data=exp_immune,aes(x=es_type,y=CD8_treg))+
    geom_boxplot()+
    stat_compare_means()+
    theme_prism()+
    labs(x="Escape",y="CD8/Treg")
  res <- list(p1,p2,p3,p4,p5,p6)
  return(res)
}

cibersort <- cibersort %>% 
  mutate(CD8_NK=`T cell CD8+_CIBERSORT-ABS`+`NK cell activated_CIBERSORT-ABS`,
         Treg=`T cell regulatory (Tregs)_CIBERSORT-ABS`,
         CD8_treg=exp(`T cell CD8+_CIBERSORT-ABS`)/exp(`T cell regulatory (Tregs)_CIBERSORT-ABS`))
res <- get_immune_plot(nes_ccf,nes_exp,cibersort,sample_len = 15,cd8_nk_names = "CD8 T + NK",treg_names = "Treg")
res[[1]] + res[[2]] + res[[3]] + res[[4]] + res[[5]] + res[[6]] + plot_layout(ncol = 2,nrow = 3)


#---------------------------


# binding_non_binding_region_ESCCF -----------------------------------------------
pro_HLAbind <- readRDS("~/data/pro_HLAbind.rds")
test <- pro_HLAbind[1:10,]
pro_HLAbind <- pro_HLAbind %>% 
  select(seqnames,pos,gene,HLA_aff_mean)
##转换为hg38
library(GenomicRanges)
pro_HLAbind$end <- pro_HLAbind$pos
pro_HLAbind_granges <- makeGRangesFromDataFrame(pro_HLAbind,
                                             keep.extra.columns=TRUE,
                                             ignore.strand=TRUE,
                                             seqinfo=NULL,
                                             seqnames.field="seqnames",
                                             start.field="pos",
                                             end.field="end",
                                             starts.in.df.are.0based=FALSE)
library(rtracklayer)
library(plyranges)
chainObject <- import.chain("~/data/hg19ToHg38.over.chain")
results <- as.data.frame(liftOver(pro_HLAbind_granges, chainObject))
results <- results %>% select(seqnames,start,gene,HLA_aff_mean)
saveRDS(results,file = "~/data/pro_HLAbinding_hg38.rds")

hla_binding <- readRDS("~/data/pro_HLAbinding_hg38.rds")
all_mut_mis_ccf <- readRDS("~/Immunoediting/data/all_mut_mis_ccf.rds")
pancancer_mutation <- readRDS("~/Immunoediting/data/pancancer_mutation.rds")
pancancer_mutation <- pancancer_mutation %>%
  mutate(index=paste(Sample_ID,chrom,start,ref,alt,sep = ":")) 

all_mut_ccf <- inner_join(
  pancancer_mutation,
  all_mut_mis_ccf %>% select(index,ccf_hat),
  by="index"
)
all_mut_ccf <- all_mut_ccf[!duplicated(all_mut_ccf$index),]
all_mut_ccf <- all_mut_ccf %>% 
  filter(!is.na(ccf_hat))

hla_binding$width <- 1
binding <- hla_binding %>% filter(HLA_aff_mean<500) %>% as_granges
binding_merge <- reduce_ranges(binding)

all_mut_ccf_o <- all_mut_ccf %>% 
  dplyr::rename(seqnames=chrom) %>% as_granges()
all_mut_ccf_o <- all_mut_ccf_o %>% mutate(n_olap = count_overlaps(., binding_merge))
all_mut_ccf_overlap <- all_mut_ccf_o %>% as.data.frame() %>% 
  mutate(neo=ifelse(n_olap==1,"yes","no")) %>% 
  dplyr::select(Sample_ID,ccf_hat,neo) %>% 
  dplyr::rename(sample=Sample_ID,ccf=ccf_hat)
saveRDS(all_mut_ccf_overlap,file = "data/hla_binding_no_binding_ccf.rds")

all_mut_ccf_overlap <- readRDS("~/Immunoediting/data/hla_binding_no_binding_ccf.rds")
samples_has_subclonal <- all_mut_ccf_overlap %>% filter(ccf<0.6) %>% select(sample) %>%
  distinct(sample)
all_mut_ccf_overlap  %>% filter(sample %in% samples_has_subclonal$sample) %>%
  group_by(sample) %>%
  summarise(c_n=sum(neo=="yes"),c_m=sum(neo=="no")) %>% filter(c_n>=1 & c_m >=1) -> summ

neo_missense <- all_mut_ccf_overlap %>% filter(sample %in% summ$sample)
cal_nes_warp <- function(dt){
  results_ccf <- vector("list",length = length(unique(dt$sample)))
  names(results_ccf) <- unique(dt$sample)
  
  cl <- makeCluster(getOption("cl.cores", 25),type="FORK")
  results_ccf <- parSapply(cl=cl,names(results_ccf),
                           function(x){
                             data <- dt %>% filter(sample == x)
                             a <- NeoEnrichment::cal_nes_new_test(dt = data,
                                                                  sample_counts = 1000,
                                                                  need_p = TRUE)
                             return(a)
                           },simplify = FALSE)
  stopCluster(cl)
  results_ccf <- Filter(function(x){length(x)>1},results_ccf)
  pancancer_nes_ccf <- bind_rows(results_ccf)
  return(pancancer_nes_ccf)
}
neo_missense <- neo_missense %>% select(sample,neo,ccf) %>% filter(!is.na(ccf))
neo_nes <- cal_nes_warp(neo_missense)
saveRDS(neo_nes,file = "data/es_hla_binding_no_binding.rds")

##sim
library(dplyr)

neo_nes <- readRDS("data/es_hla_binding_no_binding.rds")
sim_all <- readRDS("data/sim_HLA_binding_no_binding.rds")

sim_all %>%
  group_by(cancer,sim_num) %>%
  summarise(median_es=median(es)) -> summ2
#neo_nes$cancer <- EasyBioinfo::get_cancer_type(neo_nes$sample)
neo_nes_summ <- neo_nes %>% 
  group_by(cancer) %>% summarise(median_es=median(es))
sim_all %>%
  group_by(sim_num) %>%
  summarise(median_es=median(es)) -> summ
p <- WVPlots::ShadedDensity(frame = summ, 
                            xvar = "median_es",
                            threshold = median(neo_nes$es),
                            title = "",
                            tail = "left")

p$layers[[1]]$aes_params$colour <- "red"
p$layers[[1]]$aes_params$size <- 1
p$layers[[2]]$aes_params$fill <- "blue"     #geom_ribbon
p$layers[[3]]$aes_params$colour <- "black" 
p$layers[[3]]$aes_params$size <- 1
#p$layers[[4]]$aes_params$label <- "Actual median ES" #geom_text
p1 <- p + labs(x="Simulation median es")+
  theme_prism()+
  labs(y="Density")
p1

get_f1 <- function(dt,pancancer_p,dt2,median_es){
  p1 <- ggplot(data=dt,aes(x=1,y=es))+
    geom_violin(alpha=0.7,width=0.5)+
    geom_boxplot(width=0.2)+
    theme_prism()+
    labs(x=paste0("median es = ",round(median_es,digits = 3),"\n n = ",nrow(dt)),size=1)+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank())+
    theme(text = element_text(size = 7))+
    annotate(geom="text", x=1, y=1.1, label=paste0("p = ",pancancer_p),
             color="red",size=4)
  dt %>% group_by(cancer) %>%
    summarise(median_est=median(es),c=n()) %>%
    arrange(median_est) %>%
    mutate(label=paste0(cancer,"\n(n=",c,")"))-> summ1
  summ1 <- left_join(summ1,dt2 %>% select(cancer,p))
  summ1$p <- signif(summ1$p,digits = 1)
  summ1 <- summ1 %>%
    mutate(sig=case_when(
      p <= 0.05 & p > 0.01 ~ "*",
      p < 0.01 ~ "**",
      TRUE ~ "ns"
    ))
  dt <- left_join(dt,summ1)
  dt$label <- factor(dt$label,levels = summ1$label)
  df2 <- data.frame(x = 1:nrow(summ1), y = 1.1, family = summ1$sig)
  p2 <- ggplot(data=dt,aes(x=label,y=es))+
    geom_boxplot()+
    theme_prism()+
    theme(axis.title.x = element_blank())+
    theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1))+
    theme(text = element_text(size = 6))+
    geom_text(data=df2,aes(x=x,y=y,label = family))+
    geom_hline(yintercept=0,
               color = "red", size=1)
  p3 <- p1 + p2 + plot_layout(widths = c(1, 8))
  return(p1)
}
library(patchwork)
neo_nes <- neo_nes %>% select(-p)
neo_nes_summ <- neo_nes_summ %>% 
  rowwise() %>%
  mutate(p=mean(summ2$median_es[summ2$cancer==cancer] <= median_es))
p3 <- get_f1(neo_nes,pancancer_p = 0.75,dt2 = neo_nes_summ,median_es = median(neo_nes$es))
p3

#----------------------------------------------------------------------


# driver mutation cutoff --------------------------------------------------
driver <- data.table::fread("~/data/Mutation.CTAT.3D.Scores.txt")
driver_mutations <- driver %>%
  mutate(nmethod=`New_Linear (functional) flag`+`New_Linear (cancer-focused) flag`+`New_3D mutational hotspot flag`) %>% 
  filter(nmethod>=3)

all_mut_ccf <- readRDS("data/all_mut_ccf_tpm.rds")
all_mut_ccf <- all_mut_ccf %>%
  rename(ccf=ccf_hat) %>%
  mutate(neo=ifelse(neo=="neo","yes","no"))

all_mut_ccf <- all_mut_ccf %>%
  mutate(gene_protein_change=paste(Hugo_Symbol,Protein_Change,sep = "-"))
driver_mutations <- driver_mutations %>%
  mutate(gene_protein_change=paste(gene,protein_change,sep = "-"))

all_mut_ccf <- all_mut_ccf %>%
  mutate(is_driver=ifelse(gene_protein_change %in% driver_mutations$gene_protein_change,"yes","no"))
sum(all_mut_ccf$is_driver=="yes" & all_mut_ccf$neo=="yes") / sum(all_mut_ccf$neo=="yes")##0.01029878
all_mut_ccf %>% group_by(sample) %>%
  summarise(inter_gene=intersect(Hugo_Symbol[neo=="yes"],
                                 Hugo_Symbol[is_driver=="yes"])) -> aaa##701 samples
all_mut_ccf <- all_mut_ccf %>%
  mutate(sample_neo_index=paste(sample,neo,Hugo_Symbol,sep = ","))
aaa <- aaa %>% mutate(sample_neo_index=paste(sample,"yes",inter_gene,sep = ","))

samples_has_subclonal <- all_mut_ccf %>% filter(ccf<0.6) %>% select(sample) %>%
  distinct(sample)

all_mut_ccf %>%
  mutate(in_aaa = ifelse(sample_neo_index %in% aaa$sample_neo_index,"yes","no")) %>%
  group_by(sample) %>%
  summarise(need_sample=ifelse(any(in_aaa=="yes"),"no","yes")) %>%
  filter(need_sample=="yes") -> summ2
need_samples <- intersect(samples_has_subclonal$sample,summ2$sample)
all_mut_ccf  %>%
  filter(sample %in% need_samples) %>%
  group_by(sample) %>%
  summarise(c_n=sum(neo=="yes"),c_m=sum(neo=="no")) %>% filter(c_n>=1 & c_m >=1) -> summ
neo_missense <- all_mut_ccf %>% filter(sample %in% summ$sample)
neo_missense <- neo_missense %>% select(sample,neo,ccf) %>% filter(!is.na(ccf))

neo_nes_ccf06_1 <- readRDS("~/Immunoediting/data/neo_nes_ccf06_1.rds")
es_filter_1driver <- neo_nes_ccf06_1 %>% 
  filter(sample %in% unique(neo_missense$sample))

sim_all <- readRDS("~/data/sim_2000_not_filter_all_ccf.rds")
sim_filter_1driver <- sim_all %>% filter(sample %in% es_filter_1driver$sample)

sim_filter_1driver %>%
  group_by(cancer,sim_num) %>%
  summarise(median_es=median(es)) -> summ2
es_filter_1driver$cancer <- EasyBioinfo::get_cancer_type(es_filter_1driver$sample)
neo_nes_summ <- es_filter_1driver %>% 
  group_by(cancer) %>% summarise(median_es=median(es))
sim_filter_1driver %>%
  group_by(sim_num) %>%
  summarise(median_es=median(es)) -> summ
p <- WVPlots::ShadedDensity(frame = summ, 
                            xvar = "median_es",
                            threshold = median(es_filter_1driver$es),
                            title = "",
                            tail = "left")

p$layers[[1]]$aes_params$colour <- "red"
p$layers[[1]]$aes_params$size <- 1
p$layers[[2]]$aes_params$fill <- "blue"     #geom_ribbon
p$layers[[3]]$aes_params$colour <- "black" 
p$layers[[3]]$aes_params$size <- 1
#p$layers[[4]]$aes_params$label <- "Actual median ES" #geom_text
p1 <- p + labs(x="Simulation median es")+
  theme_prism()+
  labs(y="Density")
p1

neo_nes_summ <- neo_nes_summ %>%
  rowwise() %>%
  mutate(p=mean(summ2$median_es[summ2$cancer==cancer] <= median_es))

get_f1 <- function(dt,pancancer_p,dt2,median_es){
  p1 <- ggplot(data=dt,aes(x=1,y=es))+
    geom_violin(alpha=0.7,width=0.5)+
    geom_boxplot(width=0.2)+
    theme_prism()+
    labs(x=paste0("median es = ",round(median_es,digits = 3),"\n n = ",nrow(dt)),size=4,y="ESccf")+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank())+
    theme(text = element_text(size = 7))+
    annotate(geom="text", x=1, y=1.1, label=paste0("p = ",pancancer_p),
             color="red",size=4)
  dt$cancer <- get_cancer_type(dt$sample)
  dt %>% group_by(cancer) %>%
    summarise(median_est=median(es),c=n()) %>%
    arrange(median_est) %>%
    mutate(label=paste0(cancer,"\n(n=",c,")"))-> summ1
  summ1 <- left_join(summ1,dt2 %>% select(cancer,p))
  summ1 <- summ1 %>%
    mutate(sig=case_when(
      p < 0.05 & p > 0.01 ~ "*",
      p < 0.01 ~ "**",
      TRUE ~ "ns"
    ))
  dt <- left_join(dt,summ1)
  dt$label <- factor(dt$label,levels = summ1$label)
  df2 <- data.frame(x = 1:nrow(summ1), y = 1.1, family = summ1$sig)
  p2 <- ggplot(data=dt,aes(x=label,y=es))+
    geom_boxplot()+
    theme_prism()+
    labs(y="ESccf")+
    theme(axis.title.x = element_blank())+
    theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1))+
    theme(text = element_text(size = 6))+
    geom_text(data=df2,aes(x=x,y=y,label = family))+
    geom_hline(yintercept=0,
               color = "red", size=1)
  p3 <- p1 + p2 + plot_layout(widths = c(1, 8))
}

p3 <- get_f1(es_filter_1driver,pancancer_p = 0.006,dt2 = neo_nes_summ,median_es = median(es_filter_1driver$es))
p3
#-------------------------------

#--------------
# 所有癌症类型的CCF的分布以及新抗原和非新抗原的CCF分布 -------------------------------------------
all_mut_ccf <- readRDS("data/all_mut_ccf_tpm.rds") %>% as.data.frame()
all_mut_ccf <- all_mut_ccf %>%
  rename(ccf=ccf_hat) %>%
  mutate(neo=ifelse(neo=="neo","yes","no"))
es <- readRDS("~/Immunoediting/data/neo_nes_ccf06_1_remove_driver_samples_addp.rds")
neo_missense <- all_mut_ccf %>% 
  filter(sample %in% es$sample) %>% 
  select(sample,neo,ccf) %>% filter(!is.na(ccf))

neo_missense$cancer <- get_cancer_type(neo_missense$sample,parallel = T,cores = 20)

dt <- data.frame(cancer=unique(neo_missense$cancer),p=NA)
for (i in 1:nrow(dt)){
  a <- neo_missense %>% filter(cancer==dt$cancer[i])
  dt$p[i] <- ks.test(a[a$neo=="yes","ccf"],a[a$neo=="no","ccf"])$p.value
}
ggplot(data = neo_missense,aes(x=ccf,..scaled..,color=neo))+
  geom_density(size=2)+
  theme_prism()+
  labs(y="Density",x="CCF")+
  facet_wrap(. ~ cancer)

test_dep <- data.table::fread("~/data/TCGA_bai/sample-output.mosdepth.global.dist.txt",data.table = F)
test_dep <- test_dep %>% 
  rowwise() %>% 
  mutate(ll=nchar(V1)) %>% 
  filter(ll==4)
ggplot(data=test_dep,aes(x=V2,y=V3))+
  geom_line()
#--------------

#----------------
# MArty et al的方法，计算新抗原的几何平均 ; 即6个HLA最低的IC50的几何平均-----------------------------------------------
library("psych") 
library(tidyr)
files <- list.files("~/data/neo_pre/")
re <- vector("list",100)

##先对每个HLA选择其对应的肽中novlity==1的肽的最低的IC50，如果没有等于1的，就返回NA

get_minIC50 <- function(ic50,novelty){
  dt <- data.frame(x=ic50,y=novelty)
  return(ifelse(all(dt$y==0),NA,min(dt[dt$y==1,"x"])))
}
library(foreach)
library(parallel)
numCores <- detectCores()
cl <- parallel::makeCluster(15)
doParallel::registerDoParallel(cl)

all_mut <- foreach(i = files,#要合并的多个文件路径 
        .combine = 'bind_rows',.packages = c("dplyr","tidyr","psych"),
        .verbose=TRUE) %dopar% {
          mut <- readRDS(paste("~/data/neo_pre/",i,sep = ""))
          mut_summ <- mut %>% 
            tidyr::separate(gene,into = "gene",sep = ":") %>% 
            mutate(index=paste(index,gene,sep = ":")) %>% 
            group_by(index,hla) %>% 
            summarise(min_ic50=get_minIC50(ic50 = IC50,novelty = novelty)) %>% 
            group_by(index) %>% 
            summarise(arr_ic50=ifelse(all(is.na(min_ic50)),50000,harmonic.mean(min_ic50,na.rm = T))) %>% 
            mutate(neo=ifelse(arr_ic50<50,"neo","not_neo")) %>% 
            tidyr::separate(index,into = c("sample","chr","position","ref","alt","gene"),sep = ":")
          return(mut_summ)
        }

saveRDS(all_mut,file="~/data/all_mut_harmonic1.rds")

all_mut <- readRDS("~/data/all_mut_harmonic1.rds")
tpm <- readRDS("~/data/tpm_trans.rds")
tpm <- tpm[!duplicated(tpm$gene),]
all_mut$tpm_exp <- mapply(function(sample,gene){
  tpm[tpm$gene==gene,substr(sample,1,15)]
},all_mut$sample,all_mut$gene)
table(lengths(all_mut$tpm_exp))
all_mut1 <- all_mut %>%
  filter(lengths(tpm_exp)!=0)
all_mut1$tpm_exp <- as.numeric(all_mut1$tpm_exp)

all_mut1 <- all_mut1 %>%
  mutate(neo2=ifelse(neo=="neo" & tpm_exp>1,"neo","not_neo"))
saveRDS(all_mut1,file = "data/all_mut_tpm_harmonic1.rds")

results <- readRDS("data/all_mut_mis_ccf.rds")
all_mut_tpm <- readRDS("data/all_mut_tpm_harmonic1.rds")
all_mut <- all_mut_tpm %>%
  mutate(index=paste(sample,chr,position,ref,alt,sep = ":")) %>%
  dplyr::rename(ref_allele=ref,alt_allele=alt)

all_mut_ccf <- inner_join(
  all_mut,
  results %>% select(-sample),
  by="index"
)
all_mut_ccf <- all_mut_ccf[!duplicated(all_mut_ccf$index),]
saveRDS(all_mut_ccf,file = "data/all_mut_ccf_tpm_harmonic1.rds")


all_mut_ccf <- readRDS("data/all_mut_ccf_tpm_harmonic1.rds")
all_mut_ccf <- all_mut_ccf %>%
  rename(ccf=ccf_hat) %>% 
  mutate(neo=ifelse(neo2=="neo","yes","no"))
mean(all_mut_ccf$neo=="yes")
samples_has_subclonal <- all_mut_ccf %>% filter(ccf<0.6) %>% select(sample) %>%
  distinct(sample)
all_mut_ccf  %>% filter(sample %in% samples_has_subclonal$sample) %>%
  group_by(sample) %>%
  summarise(c_n=sum(neo=="yes"),c_m=sum(neo=="no")) %>% filter(c_n>=1 & c_m >=1) -> summ
neo_missense <- all_mut_ccf %>% filter(sample %in% summ$sample)
cal_nes_warp <- function(dt){
  results_ccf <- vector("list",length = length(unique(dt$sample)))
  names(results_ccf) <- unique(dt$sample)
  
  cl <- makeCluster(getOption("cl.cores", 35),type="FORK")
  results_ccf <- parSapply(cl=cl,names(results_ccf),
                           function(x){
                             data <- dt %>% filter(sample == x)
                             a <- NeoEnrichment::cal_nes_new_test(dt = data,
                                                                  sample_counts = 1000,
                                                                  need_p = FALSE)
                             return(a)
                           },simplify = FALSE)
  stopCluster(cl)
  results_ccf <- Filter(function(x){length(x)>1},results_ccf)
  pancancer_nes_ccf <- bind_rows(results_ccf)
  return(pancancer_nes_ccf)
}
neo_missense <- neo_missense %>% select(sample,neo,ccf) %>% filter(!is.na(ccf))
neo_nes <- cal_nes_warp(neo_missense)
median(neo_nes$es)
saveRDS(neo_nes,file = "data/neo_nes_ccf06_1_harmonic1.rds")


##sim
neo_nes <- readRDS("data/neo_nes_ccf06_1_harmonic1.rds")
all_mut_ccf <- readRDS("data/all_mut_ccf_tpm_harmonic1.rds")
all_mut_ccf <- all_mut_ccf %>%
  dplyr::rename(ccf=ccf_hat) %>%
  mutate(neo=ifelse(neo2=="neo","yes","no"))
mean(all_mut_ccf$neo=="yes")
neo_missense <- all_mut_ccf %>% filter(sample %in% neo_nes$sample)
neo_missense <- neo_missense %>% select(sample,neo,ccf) %>% filter(!is.na(ccf))

cal_nes_warp <- function(dt){
  results_ccf <- vector("list",length = length(unique(dt$sample)))
  names(results_ccf) <- unique(dt$sample)
  
  cl <- makeCluster(getOption("cl.cores", 35),type="FORK")
  results_ccf <- parSapply(cl=cl,names(results_ccf),
                           function(x){
                             data <- dt %>% filter(sample == x)
                             a <- NeoEnrichment::cal_nes_new_test(dt = data,
                                                                  sample_counts = 1000,
                                                                  need_p = FALSE)
                             return(a)
                           },simplify = FALSE)
  stopCluster(cl)
  results_ccf <- Filter(function(x){length(x)>1},results_ccf)
  pancancer_nes_ccf <- bind_rows(results_ccf)
  return(pancancer_nes_ccf)
}
res <- vector("list",2000)
for (i in 1:2000){
  neo_missense_sim <- neo_missense %>%
    group_by(sample) %>%
    mutate(neo_sim=sample(neo,length(neo)))
  neo_missense_sim <- neo_missense_sim %>%
    select(-neo) %>%
    rename(neo=neo_sim)
  neo_nes_sim <- cal_nes_warp(neo_missense_sim)
  neo_nes_sim$sim_num <- i
  res[[i]] <- neo_nes_sim
  print(paste0("Complete ",i," sim. "))
}
saveRDS(res,file = "~/data/sim_2000_not_filter_driver_harmonic1.rds")

neo_nes <- readRDS("data/neo_nes_ccf06_1_harmonic1.rds")
sim_all <- readRDS("~/data/sim_2000_not_filter_driver_harmonic1.rds")
sim_all <- bind_rows(sim_all)
sim_all$cancer <- EasyBioinfo::get_cancer_type(sim_all$sample,parallel = T,cores = 30)
sim_all %>%
  group_by(cancer,sim_num) %>%
  summarise(median_es=median(es)) -> summ2
neo_nes$cancer <- EasyBioinfo::get_cancer_type(neo_nes$sample)
neo_nes_summ <- neo_nes %>% 
  group_by(cancer) %>% summarise(median_es=median(es))
sim_all %>%
  group_by(sim_num) %>%
  summarise(median_es=median(es)) -> summ
p <- WVPlots::ShadedDensity(frame = summ, 
                            xvar = "median_es",
                            threshold = median(neo_nes$es),##-0.01888807
                            title = "",
                            tail = "left")

p$layers[[1]]$aes_params$colour <- "red"
p$layers[[1]]$aes_params$size <- 1
p$layers[[2]]$aes_params$fill <- "blue"     #geom_ribbon
p$layers[[3]]$aes_params$colour <- "black" 
p$layers[[3]]$aes_params$size <- 1
#p$layers[[4]]$aes_params$label <- "Actual median ES" #geom_text
p1 <- p + labs(x="Simulation median es")+
  theme_prism()+
  labs(y="Density")
p1

all_mut_ccf <- readRDS("data/all_mut_ccf_tpm_harmonic1.rds")
all_mut_ccf <- all_mut_ccf %>%
  dplyr::rename(ccf=ccf_hat) %>%
  mutate(neo=ifelse(neo2=="neo","yes","no"))
all_mut_ccf <- all_mut_ccf %>%
  mutate(gene_protein_change=paste(Hugo_Symbol,Protein_Change,sep = "-"))
driver_mutations <- readRDS("~/Immunoediting/data/driver_mutations.rds")
driver_mutations <- driver_mutations %>%
  mutate(gene_protein_change=paste(gene,protein_change,sep = "-"))

all_mut_ccf <- all_mut_ccf %>%
  mutate(is_driver=ifelse(gene_protein_change %in% driver_mutations$gene_protein_change,"yes","no"))
sum(all_mut_ccf$is_driver=="yes" & all_mut_ccf$neo=="yes") / sum(all_mut_ccf$neo=="yes")
all_mut_ccf %>% group_by(sample) %>%
  summarise(inter_gene=intersect(Hugo_Symbol[neo=="yes"],
                                 Hugo_Symbol[is_driver=="yes"])) -> aaa##701 samples
all_mut_ccf <- all_mut_ccf %>%
  mutate(sample_neo_index=paste(sample,neo,Hugo_Symbol,sep = ","))
aaa <- aaa %>% mutate(sample_neo_index=paste(sample,"yes",inter_gene,sep = ","))

samples_has_subclonal <- all_mut_ccf %>% filter(ccf<0.6) %>% select(sample) %>%
  distinct(sample)
all_mut_ccf %>%
  mutate(in_aaa = ifelse(sample_neo_index %in% aaa$sample_neo_index,"yes","no")) %>%
  group_by(sample) %>%
  summarise(need_sample=ifelse(any(in_aaa=="yes"),"no","yes")) %>%
  filter(need_sample=="yes") -> summ2
need_samples <- intersect(samples_has_subclonal$sample,summ2$sample)
all_mut_ccf  %>%
  filter(sample %in% need_samples) %>%
  group_by(sample) %>%
  summarise(c_n=sum(neo=="yes"),c_m=sum(neo=="no")) %>% filter(c_n>=1 & c_m >=1) -> summ
neo_missense <- all_mut_ccf %>% filter(sample %in% summ$sample)
neo_missense <- neo_missense %>% select(sample,neo,ccf) %>% filter(!is.na(ccf))
neo_nes <- cal_nes_warp(neo_missense)##addp
median(neo_nes$es)
saveRDS(neo_nes,file = "data/neo_nes_ccf06_1_harmonic_filter_driver1.rds")
###
neo_nes <- readRDS("data/neo_nes_ccf06_1_harmonic_filter_driver1.rds")
sim_all <- readRDS("~/data/sim_2000_not_filter_driver_harmonic.rds")
sim_all <- bind_rows(sim_all)
sim_all <- sim_all %>% dplyr::filter(sample %in% neo_nes$sample)
sim_all$cancer <- EasyBioinfo::get_cancer_type(sim_all$sample,parallel = T,cores = 30)
sim_all %>%
  group_by(cancer,sim_num) %>%
  summarise(median_es=median(es)) -> summ2
neo_nes$cancer <- EasyBioinfo::get_cancer_type(neo_nes$sample)
neo_nes_summ <- neo_nes %>% 
  group_by(cancer) %>% summarise(median_es=median(es))
sim_all %>%
  group_by(sim_num) %>%
  summarise(median_es=median(es)) -> summ
p <- WVPlots::ShadedDensity(frame = summ, 
                            xvar = "median_es",
                            threshold = median(neo_nes$es),##-0.0217571
                            title = "",
                            tail = "left")

p$layers[[1]]$aes_params$colour <- "red"
p$layers[[1]]$aes_params$size <- 1
p$layers[[2]]$aes_params$fill <- "blue"     #geom_ribbon
p$layers[[3]]$aes_params$colour <- "black" 
p$layers[[3]]$aes_params$size <- 1
#p$layers[[4]]$aes_params$label <- "Actual median ES" #geom_text
p2 <- p + labs(x="Simulation median es")+
  theme_prism()+
  labs(y="Density")
p2

library(patchwork)
p1+p2
#---------------------------------------------------------

