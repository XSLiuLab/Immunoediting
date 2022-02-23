library(dplyr)
library(ggplot2)
library(ggprism)
library(NeoEnrichment)
library(ggpubr)
library(patchwork)
es_exp <- readRDS("data/es_exp_filter_driver.rds")
es_ccf <- readRDS("data/neo_nes_ccf06_1_remove_driver_samples_addp.rds")
pancancer_subtcells <- readRDS("data/pancancer_subtcells.rds")
exp_immune <- left_join(
  es_exp,
  pancancer_subtcells
) %>%
  mutate(escape=ifelse(es<0 & p_value <0.05,"yes","no")) %>% na.omit(.)
anno <- data.table::fread("data/annotation-tcga.tsv",data.table = F)
exp_immune_anno <- left_join(
  exp_immune %>% mutate(sample=substr(sample,1,12)),
  anno %>% rename(sample=V1)
)
exp_MSI <- exp_immune_anno %>% filter(!is.na(MSI)) %>% 
  filter(nchar(MSI)>=1)

p0 <- ggplot(data=exp_immune,aes(x=escape,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="CD8 T + NK")
p01 <- ggplot(data=exp_immune,aes(x=escape,y=iTreg+nTreg))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="iTreg + nTreg")
p02 <- ggplot(data=exp_immune,aes(x=escape,y=exp(CD8_T)/exp(iTreg+nTreg)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="CD8_T / Treg")
p0+p01+p02

p1 <- ggplot(data=exp_MSI,aes(x=escape,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="CD8 T + NK")+
  facet_wrap(~MSI)
p11 <- ggplot(data=exp_MSI,aes(x=escape,y=iTreg+nTreg))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="iTreg + nTreg")+
  facet_wrap(~MSI)
p12 <- ggplot(data=exp_MSI,aes(x=escape,y=exp(CD8_T)/exp(iTreg+nTreg)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="CD8_T / Treg")+
  facet_wrap(~MSI)
p1 + p11 + p12
(p0 + p01 + p02)/(p1 + p11 + p12)

ccf_immune <- left_join(
  es_ccf,
  pancancer_subtcells
) %>%
  mutate(es_type=ifelse(es<0 & p <0.05,"yes","no")) %>% na.omit(.)
ccf_immune_anno <- left_join(
  ccf_immune %>% mutate(sample=substr(sample,1,12)),
  anno %>% rename(sample=V1)
)
ccf_MSI <- ccf_immune_anno %>% filter(!is.na(MSI)) %>% 
  filter(nchar(MSI)>=1)
p0 <- ggplot(data=ccf_immune,aes(x=es_type,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8 T + NK")
p01 <- ggplot(data=ccf_immune,aes(x=es_type,y=iTreg+nTreg))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="iTreg + nTreg")
p02 <- ggplot(data=ccf_immune,aes(x=es_type,y=exp(CD8_T)/exp(iTreg+nTreg)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8_T / Treg")
p0+p01+p02

p1 <- ggplot(data=ccf_MSI,aes(x=es_type,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8 T + NK")+
  facet_wrap(~MSI)
p11 <- ggplot(data=ccf_MSI,aes(x=es_type,y=iTreg+nTreg))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="iTreg + nTreg")+
  facet_wrap(~MSI)
p12 <- ggplot(data=ccf_MSI,aes(x=es_type,y=exp(CD8_T)/exp(iTreg+nTreg)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8_T / Treg")+
  facet_wrap(~MSI)
p1 + p11 + p12
(p0 + p01 + p02)/(p1 + p11 + p12)

##cancer type
ccf_immune$cancer <- get_cancer_type(ccf_immune$sample)
ccf_immune %>% 
  group_by(cancer,es_type) %>% 
  summarise(c=n()) %>% 
  group_by(cancer) %>% 
  summarise(need=ifelse(all(c>3) & length(c)==2,"yes","no")) %>% filter(need=="yes")-> summ1
ccf_immune_filter <- ccf_immune %>% 
  filter(cancer %in% summ1$cancer)
p0 <- ggplot(data=ccf_immune_filter,aes(x=es_type,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8 T + NK")+
  facet_wrap(~cancer)
p01 <- ggplot(data=ccf_immune_filter,aes(x=es_type,y=iTreg+nTreg))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="iTreg + nTreg")+
  facet_wrap(~cancer)
p02 <- ggplot(data=ccf_immune_filter,aes(x=es_type,y=exp(CD8_T)/exp(iTreg+nTreg)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8_T / Treg")+
  facet_wrap(~cancer)
p02

exp_immune$cancer <- get_cancer_type(exp_immune$sample)
exp_immune %>% 
  group_by(cancer,escape) %>% 
  summarise(c=n()) %>% 
  group_by(cancer) %>% 
  summarise(need=ifelse(all(c>3) & length(c)==2,"yes","no")) %>% filter(need=="yes")-> summ1
exp_immune_filter <- exp_immune %>% 
  filter(cancer %in% summ1$cancer)
p0 <- ggplot(data=exp_immune_filter,aes(x=escape,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="CD8 T + NK")+
  facet_wrap(~cancer)
p01 <- ggplot(data=exp_immune_filter,aes(x=escape,y=iTreg+nTreg))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="iTreg + nTreg")+
  facet_wrap(~cancer)
p02 <- ggplot(data=exp_immune_filter,aes(x=escape,y=exp(CD8_T)/exp(iTreg+nTreg)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="CD8_T / Treg ")+
  facet_wrap(~cancer)
p02

##更新S12-S14
p1 <- ggplot(data=ccf_immune,aes(x=es_type,y=exp(CD8_T)/exp(iTreg+nTreg)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8_T / Treg")
p2 <- ggplot(data=ccf_immune,aes(x=es_type,y=exp(NK)/exp(iTreg+nTreg)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="NK / Treg")
p3 <- ggplot(data=exp_immune,aes(x=escape,y=exp(CD8_T)/exp(iTreg+nTreg)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="CD8_T / Treg")
p4 <- ggplot(data=exp_immune,aes(x=escape,y=exp(NK)/exp(iTreg+nTreg)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="NK / Treg")
(p1 +p2)/(p3+p4)

immune <- data.table::fread("data/infiltration_estimation_for_tcga.csv",check.names = F,data.table = F) %>%
  rename(sample=cell_type)
cibersort <- immune %>% 
  select(sample,ends_with("CIBERSORT-ABS"))
colnames(immune)
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
res <- get_immune_plot(es_ccf,es_exp,cibersort,sample_len = 15,cd8_nk_names = "CD8 T + NK",treg_names = "Treg")
res[[1]] + res[[2]] + res[[3]] + res[[4]] + res[[5]] + res[[6]] + plot_layout(ncol = 2,nrow = 3)

quantiseq <- immune %>% 
  select(sample,ends_with("QUANTISEQ"))
colnames(quantiseq)
quantiseq <- quantiseq %>% 
  mutate(CD8_NK=`T cell CD8+_QUANTISEQ`+`NK cell_QUANTISEQ`,
         Treg=`T cell regulatory (Tregs)_QUANTISEQ`,
         CD8_treg=exp(`T cell CD8+_QUANTISEQ`)/exp(`T cell regulatory (Tregs)_QUANTISEQ`))
res <- get_immune_plot(es_ccf,es_exp,quantiseq,sample_len = 15,cd8_nk_names = "CD8 T + NK",treg_names = "Treg")
res[[1]] + res[[2]] + res[[3]] + res[[4]] + res[[5]] + res[[6]] + plot_layout(ncol = 2,nrow = 3)

####
neo_missense <- readRDS("../tmp/mut_cal_es_filter_driver_samples.rds")
ccf_neo <- readRDS("../tmp/ccf_neo_mut_counts.rds")
neo_mut_sig <- ccf_neo %>% filter(p<0.05)
neo_missense_summ <- neo_missense %>% 
  group_by(sample) %>% 
  summarise(neo_c=sum(neo=="yes"),not_neo_c=sum(neo=="no")) %>% 
  mutate(type=case_when(
    neo_c==1 & not_neo_c >=19 ~ "a",
    neo_c==2 & not_neo_c >=5 ~ "a",
    neo_c==3 & not_neo_c >=4 ~ "a",
    neo_c==4 & not_neo_c >=3 ~ "a",
    neo_c >= 5 & not_neo_c >=2 ~ "a",
    TRUE ~ "b"
  ))

sample_filter  <- neo_missense_summ %>% 
  filter(type=="a")
es_ccf <- readRDS("data/neo_nes_ccf06_1_remove_driver_samples_addp.rds")

es_ccf_filter <- es_ccf %>% filter(sample %in% sample_filter$sample)
median(es_ccf_filter$es)
median(es_ccf$es)

sim_all <- readRDS("../tmp/sim_2000_not_filter_all_ccf.rds")##this file is too large to push to Github
sim_filter <- sim_all %>% 
  filter(sample %in% es_ccf_filter$sample)
sim_filter %>%
  group_by(sim_num) %>%
  summarise(median_es=median(es)) -> summ
p <- WVPlots::ShadedDensity(frame = summ, 
                            xvar = "median_es",
                            threshold = median(es_ccf_filter$es),
                            title = "",
                            tail = "left")

p$layers[[1]]$aes_params$colour <- "red"
p$layers[[1]]$aes_params$size <- 1
p$layers[[2]]$aes_params$fill <- "blue"     #geom_ribbon
p$layers[[3]]$aes_params$colour <- "black" 
p$layers[[3]]$aes_params$size <- 1
#p$layers[[4]]$aes_params$label <- "Actual median ES" #geom_text
p1 <- p + labs(x="Simulation median es")+
  theme_prism()
p1


##exp
all_mut_exp <- readRDS("data/all_mut_tpm_not_filter.rds")
all_mut_exp <- all_mut_exp %>% 
  select(sample,tpm_exp,neo,chr,position,gene) %>% 
  rename(exp=tpm_exp)
es_exp <- readRDS("data/es_exp_filter_driver.rds")
exp_neo <- readRDS("../tmp/exp_neo_mut_counts.rds")
neo_mut_sig <- exp_neo %>% filter(p<0.05)
all_mut_exp_summ <- all_mut_exp %>% 
  group_by(sample) %>% 
  summarise(neo_c=sum(neo=="neo"),not_neo_c=sum(neo=="not_neo")) %>% 
  mutate(type=case_when(
    neo_c==1 & not_neo_c >=19 ~ "a",
    neo_c==2 & not_neo_c >=5 ~ "a",
    neo_c==3 & not_neo_c >=3 ~ "a",
    neo_c==4 & not_neo_c >=3 ~ "a",
    neo_c==5 & not_neo_c >=3 ~ "a",
    neo_c >= 6 & not_neo_c >=2 ~ "a",
    TRUE ~ "b"
  ))

sample_filter  <- all_mut_exp_summ %>% 
  filter(type=="a")

es_exp_filter <- es_exp %>% filter(sample %in% sample_filter$sample)
median(es_exp_filter$es)
median(es_exp$es)

sim_all <- readRDS("../tmp/sim_2000_not_filter_driver_exp_all.rds")##this file is too large to push to Github
sim_filter <- sim_all %>% 
  filter(sample %in% es_exp_filter$sample)
sim_filter %>%
  group_by(sim_num) %>%
  summarise(median_es=median(es)) -> summ
p <- WVPlots::ShadedDensity(frame = summ, 
                            xvar = "median_es",
                            threshold = median(es_ccf_filter$es),
                            title = "",
                            tail = "left")

p$layers[[1]]$aes_params$colour <- "red"
p$layers[[1]]$aes_params$size <- 1
p$layers[[2]]$aes_params$fill <- "blue"     #geom_ribbon
p$layers[[3]]$aes_params$colour <- "black" 
p$layers[[3]]$aes_params$size <- 1
#p$layers[[4]]$aes_params$label <- "Actual median ES" #geom_text
p1 <- p + labs(x="Simulation median es")+
  theme_prism()
p1

###突变数量分组 看免疫浸润
neo_missense <- readRDS("../tmp/mut_cal_es_filter_driver_samples.rds")
neo_missense %>% group_by(sample) %>% 
  summarise(neo_counts=sum(neo=="yes"),all_counts=n()) -> summ
ccf_immune <- left_join(
  ccf_immune,summ
)
quantile(ccf_immune$all_counts)
ccf_immune <- ccf_immune %>% 
  mutate(mut_group=case_when(
    all_counts <= 22 ~ "I",
    all_counts > 22 & all_counts <= 41 ~ "II",
    all_counts > 41 & all_counts <= 85 ~ "III",
    all_counts > 85  ~ "IV"
  ))
p1 <- ggplot(data=ccf_immune,aes(x=es_type,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8 T + NK")+
  facet_wrap(~mut_group)
p2 <- ggplot(data=ccf_immune,aes(x=es_type,y=iTreg+nTreg))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="iTreg + nTreg")+
  facet_wrap(~mut_group)
p3 <- ggplot(data=ccf_immune,aes(x=es_type,y=exp(CD8_T)/exp(iTreg+nTreg)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8 T / Treg")+
  facet_wrap(~mut_group)
p1 / p2 / p3

all_mut_exp <- readRDS("data/all_mut_tpm_not_filter.rds")
all_mut_exp %>% group_by(sample) %>% 
  summarise(neo_counts=sum(neo=="neo"),all_counts=n()) -> summ
exp_immune <- left_join(
  exp_immune,summ
)
quantile(exp_immune$all_counts)
exp_immune <- exp_immune %>% 
  mutate(mut_group=case_when(
    all_counts <= 21 ~ "I",
    all_counts > 21 & all_counts <= 41 ~ "II",
    all_counts > 41 & all_counts <= 92 ~ "III",
    all_counts > 92  ~ "IV"
  ))
p4 <- ggplot(data=exp_immune,aes(x=escape,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="CD8 T + NK")+
  facet_wrap(~mut_group)
p5 <- ggplot(data=exp_immune,aes(x=escape,y=iTreg+nTreg))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="iTreg + nTreg")+
  facet_wrap(~mut_group)
p6 <- ggplot(data=exp_immune,aes(x=escape,y=exp(CD8_T)/exp(iTreg+nTreg)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="CD8 T / Treg")+
  facet_wrap(~mut_group)
(p1 + p4) / (p2 + p5) / (p3 + p6)

###
all_samples <- data.table::fread("../tmp/samples_hla.txt",data.table = F,header = F,sep = " ")
com1 <- data.table::fread("../tmp/completed_samples",header = F,data.table = F) %>% 
  mutate(sample=gsub(".all_epitopes.tsv","",V1))
com2 <- data.table::fread("../tmp/completed_samples2",header = F,data.table = F) %>% 
  mutate(sample=gsub(".all_epitopes.tsv","",V1))
remain <- all_samples %>% 
  filter(!(V2 %in% c(com1$sample,com2$sample)))
remain <- remain %>% select(V1)
write.table(remain,file = "../tmp/remain_files",row.names = F,col.names = F,quote = F)

###exp mhcfurry
neo_nes <- readRDS("../tmp/es_exp_filter_driver_mhcfurry.rds")
sim_all <- readRDS("../tmp/sim_2000_filter_all_exp_mhcfurry.rds")
sim_all %>%
  group_by(cancer,sim_num) %>%
  summarise(median_es=median(es)) -> summ2
neo_nes$cancer <- get_cancer_type(neo_nes$sample)
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
p1 <- p + labs(x="Simulation median es")+
  theme_prism()
p1

###更新exp 癌症类型
res <- vector("list",30)
for (i in 1:30){
  dt <- sim_all %>% filter(cancer==neo_nes_summ$cancer[i])
  dt_summ <- dt %>%
    group_by(sim_num) %>%
    summarise(median_es=median(es))
  p <- WVPlots::ShadedDensity(frame = dt_summ,
                              xvar = "median_es",
                              threshold = neo_nes_summ$median_es[i],
                              title = neo_nes_summ$cancer[i],
                              tail = "left")
  p$layers[[1]]$aes_params$colour <- "red"
  p$layers[[1]]$aes_params$size <- 1
  p$layers[[2]]$aes_params$fill <- "blue"     #geom_ribbon
  p$layers[[3]]$aes_params$colour <- "black"
  p$layers[[3]]$aes_params$size <- 1
  p1 <- p + labs(x="Simulation median es")+
    theme_prism()
  res[[i]] <- p1
}
library(cowplot)
plot_grid(plotlist = res)
neo_nes_summ <- neo_nes_summ %>%
  rowwise() %>%
  mutate(p=mean(summ2$median_es[summ2$cancer==cancer] <= median_es))
saveRDS(neo_nes_summ,file = "data/neo_nes_summ_exp.rds")

