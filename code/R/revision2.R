library(dplyr)
library(ggplot2)
library(ggprism)
library(NeoEnrichment)
library(ggpubr)
library(patchwork)
# Rank2 -------------------------------------------------------------------
neo_nes <- readRDS("data/neo_es_ccf_IC50_Rank2.rds")
sim_all <- readRDS("../tmp/sim_2000_filter_driver_IC50_Rank2_all_ccf.rds")
neo_nes$cancer <- NeoEnrichment::get_cancer_type(neo_nes$sample)
neo_nes_summ <- neo_nes %>% 
  group_by(cancer) %>% summarise(median_es=median(es))

#The following code can be used to plot cancer type simulation which showd in FigS4
#res <- vector("list",30)
# for (i in 1:30){
#   dt <- sim_all %>% filter(cancer==neo_nes_summ$cancer[i])
#   dt_summ <- dt %>% 
#     group_by(sim_num) %>%
#     summarise(median_es=median(es))
#   p <- WVPlots::ShadedDensity(frame = dt_summ, 
#                               xvar = "median_es",
#                               threshold = neo_nes_summ$median_es[i],
#                               title = neo_nes_summ$cancer[i],
#                               tail = "left")
#   p$layers[[1]]$aes_params$colour <- "red"
#   p$layers[[1]]$aes_params$size <- 1
#   p$layers[[2]]$aes_params$fill <- "blue"     #geom_ribbon
#   p$layers[[3]]$aes_params$colour <- "black" 
#   p$layers[[3]]$aes_params$size <- 1
#   p1 <- p + labs(x="Simulation median es")+
#     theme_prism()
#   res[[i]] <- p1
# }
# 
# plot_grid(plotlist = res)

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
  theme_prism()
p1

##rna
neo_nes <- readRDS("data/nes_exp_Rank2_IC50.rds")
sim_all <- readRDS("../tmp/sim_2000_es_exp_Rank2_IC50_all.rds")
neo_nes$cancer <- NeoEnrichment::get_cancer_type(neo_nes$sample)
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
p2 <- p + labs(x="Simulation median es")+
  theme_prism()
p2
p1/p2

# 考虑HLA-LOH ---------------------------------------------------------------
library(dplyr)
nes_ccf <- readRDS("data/neo_nes_ccf06_1_remove_driver_samples_addp.rds")
nes_exp <- readRDS("data/es_exp_filter_driver.rds")
checkpoint_over <- readRDS("data/checkpoint_over.rds")
apm_mut_sample <- readRDS("data/apm_mut_sample.rds")
loh <- readxl::read_xlsx("data/TCGA_HLA_benchmark_20200810.xlsx",sheet = 2)
loh <- loh %>%
  filter(LOH_cal == TRUE) %>% 
  mutate(LOH_status = ifelse(nchar(LossAllele)==1,"no","yes"))
table(loh$LOH_status)
loh <- loh %>% select(Sample_Barcode,LOH_status)
saveRDS(loh,file = "data/HLA_loh.rds")
LOH <- readRDS("data/HLA_loh.rds")

apm_mut_sample <- apm_mut_sample %>% 
  mutate(sample=substr(sample,1,12))
checkpoint_over <- checkpoint_over %>% 
  mutate(sample=substr(sample,1,12)) %>% 
  mutate(over_exp=ifelse(PDL1_over == "yes" | CTLA4_over == "yes" ,"yes","no")) %>% 
  select(sample,cancer,over_exp) %>% 
  distinct_all(.keep_all = T)
checkpoint_over <- checkpoint_over[!duplicated(checkpoint_over$sample),]
nes_exp <- nes_exp %>%
  mutate(sample=substr(sample,1,12))
LOH <- LOH %>% rename(sample=Sample_Barcode)
all_escape <- Reduce(left_join,list(apm_mut_sample,checkpoint_over,nes_exp,LOH))
all_escape <- na.omit(all_escape)
all_escape <- all_escape %>% 
  mutate(es_down=ifelse(es<0 & p_value<0.05,"yes","no"))
saveRDS(all_escape,file = "data/escape_all.rds")
dt <- all_escape %>% 
  group_by(cancer) %>% 
  summarise(`APM mutation` = mean(apm_mut == TRUE),
         `Checkpoint over-expression` = mean(over_exp=="yes"),
         `Rna downregulation` = mean(es_down=="yes"),
         `LOH`=mean(LOH_status=="yes")) %>% as.data.frame()
rownames(dt) <- dt$cancer
dt <- dt %>% select(-cancer)

library(ComplexHeatmap)
library(circlize)
dt <- dt * 100
col_fun = colorRamp2(c(0, 50,100), c("green", "white", "red"))
Heatmap(dt, name = "Sample proportion (%)",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", dt[i, j]), x, y, gp = gpar(fontsize = 10))
        },cluster_rows = F,cluster_columns = F)
###########

# 探索免疫浸润 ------------------------------------------------------------------
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

ggplot(data=ccf_immune,aes(x=es_type,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8 T + NK")+
  facet_wrap(~cancer)

ggplot(data=ccf_immune,aes(x=test,y=CD8_T+NK))+
  geom_point()+
  geom_smooth(method = lm)+
  stat_cor()+
  theme_prism()+
  labs(x="ES CCF",y="CD8 T + NK")+
  facet_wrap(~cancer)

ccf_immune_escape <- left_join(
  ccf_immune %>% mutate(sample=substr(sample,1,12)),
  escape_all %>% 
    mutate(escape=ifelse(apm_mut==TRUE | over_exp=="yes" | LOH_status=="yes" | es_down=="yes","yes","no")) %>% 
    select(sample,escape)
) %>% 
  filter(!is.na(escape))

exp_immune <- left_join(
  escape_all,
  pancancer_subtcells %>% mutate(sample=substr(sample,1,12))
) %>% distinct(sample,.keep_all = T) %>% 
  mutate(escape=ifelse(apm_mut==TRUE | over_exp=="yes" | LOH_status=="yes" | es_down=="yes","yes","no")) 

escape_samples <- ccf_immune_escape %>% 
  filter(escape=="no")
ggplot(data=escape_samples,aes(x=es_type,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8 T + NK")

msi <- readxl::read_xlsx("../tmp/1-s2.0-S1672022920300218-mmc21.xlsx",skip = 2)

ccf_immune_msi <- left_join(
  ccf_immune,
  msi %>% 
    select(`Tumor sample used in this study`,`MSI status (determined by MSI-PCR)`) %>% 
    rename(sample=`Tumor sample used in this study`,MSI=`MSI status (determined by MSI-PCR)`) %>% 
    mutate(sample=substr(sample,1,16))
) %>% filter(!is.na(MSI))

exp_immune_msi <- left_join(
  exp_immune,
  msi %>% 
    select(`Tumor sample used in this study`,`MSI status (determined by MSI-PCR)`) %>% 
    rename(sample=`Tumor sample used in this study`,MSI=`MSI status (determined by MSI-PCR)`) %>% 
    mutate(sample=substr(sample,1,12))
) %>% filter(!is.na(MSI)) %>% filter(!is.na(CD8_T))

ggplot(data=ccf_immune_msi,aes(x=es_type,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Elimination",y="CD8 T + NK")+
  facet_wrap(~MSI)

ggplot(data=ccf_immune_msi,aes(x=es,y=CD8_T+NK))+
  geom_point()+
  geom_smooth(method = lm)+
  stat_cor()+
  theme_prism()+
  labs(x="ES CCF",y="CD8 T + NK")+
  facet_wrap(~MSI)

ggplot(data=exp_immune_msi,aes(x=escape,y=CD8_T+NK))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="CD8 T + NK")+
  facet_wrap(~MSI)
ggplot(data=exp_immune_msi,aes(x=escape,y=iTreg+nTreg))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="iTreg + nTreg")+
  facet_wrap(~MSI)
ggplot(data=exp_immune_msi,aes(x=escape,y=(CD8_T+NK)/(iTreg+nTreg)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="Escape",y="CD8_T+NK/iTreg+nTreg")+
  geom_hline(yintercept=1)+
  facet_wrap(~MSI)

p1 <- ggplot(data=exp_immune_msi,aes(x=es,y=iTreg+nTreg))+
  geom_point()+
  geom_smooth(method = lm)+
  stat_cor()+
  theme_prism()+
  labs(x="ES EXP",y="iTreg+nTreg")
p2 <- ggplot(data=exp_immune_msi,aes(x=es,y=CD8_T+NK))+
  geom_point()+
  geom_smooth(method = lm)+
  stat_cor()+
  theme_prism()+
  labs(x="ES EXP",y="CD8_T+NK")
p3 <- ggplot(data=exp_immune_msi,aes(x=es,y=iTreg+nTreg))+
  geom_point()+
  geom_smooth(method = lm)+
  stat_cor()+
  theme_prism()+
  labs(x="ES EXP",y="iTreg+nTreg")+
  facet_wrap(~MSI)
p4 <- ggplot(data=exp_immune_msi,aes(x=es,y=CD8_T+NK))+
  geom_point()+
  geom_smooth(method = lm)+
  stat_cor()+
  theme_prism()+
  labs(x="ES EXP",y="CD8_T+NK")+
  facet_wrap(~MSI)
(p1 + p2) / (p3 + p4)

ggplot(data=exp_immune_msi,aes(x=es,y=(CD8_T+NK)/(iTreg+nTreg)))+
  geom_point()+
  geom_smooth(method = lm)+
  stat_cor()+
  theme_prism()+
  labs(x="ES EXP",y="CD8 T + NK/ iTreg + nTreg")+
  facet_wrap(~MSI)

##批量画散点图
ccf_immune <- ccf_immune %>% 
  mutate(cd8_nk = CD8_T + NK,Treg=iTreg+nTreg)
cells <- colnames(ccf_immune)[5:28]
plot_list <- vector("list",24)
for (i in seq_along(plot_list)){
  dt <- ccf_immune %>% 
    select(es,cells[i])
  colnames(dt)[2] <- "immune"
  p <- ggplot(data=dt,aes(x=es,y=immune))+
    geom_point()+
    labs(y=cells[i])
  plot_list[[i]] <- p
}
library(cowplot)
plot_grid(plotlist = plot_list)

##cancer_type 
cancers <- unique(ccf_immune$cancer)
plot_list <- vector("list",30)
for (i in seq_along(plot_list)){
  dt <- ccf_immune %>%
    filter(cancer==cancers[i]) %>% 
    select(es,cd8_nk)
  colnames(dt)[2] <- "immune"
  p <- ggplot(data=dt,aes(x=es,y=immune))+
    geom_point()+
    labs(y="CD8 + NK",title=cancers[i])
  plot_list[[i]] <- p
}
plot_grid(plotlist = plot_list)


ggplot(data=ccf_immune_escape,aes(x=es_type,y=cd8_nk))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()+
  labs(x="ES CCF",y="CD8 + NK")+
  facet_wrap(~escape)

sig <- es_ccf %>% filter(p<0.05)
ggplot(data=es_ccf,aes(x=es,y=p))+
  geom_point()+
  geom_point(data=sig,aes(x=es,y=p),color="red")


neo_missense <- readRDS("../tmp/mut_cal_es_filter_driver_samples.rds")
neo_missense %>% group_by(sample) %>% 
  summarise(neo_counts=sum(neo=="yes"),all_counts=n()) -> summ
ccf_immune <- left_join(
  ccf_immune,summ
)
ggplot(data=ccf_immune,aes(x=log(all_counts+1),y=CD8_T+NK))+
  geom_point()

###按照mutation burden 分组
es_ccf <- left_join(
  es_ccf,summ
)
quantile(es_ccf$ratio)
es_ccf <- es_ccf %>% 
  mutate(mut_group=case_when(
    ratio <=0.061538462 ~ "I",
    ratio > 0.061538462 & ratio <= 0.093750000 ~ "II",
    ratio > 0.093750000 & ratio <= 0.134723423 ~ "III",
    ratio > 0.134723423  ~ "IV"
  ))
es_ccf %>% 
  group_by(mut_group) %>% 
  summarise(e=median(es))
 
sim_all <- readRDS("../tmp/sim_2000_not_filter_all_ccf.rds")##this file is too large to push to Github
get_plot <- function(x){
  dt <- es_ccf %>% filter(mut_group==x)
  sim_filter <- sim_all %>% 
    filter(sample %in% dt$sample)
  sim_filter %>%
    group_by(sim_num) %>%
    summarise(median_es=median(es)) -> summ
  p <- WVPlots::ShadedDensity(frame = summ, 
                              xvar = "median_es",
                              threshold = median(dt$es),
                              title = "",
                              tail = "left")
  
  p$layers[[1]]$aes_params$colour <- "red"
  p$layers[[1]]$aes_params$size <- 1
  p$layers[[2]]$aes_params$fill <- "blue"     #geom_ribbon
  p$layers[[3]]$aes_params$colour <- "black" 
  p$layers[[3]]$aes_params$size <- 1
  p1 <- p + labs(x="Simulation median es")+
    theme_prism()+
    labs(title = x)
  cat(unique(dt$sample) %>% length())
  return(p1)
}

p1 <- get_plot("I")
p2 <- get_plot("II")
p3 <- get_plot("III")
p4 <- get_plot("IV")
p1 + p2 + p3 + p4

ggplot(data=es_ccf,aes(x=log(all_counts+1),y=p))+
  geom_point()

ccf_escape <- left_join(
  es_ccf %>% mutate(sample=substr(sample,1,12)),
  escape_all %>% 
    mutate(escape=ifelse(apm_mut==TRUE | over_exp=="yes" | LOH_status=="yes" | es_down=="yes","yes","no")) %>% 
    select(sample,escape)
) %>% 
  filter(!is.na(escape))
ggplot(data=ccf_escape,aes(x=escape,y=log(all_counts+1)))+
  geom_boxplot()+
  stat_compare_means()+
  theme_prism()
no_escape <- ccf_escape %>% filter(escape=="no")
quantile(no_escape$all_counts)

dt <- ccf_escape %>% filter(escape=="yes")
sim_filter <- sim_all %>% 
  filter(substr(sample,1,12) %in% dt$sample)
sim_filter %>%
  group_by(sim_num) %>%
  summarise(median_es=median(es)) -> summ
p <- WVPlots::ShadedDensity(frame = summ, 
                            xvar = "median_es",
                            threshold = median(dt$es),
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

ggplot(data=ccf_immune,aes(x=es,y=(CD8_T+NK)))+
  geom_point()+
  geom_smooth(method = lm)+
  stat_cor()+
  theme_prism()+
  labs(x="ES CCF",y="CD8 T + NK")+
  facet_wrap(~mut_group)

pcor.test(ccf_immune$es,ccf_immune$CD8_T,ccf_immune$neo_counts)







