all_clinical <- readRDS("data/Immunotherapy/all_clinical.rds")
neo_nes <- readRDS("data/Immunotherapy/nes_immunetherapy.rds")
all_mut_ici <- readRDS("../tmp/all_mut_ici_withrefalt.rds")
neo_nes <- left_join(neo_nes,all_clinical,by="sample")
neo_nes <- neo_nes %>% filter(!is.na(response2))
all_ccf <- readRDS("data/Immunotherapy/all_ccf.rds") %>% 
  filter(!is.na(cancer_cell_frac))
#clonal TMB
all_mut_ccf_ici <- readRDS("data/Immunotherapy/all_mut_ccf_ici.rds")
all_mut_ccf_ici %>%
  group_by(sample) %>%
  summarise(all_tmb=n()/38,clonal_tmb=sum(cancer_cell_frac>=0.9)/38) -> tmb

neo_nes <- left_join(neo_nes,tmb)

##indel tmb
all_indel_counts <- readRDS("data/Immunotherapy/all_indel_counts.rds")
neo_nes <- left_join(neo_nes,all_indel_counts)

##escape NMD indel TMB
nes_ici_escapeNMD_indel <- readRDS("data/Immunotherapy/nes_ici_escapeNMD_indel.rds")
neo_nes <- left_join(neo_nes,nes_ici_escapeNMD_indel %>% select(-es))

##signature
#smoke sig4
sigs <- readRDS("data/sig_ici_need.rds")
sigs <- sigs %>% 
  dplyr::rename(smoke=Signature.4,UV=Signature.7) %>%
  mutate(APOBEC=Signature.2+Signature.13) %>% 
  dplyr::select(-Signature.13,-Signature.2)
neo_nes <- left_join(
  neo_nes,
  sigs
)

###SERPINB3 mut
SERPINB3_mut <- all_mut_ici %>%
  group_by(sample) %>% 
  summarise(SERPINB3_mut=ifelse("SERPINB3" %in% gene,"yes","no"))
neo_nes <- left_join(
  neo_nes,
  SERPINB3_mut
)

##gene exp
need_gene_exp_ici <- readRDS("data/Immunotherapy/need_gene_exp_ici.rds")
neo_nes <- left_join(
  neo_nes,
  need_gene_exp_ici
)

###gene exp signature
#Ayers IMPRES POPLAR CYT
ayer_score_ici <- readRDS("data/Immunotherapy/ayer_score_ici.rds")
IMPRES_ici <- readRDS("data/Immunotherapy/IMPRES_ici.rds")
POPLAR_score_ici <- readRDS("data/Immunotherapy/POPLAR_score_ici.rds")
cyt_ici <- readRDS("data/Immunotherapy/cyt_ici.rds")

neo_nes <- left_join(
  neo_nes,ayer_score_ici
) %>% left_join(.,IMPRES_ici) %>% left_join(.,POPLAR_score_ici) %>% 
  left_join(.,cyt_ici)
colnames(neo_nes)

saveRDS(neo_nes,file = "data/Immunotherapy/all_factors_ici.rds")

##比较
neo_nes <- readRDS("data/Immunotherapy/all_factors_ici.rds")
library(ezcox)
show_forest(neo_nes,covariates = c("es",colnames(neo_nes)[7:10],
                                   colnames(neo_nes)[14:24],"UV"),
            time = "OS.time",status = "OS",merge_models = T)
show_forest(neo_nes,covariates = c("es",colnames(neo_nes)[7:24]),
            time = "OS.time",status = "OS",merge_models = T)
#show_forest(neo_nes,covariates = c("UV","APOBEC"),
#            time = "OS.time",status = "OS",merge_models = T)

#############
neo_nes <- neo_nes %>%
  mutate(obj=ifelse(response2=="response",1,0))
res_plot <- vector("list",14)
names(res_plot) <- c("es",colnames(neo_nes)[7:10],
                     colnames(neo_nes)[15:20],"UV","POPLAR_score","cyt")
for (i in 1:14){
  dt <- neo_nes %>% select(obj,names(res_plot)[i])
  colnames(dt)[2] <- "tt"
  model2 <- glm(obj ~ tt, 
                data = dt, family = "binomial")
  a <- performance::performance_hosmer(model2)
  p <- ggplot(dt, aes(x=tt, y=obj)) +
    geom_point(alpha=.5) +
    stat_smooth(method="glm", se=FALSE, fullrange=TRUE,
                method.args = list(family=binomial),color="red")+
    theme_bw()+
    theme(axis.title.y = element_blank())+
    labs(x=names(res_plot)[i],title = paste0("H-L test P value = ",round(a$p.value,digits = 3)))
  res_plot[[i]] <- p
}
library(cowplot)
plot_grid(plotlist = res_plot)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggprism)
ggplot(data=neo_nes,aes(x=response2,y=all_tmb))+
  geom_boxplot()+
  stat_compare_means()

model2 <- glm(obj ~ es, data = neo_nes, family = "binomial")
dt <- neo_nes %>% select(es)
predict(model2,newdata = dt)

aa <- c("es",colnames(neo_nes)[7:10],
        colnames(neo_nes)[15:19], colnames(neo_nes)[21:24],"UV")
res <- vector("list",15)
for (i in 1:length(aa)){
  dt <- neo_nes %>% 
    select(response2,aa[i])
  colnames(dt)[2] <- "tt"
  p <- ggplot(data=dt,aes(x=response2,y=tt))+
    geom_boxplot()+
    stat_compare_means()+
    labs(y=aa[i])+
    theme_prism()+
    theme(axis.title.x = element_blank())
  res[[i]] <- p
}
library(cowplot)
plot_grid(plotlist = res,nrow = 4,ncol = 4)


model2 <- glm(obj ~ CD8A, 
              data = neo_nes, family = "binomial")
summary(model2)
exp(coef(model2))

#########
neo_nes <- neo_nes %>% 
  mutate(cohort=gsub("_.+","",sample))
liu <- neo_nes %>% filter(cohort=="liu")
nadeem <- neo_nes %>% filter(cohort=="nadeem")
willy <- neo_nes %>% filter(cohort=="willy")

model2 <- glm(obj ~ CXCL9, 
              data = neo_nes, family = "binomial")
summary(model2)
exp(coef(model2))
round(exp(cbind(coef(model2), confint(model2))), digits=3) %>% 
  as.data.frame()

