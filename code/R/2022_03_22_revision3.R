files <- list.files("../../evolution/immunetherapy/")
files <- files[grepl("neoantigens.unfiltered.txt",files)]
re <- vector("list",25)
for (i in 1:25){
  mut <- read_table(paste("../../evolution/immunetherapy/",files[i],sep = ""),col_names = paste0(rep("V",27),c(1:27)))
  
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
    filter(IC50<50 & bindlevel=="SB" & novelty==1) %>%
    distinct(index,.keep_all = T)
  
  mut <- mut %>%
    distinct(index,.keep_all = T) %>%
    mutate(neo = ifelse(index %in% neo$index , "neo","not_neo"))
  mut <- mut %>% 
    select(sample,chr,position,gene,exp,neo) %>% 
    mutate(gene=gsub("\\:.+","",gene))
  mut$sample <- paste(gsub("_.+","",files[i]),mut$sample,sep = "_")
  mut <- mut %>% filter(!is.na(exp))
  re[[i]] <- mut
}

all_mut <- bind_rows(re)
saveRDS(all_mut,file = "data/Immunotherapy/all_mut_exp.rds")

##在服务器上计算了esrna
all_clinical <- readRDS("data/Immunotherapy/all_clinical.rds")

neo_nes <- readRDS("data/Immunotherapy/nes_immunetherapy_exp.rds")
neo_nes <- left_join(neo_nes,all_clinical,by="sample")
neo_nes <- neo_nes %>% filter(!is.na(response2))

##cox analysis
p1 <- show_forest(neo_nes,covariates = "es",time = "OS.time",status = "OS")
p1

neo_nes <- neo_nes %>%
  mutate(obj=ifelse(response2=="response",1,0))
model <- glm(obj ~ es, data = neo_nes, family = "binomial")
summary(model)
performance::performance_hosmer(model)

ggplot(neo_nes, aes(x=es, y=obj)) +
  geom_point(alpha=.5) +
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE,
              method.args = list(family=binomial),color="red")+
  theme_bw()+
  theme(axis.title.y = element_blank())+
  labs(x="ES",title = "H-L test P value = 0.085")
###逃逸 生存
pancancer_survial <- readRDS("data/pancancer_survial.rds")
all_escape <- readRDS("data/escape_all.rds")
all_escape <- left_join(all_escape,pancancer_survial)
all_escape$apm_mut <- as.character(all_escape$apm_mut)
show_forest(all_escape,covariates = "apm_mut",time = "OS.time",status = "OS")
show_forest(all_escape,covariates = "over_exp",time = "OS.time",status = "OS")
show_forest(all_escape,covariates = "LOH_status",time = "OS.time",status = "OS")
show_forest(all_escape,covariates = "es_down",time = "OS.time",status = "OS")

all_escape <- all_escape %>% 
  mutate(apm_mut = ifelse(apm_mut == "TRUE","yes","no")) %>% 
  rename(`APM mutation`=apm_mut,
         `Checkpoint over-expression` = over_exp,
         `HLA LOH` = LOH_status,
         `Rna downregulation` = es_down)
show_forest(all_escape,covariates = c("APM mutation","Checkpoint over-expression",
                                      "HLA LOH","Rna downregulation"),
            time = "OS.time",status = "OS",merge_models = T)

###结合
neo_nes_exp <- readRDS("data/Immunotherapy/nes_immunetherapy_exp.rds")
neo_nes_ccf <- readRDS("data/Immunotherapy/neo_nes_ccf_addp.rds")
neo_nes_ccf <- left_join(
  neo_nes_ccf %>% rename(es_ccf=es),
  neo_nes_exp %>% rename(es_exp=es)
)
all_clinical <- readRDS("data/Immunotherapy/all_clinical.rds")
neo_nes_ccf <- left_join(neo_nes_ccf,all_clinical,by="sample")

neo_nes_ccf <- neo_nes_ccf %>% 
  mutate(ccf_type=ifelse(es_ccf<0 & p<0.05,"yes","no"),
         exp_type=ifelse(es_exp<0 & p_value<0.05,"yes","no")) %>% 
  mutate(type=case_when(
    ccf_type == "yes" & exp_type == "yes" ~ "I",
    ccf_type == "yes" & exp_type == "no" ~ "II",
    ccf_type == "no" & exp_type == "yes" ~ "III",
    ccf_type == "no" & exp_type == "no" ~ "IV"
  ))
table(neo_nes_ccf$type)
neo_nes_ccf <- neo_nes_ccf %>% filter(!is.na(response2))
show_forest(neo_nes_ccf,covariates = "type",time = "OS.time",status = "OS")
