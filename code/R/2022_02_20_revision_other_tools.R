files <- list.files("/public/slst/home/wutao2/TCGA_pvacseq/out_reorganized",full.names = T)
re <- vector("list",8044)

for (i in 1:8044){
  mut <- data.table::fread(files[i],data.table = F)
  mut <- mut %>% mutate(
    index = paste(V1,V2,V3,V4,V5,sep = ":")
  )
  neo_mhcfurry <- mut %>%
    filter(V8<50 & V9 <0.5) %>%
    distinct(index,.keep_all = T)
  neo_mhcnuggets <- mut %>%
    filter(V10 < 50) %>%
    distinct(index,.keep_all = T)
  mut <- mut %>%
    distinct(index,.keep_all = T) %>%
    mutate(neo_MHCfurry = ifelse(index %in% neo_mhcfurry$index , "neo","not_neo")) %>% 
    mutate(neo_MHCnuggets = ifelse(index %in% neo_mhcnuggets$index , "neo","not_neo"))
  mut <- mut %>% select(V6,V11,index,neo_MHCfurry,neo_MHCnuggets)
  re[[i]] <- mut
}

all_mut <- bind_rows(re)
saveRDS(all_mut,file = "/public/slst/home/wutao2/TCGA_pvacseq/res/all_mut_neo.rds")