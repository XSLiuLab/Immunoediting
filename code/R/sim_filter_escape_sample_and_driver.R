library(dplyr)
library(parallel)
library(NeoEnrichment)
neo_nes <- readRDS("data/filter_escape.rds")
all_mut_ccf <- readRDS("data/all_mut_ccf_tpm.rds")
all_mut_ccf <- all_mut_ccf %>%
  rename(ccf=ccf_hat) %>%
  mutate(neo=ifelse(neo=="neo","yes","no"))
neo_missense <- all_mut_ccf %>% filter(sample %in% neo_nes$sample)
neo_missense <- neo_missense %>% select(sample,neo,ccf) %>% filter(!is.na(ccf))

cal_nes_warp <- function(dt){
  results_ccf <- vector("list",length = length(unique(dt$sample)))
  names(results_ccf) <- unique(dt$sample)
  
  cl <- makeCluster(getOption("cl.cores",32),type="FORK")
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
saveRDS(res,file = "data/sim_filter_excapeAnddriver.rds")
