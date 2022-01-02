library(dplyr)
library(parallel)
library(NeoEnrichment)

neo_nes <- readRDS("~/Immunoediting/data/neo_nes_ccf06_1_harmonic1.rds")
all_mut_ccf <- readRDS("~/Immunoediting/data/all_mut_ccf_tpm_harmonic1.rds")
all_mut_ccf <- all_mut_ccf %>%
  dplyr::rename(ccf=ccf_hat) %>%
  mutate(neo=ifelse(neo2=="neo","yes","no"))
#mean(all_mut_ccf$neo=="yes")
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
