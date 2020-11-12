library(dplyr)
cales <- function(test,neoantigen_list){
  a <- sum(test[test$mut %in% neoantigen_list,"rank"])
  b <- (nrow(test)-length(neoantigen_list))
  
  test$re <- ifelse(test$mut %in% neoantigen_list,((test$rank)/a),-(1/b))
  test$cum_re <- cumsum(test$re)
  es <- max(0,test$cum_re)-abs(min(test$cum_re,0))
}

cal_p_and_normalized <- function(es,neo_list,test){
  sample_res <- data.frame(res=rep(1,1000))
  for (i in 1:1000) {
    neoantigen_list <- sample(test$mut,length(neo_list),replace = F)
    sample_res$res[i] <- cales(test,neoantigen_list)
  }
  p <- ifelse(es<0,mean(sample_res$res<es),mean(sample_res$res>es))
  nes <- ifelse(es<0,es/abs(mean(sample_res$res[sample_res$res<0])),
                es/mean(sample_res$res[sample_res$res>0]))
  es_p <- list(es=es,nes=nes,p=p,sample_res=sample_res)
  return(es_p)
}


cales_t <- function(file,calp=FALSE){
  a <- (nrow(file)/2)
  test <- file %>%
    dplyr::arrange(desc(ccf)) %>% 
    dplyr::mutate(rank = rank(ccf)) %>%
    dplyr::mutate(rank=abs(a-rank)+1)
  
  neo_list <- test[test$mut_neo=="true","mut"]
  if(length(neo_list)==0){
    r <- list(es="no",nes="no",p="no",neo_list=length(neo_list))
  }else{
    b <- length(neo_list)
    es <- cales(test,neo_list)
    if(calp==T){
      r <- cal_p_and_normalized(es,neo_list,test)
      r$neo_list <- b   
    }else{r <- es}
  }
  return(r)
}
######
dt <- data.frame(id=1:200,es=NA,nes=NA,p=NA,neo_counts=NA,mt_counts=NA,neo_list=NA)

#set.seed(<seed>)
for (i in 1:200){
  t1 <- read.table(paste("/public/slst/home/wutao2/julia_simulation/out/minus_<dir>/mutsumm_",i,".txt",sep = ""),sep = ",",header = T)
  vaf1 <- read.table(paste("/public/slst/home/wutao2/julia_simulation/out/minus_<dir>/vaf_preIT_",i,".txt",sep = ""))
  mt_counts <- vaf1 %>% filter(V2>5) %>% nrow()
  vaf1 <- vaf1 %>% 
    rename(mut=V1,cells=V2) %>% 
    filter(cells>5) %>% 
    mutate(ccf=cells/1e5) %>% 
    left_join(.,t1 %>% distinct(mut,.keep_all=T),
              by="mut")
  neo_counts <- sum(vaf1$mut_neo=="true")
  a <- cales_t(file = vaf1,calp = T)
  dt[i,]$es <- a$es
  dt[i,]$nes <- a$nes
  dt[i,]$p <- a$p
  dt[i,]$neo_list <- a$neo_list
  dt[i,]$neo_counts <- neo_counts
  dt[i,]$mt_counts <- mt_counts
}

saveRDS(dt,file="/public/slst/home/wutao2/cal_es_simulation/out/dt_<dir>.rds")