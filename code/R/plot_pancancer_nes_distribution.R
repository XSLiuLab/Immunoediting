### Plot pancancer nes distribution
library(tidyverse)
library(NeoEnrichment)

plot_pancancer <- function(mutation, neutral, text_position = 5.2){
  
  mutation <- readRDS(mutation)
  neutral <- readRDS(neutral)
  
  neutral <- neutral %>% 
    mutate(type = 'Neutral') %>% 
    select(nes, cancer_type, type)
  
  mutation <- mutation %>% 
    mutate(
      cancer_type = sapply(sample, NeoEnrichment::getcancer_type),
      type = 'Tumor'
    ) %>% 
    select(nes, cancer_type, type) 
  
  
  total <- rbind(mutation, neutral)
  
  ##p value
  ccf_p_n <- as.data.frame(table(mutation$cancer_type))
  ccf_p_n$p <- NA
  colnames(ccf_p_n) <- c("cancer_type","n","p")
  
  for(i in 1:length(ccf_p_n$cancer_type)){
    a = mutation %>%
      filter(cancer_type == ccf_p_n$cancer_type[i])
    cn <- neutral %>% 
      filter(cancer_type == ccf_p_n$cancer_type[i])
    p = wilcox.test(a$nes, mu = mean(cn$nes, na.rm = T), alternative = "less")$p.value
    ccf_p_n$p[i]=p
  }
  ccf_p_n <- ccf_p_n %>% 
    arrange(desc(cancer_type)) %>% 
    mutate(p = signif(p, 3),
           pn = paste0("P=",p,",n=",n),
           position = text_position,
           type = 'Tumor')
  mutation <- mutation %>% 
    group_by(cancer_type) %>% 
    mutate(
      order = median(nes)
    ) %>% 
    arrange(desc(order)) %>% 
    ungroup()
  cancer_list <- unique(mutation$cancer_type)
  
  total <- total %>%
    mutate(
      cancer_type = factor(cancer_type, levels = cancer_list)
    )
  
  ggplot(total, 
         mapping = aes(x = cancer_type, 
                       y = nes, 
                       fill = type)) + 
    geom_boxplot() +
    theme_classic()+
    geom_text(data = ccf_p_n, 
              mapping = aes(x = cancer_type,
                            y = position,
                            label = pn))+
    theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "transparent",colour = NA)
    ) + 
    labs(x = NULL) + 
    guides(fill = F) + 
    coord_flip() 
  
}