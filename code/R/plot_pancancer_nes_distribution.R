### Plot pancancer nes distribution
library(tidyverse)
library(NeoEnrichment)

plot_pancancer <- function(mutation, neutral, text_position = 5.2,size){
  
  neutral <- neutral %>% 
    mutate(type = 'Neutral') %>% 
    select(nes, cancer, type)
  
  mutation <- mutation %>% 
    mutate(
      type = 'Tumor'
    ) %>% 
    select(nes, cancer, type) 
  
  
  total <- rbind(mutation, neutral)
  
  ##p value
  ccf_p_n <- as.data.frame(table(mutation$cancer))
  ccf_p_n$p <- NA
  colnames(ccf_p_n) <- c("cancer_type","n","p")
  
  for(i in 1:length(ccf_p_n$cancer_type)){
    a = mutation %>%
      filter(cancer == ccf_p_n$cancer_type[i])
    cn <- neutral %>% 
      filter(cancer == ccf_p_n$cancer_type[i])
    p = wilcox.test(a$nes, cn$nes, alternative = "less")$p.value
    ccf_p_n$p[i]=p
  }
  ccf_p_n <- ccf_p_n %>% 
    arrange(desc(cancer_type)) %>% 
    mutate(p = formatC(p, format = "e", digits = 2),
           pn = paste0("P=",p,",n=",n),
           position = text_position,
           type = 'Tumor')
  mutation <- mutation %>% 
    group_by(cancer) %>% 
    mutate(
      order = median(nes)
    ) %>% 
    arrange(desc(order)) %>% 
    ungroup()
  cancer_list <- unique(mutation$cancer)
  
  total$cancer <- factor(total$cancer,levels = cancer_list)
  
  ggplot(total, 
         mapping = aes(x = cancer, 
                       y = nes, 
                       fill = type)) + 
    scale_fill_manual(values=c("white", "#FC8D62"))+
    geom_boxplot(outlier.size = 1,outlier.alpha = 0.3) +
    theme_classic()+
    geom_text(data = ccf_p_n, 
              mapping = aes(x = cancer_type,
                            y = position,
                            label = pn),size=size)+
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
