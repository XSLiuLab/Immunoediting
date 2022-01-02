library(GenomicDataCommons)

all_mut_ccf_tpm <- readRDS("~/Immunoediting/data/all_mut_ccf_tpm.rds")
all_mut_ccf_tpm$cancer <- EasyBioinfo::get_cancer_type(all_mut_ccf_tpm$sample,parallel = T,cores = 10)
need_cancer <- unique(all_mut_ccf_tpm$cancer)

res <- vector("list",30)
for (i in 1:30){
  q = files() %>% 
    GenomicDataCommons::filter(~ cases.project.project_id == paste0("TCGA-",need_cancer[i]) 
                               & data_type == 'Aligned Reads' 
                               & experimental_strategy == 'WXS' 
                               & data_format == 'BAM') %>% 
    GenomicDataCommons::select('file_id') %>% 
    expand('analysis.metadata.read_groups') 
  
  file_ids = ids(q) 
  z = results_all(q)
  read_length_list = sapply(z$analysis$metadata$read_groups,'[[','read_length')
  
  library(dplyr)
  t <- z$analysis$metadata$read_groups
  names(t) <- z$file_id
  a <- t %>% bind_rows(.id = "file_ids") %>% as_tibble()
  res[[i]] <- a
}

re <- bind_rows(res)
saveRDS(re,file = "~/data/TCGA_read_lengths.rds")

re <- re %>% 
  dplyr::select(experiment_name,read_length) %>% 
  mutate(sample=substr(experiment_name,1,16)) %>% 
  distinct(sample,read_length,.keep_all = T) %>% 
  dplyr::filter(grepl("TCGA",sample)) %>% 
  dplyr::filter(as.numeric(substr(sample,14,15)) <= 9)
re$cancer <- EasyBioinfo::get_cancer_type(re$sample)
saveRDS(re,file = "~/data/TCGA_read_lengths_re.rds")
###
#get BAM file manifest
manifest = GenomicDataCommons::files() %>%   
  GenomicDataCommons::filter(~ experimental_strategy == "WXS" &
                               data_format == "BAM") %>%   
  GenomicDataCommons::manifest()

dt <- manifest %>% 
  mutate(tmp=gsub(".+[.]TCGA","TCGA",filename)) %>% 
  mutate(sample=substr(tmp,1,16)) %>% 
  dplyr::filter(!grepl("hg19",filename)) %>% 
  dplyr::filter(sample %in% re$sample)
#get BAM BAI file manifest
res <- vector("list",length(unique(dt$id)))
for (i in 5458:length(res)){
  con = curl::curl(paste0("https://api.gdc.cancer.gov/files/", dt$id[i], "?pretty=true&expand=index_files"))
  tbl = jsonlite::fromJSON(con)
  bai = data.frame(id = tbl$data$index_files$file_id,
                   filename = tbl$data$index_files$file_name,
                   md5 = tbl$data$index_files$md5sum,
                   size = tbl$data$index_files$file_size,
                   state = tbl$data$index_files$state)
  res[[i]] <- bai
  cat("complete",i,"\n")
}

re <- bind_rows(res)
write.table(re,file = "~/data/TCGA_bai_manifest.txt",sep = "\t",row.names = F,quote = F)
saveRDS(re,file="~/data/TCGA_bai_manifest.rds")

###计算
TCGA_read_lengths_re <- readRDS("~/data/TCGA_read_lengths_re.rds")
which(duplicated(TCGA_read_lengths_re$sample))
test <- TCGA_read_lengths_re %>% filter(sample %in% TCGA_read_lengths_re$sample[which(duplicated(TCGA_read_lengths_re$sample))])

files <- list.files("~/data/TCGA_bai/stat_files/")
dt <- data.frame(file=files)
dt$samples <- gsub(".+[.]TCGA","TCGA",files) %>% 
  gsub("[.][0-9].+[_gdc_realn.txt]","",.) %>% 
  substr(.,1,16)
TCGA_need <- TCGA_read_lengths_re %>% 
  filter(sample %in% dt$samples)
TCGA_need <- left_join(
  TCGA_need,dt %>% rename(sample=samples)
)
TCGA_need <- TCGA_need %>% 
  group_by(file) %>% 
  summarise(max_len=max(read_length))

TCGA_need$depth <- NA
TCGA_need$samples <- gsub(".+[.]TCGA","TCGA",TCGA_need$file) %>% 
  gsub("[.][0-9].+[_gdc_realn.txt]","",.) %>% 
  substr(.,1,16)
TCGA_need <- TCGA_need %>% 
  mutate(is_wgs=ifelse(substr(samples,1,15) %in% pcawg_tcga_type$tcga_sample,"yes","no"))
for (i in 1:nrow(TCGA_need)){
  a <- data.table::fread(paste0("~/data/TCGA_bai/stat_files/",TCGA_need$file[i]),data.table = F)
  gl <- ifelse(TCGA_need$is_wgs[i]=="yes",3000000000,38000000)
  depth <- ((sum(a$V3))/(gl)) * TCGA_need$max_len[i]
  TCGA_need$depth[i] <- depth
}

TCGA_need$cancer <- EasyBioinfo::get_cancer_type(TCGA_need$samples)
saveRDS(TCGA_need,file = "data/TCGA_sample_depth.rds")
####




