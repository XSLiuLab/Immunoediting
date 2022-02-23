library(DoAbsolute )

meta <- read.table("~/neo_dep/immunethearpy_dataset/Nadeem/On_sample_absolute/on_norm_pat",sep = ",")
seg <- data.table::fread("~/neo_dep/immunethearpy_dataset/Nadeem/On_sample_absolute/seg/SRR5134744.cr.igv.seg")

files <- list.files(path = "~/neo_dep/immunethearpy_dataset/Nadeem/On_sample_absolute/maf/")
maf <- do.call(rbind,
               lapply(paste0("~/neo_dep/immunethearpy_dataset/Nadeem/On_sample_absolute/maf/",files), data.table::fread))
colnames(maf)[6] <- "Start_position"

files <- list.files(path = "~/neo_dep/immunethearpy_dataset/Nadeem/On_sample_absolute/seg/")
seg <- do.call(rbind,
               lapply(paste0("~/neo_dep/immunethearpy_dataset/Nadeem/On_sample_absolute/seg/",files), data.table::fread))
DoAbsolute(Seg = seg, Maf = maf, platform = "Illumina_WES", copy.num.type = "total",
           results.dir = "~/neo_dep/immunethearpy_dataset/Nadeem/On_sample_absolute/out/", nThread = 10, keepAllResult = TRUE, verbose = TRUE)
###process
files <- list.files("~/neo_dep/immunethearpy_dataset/Nadeem/On_sample_absolute/out/seg/")
ccf <-  do.call(rbind,
                lapply(paste0("~/neo_dep/immunethearpy_dataset/Nadeem/On_sample_absolute/out/seg/",files), data.table::fread))