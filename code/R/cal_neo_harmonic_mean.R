cal_harmonic <- function(file, sep, header = T, 
                        stringsAsFactors = F, 
                        out = '~/harmonic.rds'){
    df <- data.table::fread(file = file, sep = sep, 
                        stringsAsFactors = stringsAsFactors, 
                        header = header)
    df <- df %>% 
        filter(`Gene Name` != "") %>%
    mutate(
        index = paste(Tumor_Barcode, Chromosome, Stop, sep = ':')
    ) %>% 
    group_by(index, `HLA Allele`) %>% 
    summarise(
        MT_aff = min(`Best MT Score`),
        WT_aff = `Corresponding WT Score`[which(`Best MT Score` == min(`Best MT Score`))][1]
    ) %>% 
    group_by(index) %>% 
    summarise(
        MT_har_mean = harmonic_mean(MT_aff),
        WT_har_mean = harmonic_mean(WT_aff)
    ) %>% 
    ungroup()

    saveRDS(df, file = out)
}