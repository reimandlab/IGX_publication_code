message('loading libraries')
library(data.table)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

pths <- list.files(paste0(projDir, '/input/tcga_cna_gistic'), pattern = 'all_lesions.conf_99.txt', recursive = T, full.names = T)

dt <- do.call(rbind, lapply(pths, function(pth){
    
    cna <- fread(pth)
    cna <- cna[!grepl(' - CN values', `Unique Name`)]
    cna$feature <- paste0(ifelse(grepl('Deletion', cna$`Unique Name`), 'Del_', 'Amp_'), cna$Descriptor)
    
    # in case of duplicate events, keep the one with smaller q value
    cna$id <- 1:nrow(cna)
    select_ids <- cna[, list(id = id[which.min(`q values`)]), by='feature']$id
    cna <- cna[id %in% select_ids, -'id']    
    cna <- t(data.frame(cna[, grep('TCGA-', colnames(cna), value=T), with=F], row.names = cna$feature, check.names = F))
    cna[cna > 1] <- 1
    cna <- data.table(cna, keep.rownames = 'sample_id')
    
    # retain only samples of primary solid tumors
    cna <- cna['01' == substr(cna$sample_id, 14, 15)]

    # sample_id --> patient_id
    setnames(cna, old = 'sample_id', new = 'patient_id')
    cna$patient_id <- substr(cna$patient_id, 1, 12)

    # make sure each patient has only one sample
    stopifnot(all(!duplicated(cna$patient_id)))

    cna <- melt(cna, id.vars = 'patient_id', variable.name = 'feature', variable.factor = F, value.name = 'event')
    
    return(cna)
}))

dt <- dcast(unique(dt), patient_id ~ feature, value.var = 'event', fill = NA)

# this is a simple check ...
tss2study <- merge(fread(paste0(projDir, '/input/tcga_code_tables/tissueSourceSite.tsv')),
                   fread(paste0(projDir, '/input/tcga_code_tables/diseaseStudy.tsv')), 
                   by = 'Study Name', all=TRUE)
tss2study <- tss2study[0 == rowSums(is.na(tss2study))]
tss2study <- setNames(tss2study[['Study Abbreviation']], tss2study[['TSS Code']])
n <- table(tss2study[substr(dt$patient_id, 6, 7)])
stopifnot(all(sapply(names(n), function(s){
    identical(names(table(colSums(is.na(dt[tss2study[substr(dt$patient_id, 6, 7)] == s])))), as.character(c(0, n[s])))
})))

saveRDS(dt, paste0(projDir, '/intermediate/cna_features.rds'))
