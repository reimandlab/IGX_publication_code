message('loading libraries')
library(data.table)
library(readxl)
library(maftools)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

# 299 driver genes
pth <- paste0(projDir, '/input/NIHMS948705-supplement-8.xlsx')
driver_genes <- read_xlsx(pth, sheet = 'Table S1', .name_repair = 'minimal')
driver_genes <- setdiff(driver_genes[[1]], c(NA, 'Gene'))

# TCGA SNVs
message('loading mc3.v0.2.8.PUBLIC.maf')
snvs <- read.maf(maf = paste0(projDir, '/input/mc3.v0.2.8.PUBLIC.maf'))
snvs <- snvs@data
setnames(snvs, old = c('Tumor_Sample_Barcode', 'Hugo_Symbol'), new = c('sample_id', 'gene'))

# check that all driver_genes exist in the maf
stopifnot(all(driver_genes %in% snvs$gene))

select_variants <- c('Translation_Start_Site',
                     'Missense_Mutation',
                     'Nonsense_Mutation',
                     'Nonstop_Mutation',
                     'Frame_Shift_Ins',
                     'Frame_Shift_Del',
                     'In_Frame_Ins',
                     'In_Frame_Del')

# make the table of binary mutation features, where events are based on selected variants driver genes
indx <- (snvs$Variant_Classification %in% select_variants) & (snvs$gene %in% driver_genes)
dt1 <- cbind(unique(snvs[indx, c('sample_id', 'gene')]), mut = 1)
dt0 <- expand.grid(sample_id = unique(snvs$sample_id), gene = driver_genes, mut = 0, stringsAsFactors = F)
dt <- rbind(dt1, dt0)
dt <- dt[, list(mut = sum(mut)), by=c('sample_id', 'gene')]
dt <- dcast(dt, 'sample_id ~ gene', value.var = 'mut')

# retain only samples of primary solid tumors
dt <- dt['01' == substr(dt$sample_id, 14, 15)]

# sample_id --> patient_id
setnames(dt, old = 'sample_id', new = 'patient_id')
dt$patient_id <- substr(dt$patient_id, 1, 12)

# make sure each patient has only one sample
stopifnot(all(!duplicated(dt$patient_id)))

saveRDS(dt, paste0(projDir, '/intermediate/tcga_mutation_features.rds'))
