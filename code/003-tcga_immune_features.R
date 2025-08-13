message('loading libraries')
library(data.table)
library(readxl)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

dt <- fread(paste0(projDir, '/input/TCGA.Kallisto.fullIDs.cibersort.relative.tsv'))

dt <- dt[, -c('CancerType', 'P.value', 'Correlation', 'RMSE')]

colnames(dt) <- trimws(gsub('\\.+', '_', colnames(dt)), whitespace='_')

setnames(dt, old = 'SampleID', new = 'sample_id')
dt$sample_id <- gsub('\\.', '-', dt$sample_id)

# retain only samples of primary solid tumors
dt <- dt['01' == substr(dt$sample_id, 14, 15)]

# sample_id --> patient_id
setnames(dt, old = 'sample_id', new = 'patient_id')
dt$patient_id <- substr(dt$patient_id, 1, 12)

# make sure each patient has only one sample
dt <- dt[patient_id %in% names(which(1 == table(dt$patient_id)))]

saveRDS(dt, paste0(projDir, '/intermediate/tcga_immune_features.rds'))
