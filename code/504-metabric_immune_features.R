# After running CIBERSORTx for "metabric_lm22expr.tsv", save the results in "intermediate/metabric_cibersort_results.txt".

# This script makes the table of immune features for METABRIC samples based on the CIBERSORTx results.

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

message('loading libraries')
library(data.table)

imm <- fread(paste0(projDir, '/intermediate/metabric_cibersort_results.txt'))

imm <- imm[,-c('P-value','Correlation','RMSE')]
colnames(imm) <- gsub(' ', '_', gsub('\\(', '', gsub('\\)', '', colnames(imm))))
imm$Mixture <- gsub('\\.','-',imm$Mixture)

setnames(imm, 'Mixture', 'PATIENT_ID')

saveRDS(imm, paste0(projDir, '/intermediate/metabric_immune_features.rds'))
