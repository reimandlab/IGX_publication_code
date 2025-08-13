# Plot the expression of MEN1 vs 11q13.1 loss in the two datasets.

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

message('loading libraries')
library(data.table)
library(ggplot2)
library(limma)

message('loading data_mrna_illumina_microarray.txt')
expr_dt <- fread(paste0(projDir, '/input/brca_metabric/data_mrna_illumina_microarray.txt'))

gene_info <- expr_dt[, c('Hugo_Symbol', 'Entrez_Gene_Id')]
gene_info$uid <- paste0('g', 1:nrow(gene_info))

expr_mtx <- as.matrix(expr_dt[, -c('Hugo_Symbol', 'Entrez_Gene_Id')])
rownames(expr_mtx) <- gene_info$uid

dt <- readRDS(paste0(projDir, '/intermediate/metabric_two_datasets.rds'))[['LumA']]
sample_info <- data.frame(dt[, c('G', 'I')], row.names = dt$PATIENT_ID)
stopifnot(all(rownames(sample_info) %in% colnames(expr_mtx)))
expr_mtx <- expr_mtx[, rownames(sample_info)]

DEA <- function(design.formula, coefficient){    
    design <- model.matrix(as.formula(design.formula), data = sample_info)
    fit <- lmFit(expr_mtx, design)
    fit <- eBayes(fit)
    res <- data.table(topTable(fit, coef = coefficient, number = Inf), keep.rownames = 'uid')
    res <- merge(gene_info, res, by = 'uid', all = T)[, -'uid'][order(P.Value)]
    stopifnot(0 == sum(is.na(res)))
    return(res)
}
message('differential expression analysis')
res.list <- list('IGX' = DEA('~G*I', 'G:I'), 'G' = DEA('~G', 'G'))
saveRDS(res.list, paste0(projDir, '/intermediate/metabric_diff_expr.rds'))
