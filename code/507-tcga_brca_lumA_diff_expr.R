# Differential gene expression analysis in TCGA-BRCA-LumA with respect to the IGX "11q13.1 loss & Neutrophils low" and the genomic feature "11q13.1 loss" alone.

message('loading libraries')
suppressMessages({
    library(data.table)
    library(DESeq2)
})

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

message('loading expression data')

dt <- readRDS(paste0(projDir, '/intermediate/tcga_all_datasets.rds'))[['BRCA__LumA']][['cna']]$ft
dt$G <- as.factor(as.numeric(dt$Del_11q13.1)-1)
dt$I <- as.factor(as.numeric(dt$Neutrophils <= median(dt$Neutrophils)))
sample_info <- data.frame(dt[, c('G', 'I')], row.names = dt$patient_id)

counts <- readRDS(paste0(projDir, '/input/expr_counts_tcga_brca_lumA.rds'))
counts <- as.data.table(counts)
gene_info <- counts[, c('gene_id', 'gene_name')]
gene_info$uid <- paste0('g', 1:nrow(gene_info))
counts_mtx <- as.matrix(counts[, -c('gene_id', 'gene_name')])
rownames(counts_mtx) <- gene_info$uid

stopifnot(all(rownames(sample_info) %in% colnames(counts_mtx)))
counts_mtx <- counts_mtx[, rownames(sample_info)]

stopifnot(0 == sum(is.na(counts_mtx)))

DEA <- function(design.formula, coefficient){
    dds <- DESeqDataSetFromMatrix(countData = counts_mtx, colData = sample_info, design = as.formula(design.formula))
    dds <- DESeq(dds, quiet = T)
    res <- data.table(data.frame(results(dds, name = coefficient)), keep.rownames = 'uid')
    stopifnot(setequal(res$uid, gene_info$uid))
    res <- merge(gene_info, res, by = 'uid', all = T)[, -'uid'][order(pvalue)]
    return(res)
}

message('differential expression analysis')
res.list <- list('IGX' = DEA('~G*I', 'G1.I1'), 'G' = DEA('~G', 'G_1_vs_0'))
saveRDS(res.list, paste0(projDir, '/intermediate/tcga_brca_lumA_diff_expr.rds'))
