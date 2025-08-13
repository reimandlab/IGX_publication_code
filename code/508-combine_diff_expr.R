# Make a unified table for differential expression analyses of the IGX in TCGA-BRCA-LumA and METABRIC-LumA.
# Retain only genes that are shared between the two datasets
# Also exclude the small number of genes that have duplicated names in either of the datasets.

# Tag LM22 genes, as they should be removed in the pathway enrichment analysis.

message('loading libraries')
library(data.table)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

# load LM22 genes and add as many aliases as you could find.
LM22_genes <- fread(paste0(projDir, '/input/LM22.txt'))[['Gene symbol']]
LM22_genes <- union(LM22_genes, 
                    c('PGGHG', # ATHL1
                      'RADX', # CXorf57
                      'ADGRE1', # EMR1
                      'ADGRE2', # EMR2
                      'ADGRE3', # EMR3
                      'FCMR', # FAIM3
                      'GASK1B', # FAM198B
                      'INKA2', # FAM212B
                      'RIPOR2', # FAM65B
                      'ADGRG3', # GPR97
                      'GUSB', # GUSBP11
                      'H2AC8', # HIST1H2AE
                      'H2BC8', # HIST1H2BG
                      'RUBCNL', # KIAA0226L
                      'MACF1', # KIAA0754
                      'IRAG2', # LRMP
                      'MARCHF3', # MARCH3
                      'SEPTIN5', # SEPT5
                      'SEPTIN8' # SEPT8
                     ))

combine_dea <- function(tcga_dea, metabric_dea){
    
    tcga_dea <- cbind(data='TCGA', tcga_dea[, c('gene_name', 'log2FoldChange', 'pvalue')])
    metabric_dea <- cbind(data='METABRIC', metabric_dea[, c('Hugo_Symbol', 'logFC', 'P.Value')])
    
    colnames(tcga_dea) <- c('data', 'gene', 'log2FC', 'P')
    colnames(metabric_dea) <- c('data', 'gene', 'log2FC', 'P')

    # For any row with P = NA, set P to 1 and log2FC to 0.
    i <- is.na(tcga_dea$P); tcga_dea[i]$P <- 1; tcga_dea[i]$log2FC <- 0
    i <- is.na(metabric_dea$P); metabric_dea[i]$P <- 1; metabric_dea[i]$log2FC <- 0
    
    # retain only rows with unique gene name in each dataset
    tcga_dea <- tcga_dea[gene %in% names(which(1 == table(tcga_dea$gene)))]
    metabric_dea <- metabric_dea[gene %in% names(which(1 == table(metabric_dea$gene)))]

    # concatenate the tables
    dea <- rbind(tcga_dea, metabric_dea)

    # retain only genes that are shared between the two datasets
    dea <- dea[gene %in% names(which(2 == table(dea$gene)))]
    
    # tag genes that are in LM22
    dea$LM22 <- dea$gene %in% LM22_genes

    # No NA exists
    stopifnot(0 == sum(is.na(dea)))
    
    return(dea)
}


tcga_dea_list <- readRDS(paste0(projDir, '/intermediate/tcga_brca_lumA_diff_expr.rds'))
metabric_dea_list <- readRDS(paste0(projDir, '/intermediate/metabric_diff_expr.rds'))

res <- rbind(cbind(analysis = 'IGX', combine_dea(tcga_dea_list$IGX, metabric_dea_list$IGX)),
             cbind(analysis = 'G', combine_dea(tcga_dea_list$G, metabric_dea_list$G)))

saveRDS(res, paste0(projDir, '/intermediate/combined_tcga_and_metabric_diff_expr.rds'))
