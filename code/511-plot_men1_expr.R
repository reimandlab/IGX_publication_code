# Plot MEN1 expression with respect to "11q13.1 loss" in the two datasets

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

message('loading libraries')
library(data.table)
library(ggplot2)

message('load expression data for TCGA')
tpm <- readRDS(paste0(projDir, '/input/expr_tpm_20_genes_in_tcga.rds'))
tpm <- as.data.table(tpm)
tpm <- t(tpm[gene_name == 'MEN1', -'gene_name'])
dt_tcga <- readRDS(paste0(projDir, '/intermediate/tcga_all_datasets.rds'))[['BRCA__LumA']][['cna']]$ft
dt_tcga$G <- as.numeric(dt_tcga$Del_11q13.1)-1
stopifnot(all(dt_tcga$patient_id %in% rownames(tpm)))
dt_tcga$MEN1_expr <- tpm[dt_tcga$patient_id, ]
dt_tcga$MEN1_log2expr <- log2(dt_tcga$MEN1_expr)
dt_tcga$data <- 'TCGA'

message('load expression data for METABRIC')
expr <- fread(paste0(projDir, '/input/brca_metabric/data_mrna_illumina_microarray.txt'))
expr <- t(expr[Hugo_Symbol == 'MEN1', -c('Hugo_Symbol', 'Entrez_Gene_Id')])
dt_metabric <- readRDS(paste0(projDir, '/intermediate/metabric_two_datasets.rds'))[['LumA']]
stopifnot(all(dt_metabric$PATIENT_ID %in% rownames(expr)))
dt_metabric$MEN1_log2expr <- expr[dt_metabric$PATIENT_ID, ]
dt_metabric$data <- 'METABRIC'

# combine data
dt <- rbind(dt_tcga[, c('data', 'G', 'MEN1_log2expr')], dt_metabric[, c('data', 'G', 'MEN1_log2expr')])
dt$data <- relevel(factor(dt$data), 'TCGA')
levs <- c('11q13.1 no loss', '11q13.1 loss')
dt$G <- levs[1+dt$G]
dt$G <- factor(dt$G, levels = levs)

ast <- Vectorize(function(x){
    if(x <= 0.001) return('***')
    if(x <= 0.01) return('**')
    if(x <= 0.05) return('*')
    return('ns')
})
DEA <- readRDS(paste0(projDir, '/intermediate/combined_tcga_and_metabric_diff_expr.rds'))
DEA <- DEA[analysis == 'G' & gene == 'MEN1']
DEA$significance <- ast(DEA$P)
sig.dt <- merge(dt[, list(x=1.5, y = max(MEN1_log2expr)), by='data'], DEA, by='data')
sig.dt$data <- factor(sig.dt$data, levels = levels(dt$data))

plt <- ggplot()+theme_bw()+
    geom_boxplot(data=dt, aes(x=G, y=MEN1_log2expr, fill=G))+
    scale_fill_manual(values = c('gray', '#f25550'), name=NULL)+
    geom_text(data=sig.dt, aes(x=x, y=1.03*y, label=significance), color='red', size=4)+
    facet_wrap('.~data', scales = 'free')+
    ylab('Gene expression, log2')+
    theme(strip.text.x = element_text(size = 11),
          panel.grid = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line())

pdf(paste0(projDir, '/output/Fig-5C.pdf'), 4, 3.5)
plot(plt)
invisible(dev.off())
