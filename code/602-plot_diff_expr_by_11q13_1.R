# Visualize the differential expression analysis results based on 11q13.1 loss in the two datasets.

message('loading libraries')
library(data.table)
library(ggplot2)
library(ggrepel)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------


# load the list of Cancer Gene Census
CGC_genes <- fread(paste0(projDir, '/input/Cosmic_CancerGeneCensus_v102_GRCh37.tsv'))[['GENE_SYMBOL']]

DEA <- readRDS(paste0(projDir, '/intermediate/combined_tcga_and_metabric_diff_expr.rds'))

# focus on diff expr analysis with respect to "G" effect
DEA <- DEA[analysis == 'G']

# list the genes overlapping with "11q13.1 loss" event as in GISTIC2
dr <- paste0(projDir, '/input/tcga_cna_gistic/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0')
f <- paste0(dr, '/del_genes.conf_99.txt')
select_genes <- fread(f)
select_genes <- select_genes[, which(select_genes[1] == '11q13.1'), with=F]
select_genes <- setdiff(unlist(select_genes[-(1:4)]), '')

# focus on these genes and then apply FDR correction on the diff expr P-values
dt <- DEA[gene %in% select_genes]
dt$FDR <- p.adjust(dt$P, method = 'BH')
dt$nlog10FDR <- -log10(dt$FDR)

# retain genes that are significant (FDR < 0.05) in both datasets
dt <- dt[gene %in% dt[, all(FDR < 0.05), by='gene'][V1==TRUE]$gene]

reshape <- function(dt, v){ dt$v <- paste0(dt$data, '_', v); dcast(dt, 'gene ~ v', value.var = v) }
dt <- merge(reshape(dt, 'nlog10FDR'), reshape(dt, 'log2FC'), by='gene')

# mean of log2FC in the two datasets is used for the size of points
dt$mean_log2FC <- apply(dt[, c('METABRIC_log2FC', 'TCGA_log2FC')], 1, mean)
sz_brks <- seq(ceiling(10*min(dt$mean_log2FC))/10, floor(10*max(dt$mean_log2FC))/10, length.out = 4)

# color points (genes) by being in CGC or not
dt$is_CGC <- relevel(factor(ifelse(dt$gene %in% CGC_genes, 'yes', 'no')), 'yes')

plt <- ggplot(data=dt, aes(x=TCGA_nlog10FDR, y=METABRIC_nlog10FDR))+theme_bw()+
    geom_point(aes(size=abs(mean_log2FC), fill=is_CGC), shape=21)+
    geom_text_repel(aes(label=gene), size=2)+
    scale_size_continuous(breaks = abs(rev(sz_brks)), labels=rev(sz_brks), name = 'Mean\nlog2FC')+
    guides(size = guide_legend(reverse = T))+
    scale_fill_manual(values=c('yes'='red', 'no'='#3fcbf2'), name='Cancer\nGene')+
    coord_fixed()+
    theme(panel.grid = element_blank())

pdf(paste0(projDir, '/output/Fig-S5.pdf'), 7, 4)
plot(plt)
invisible(dev.off())
