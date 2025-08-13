# Plot expression of CNA affected genes

message('loading libraries')
library(data.table)
library(ggplot2)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

# load the catalog
all_datasets <- readRDS(paste0(projDir, '/intermediate/tcga_all_datasets.rds'))

# load expression data (TPM)
tpm <- readRDS(paste0(projDir, '/input/expr_tpm_20_genes_in_tcga.rds'))
tpm <- data.table(tpm)
tpm <- t(data.frame(tpm[,-'gene_name'], row.names = tpm[['gene_name']], check.names = F))
tpm <- data.table(tpm, keep.rownames = 'patient_id')

# hard-code - to match to Figure 2H
info <- list(list(cohort = 'SARC', G = 'Amp_1p32.1', gene = 'JUN'), 
             list(cohort = 'KIRC', G = 'Del_1p36.13', gene = 'ARID1A'),
             list(cohort = 'COAD', G = 'Del_3p26.1', gene = 'VHL'),
             list(cohort = 'KIRC', G = 'Del_4q34.3', gene = 'FBXW7'),
             list(cohort = 'BRCA__LumB', G = 'Del_4q35.2', gene = 'CASP3'),
             list(cohort = 'LUAD', G = 'Amp_5q35.1', gene = 'NPM1'),
             list(cohort = 'CESC', G = 'Amp_7p11.2', gene = 'EGFR'),
             list(cohort = 'SARC', G = 'Amp_8q24.21', gene = 'MYC'),
             list(cohort = 'KIRC', G = 'Amp_8q24.22', gene = 'COX6C'),
             list(cohort = 'BRCA__LumA', G = 'Del_11q13.1', gene = 'MEN1'),
             list(cohort = 'KIRP', G = 'Del_14q24.2', gene = 'HIF1A'),
             list(cohort = 'COAD', G = 'Del_17p12', gene = 'TP53'),
             list(cohort = 'COAD', G = 'Amp_17q12', gene = 'ERBB2'),
             list(cohort = 'COAD', G = 'Amp_17q24.1', gene = 'AXIN2'))

# extract the gene expression and patient-level annoations
dt <- do.call(rbind, lapply(info, function(f){
    dt <- all_datasets[[f$cohor]][['cna']]$ft
    dt$G <- dt[[f$G]]
    tpm$expr <- tpm[[f$gene]]
    dt <- merge(dt, tpm, by='patient_id')
    dt$cna <- factor('Balanced', levels = c('Loss', 'Gain', 'Balanced'))
    dt[G == 'yes']$cna <- ifelse(grepl('^Del_', f$G), 'Loss', 'Gain')
    dt <- dt[, c('G', 'cna', 'expr')]
    stopifnot(0 == sum(is.na(dt)))    
    dt$label <- paste0(f$gene, '\n', gsub('Amp_', '', gsub('Del_', '', f$G)), '\n', 
                       gsub('__', '.', f$cohort))
    return(dt)
}))
dt$label <- factor(dt$label, levels = unique(dt$label))

# apply U test between two groups 
ast <- Vectorize(function(x){
    if(x <= 0.001) return('***')
    if(x <= 0.01) return('**')
    if(x <= 0.05) return('*')
    return('ns')
})
p.dt <- dt[, list(P = wilcox.test(expr[G=='no'], expr[G=='yes'])$p.value), by = 'label']
p.dt$FDR <- p.adjust(p.dt$P, method = 'BH')
p.dt$sig <- ast(p.dt$FDR)

# plot
plt <- ggplot()+theme_bw()+
    geom_boxplot(data=dt, aes(x=G, y=log1p(expr), fill=cna), notch=T, linewidth=0.3)+
    geom_text(data=p.dt, aes(x=1.5, y=1.05*log1p(max(dt$expr)), label=sig))+
    facet_wrap(~label, nrow=1)+
    scale_fill_manual(values = c('#00BFC4', '#F8766D', 'grey'), name='CNA')+
    xlab('Samples grouped by CNA status of IGXs')+
    ylab('Gene expression (log1p TPM)')+
    theme(strip.background.x = element_blank(),
          panel.grid = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())


pdf(paste0(projDir, '/output/Fig-2H.pdf'), 13, 3.5)
plot(plt)
invisible(dev.off())
