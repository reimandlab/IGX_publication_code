# Barplots logHRs per I or G featrures in each IGX.

message('loading libraries')
library(data.table)
library(ggplot2)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

cohort.colors <- 
    c('BLCA'='#E9A22E', 
      'BRCA.LumA'='#C25A84', 
      'BRCA.LumB'='#F176A7', 
      'CESC'='#73C5C5', 
      'COAD'='#002CCF', 
      'GBM'='#000000', 
      'HNSC'='#7D2322', 
      'KIRC'='#F9431D', 
      'KIRP'='#EF7727', 
      'LGG'='#6F6F6F', 
      'LIHC'='#145713', 
      'LUAD'='#92905A', 
      'LUSC'='#EBE45C', 
      'OV'='#1F9293', 
      'PAAD'='#B05920', 
      'READ'='#29B376', 
      'SARC'='#B7C98D', 
      'STAD'='#B9ECFD', 
      'UCEC'='#8527B7')

# load the catalog
catalog <- readRDS(paste0(projDir, '/intermediate/catalog.rds'))

# This matches with the info in Figure 2H.
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

# make the data.tables for barplots
dt <- data.table(sapply(c('cohort', 'G_label', 'I_label'), function(x){ sapply(catalog, '[[', x) }))
dt$logHR_IGX <- sapply(catalog, function(ct){ ct$hr.dt[facet == 'IGX'][v == 'IGX']$logHR })
dt$cohort <- gsub('__', '.', dt$cohort)
dt$I <- gsub(' low$', '', gsub(' high$', '', dt$I_label))
dt$G <- sapply(strsplit(dt$G_label, ' '), '[', 1)
dt$G.type <- sapply(strsplit(dt$G_label, ' '), '[', 2)
for(f in info){
    cna <- unlist(strsplit(f$G, '_'))[2]
    new_cna <- paste0(cna, ' (', f$gene, ')')
    this_cohort <- f$cohort
    dt[cohort == this_cohort]$G <- gsub(cna, new_cna, dt[cohort == this_cohort]$G)
}
ttl_lvls <- c('Immune cell levels (ICL)', 
              'Driver mutations (SNV, indel)',
              'Copy number losses',
              'Copy number gains')
dt$I.type <- 'I'
dt1 <- cbind(ttl = ttl_lvls[1], dt[I.type == 'I', c('cohort', 'I', 'logHR_IGX')]); colnames(dt1)[3] <- 'x'
dt2 <- cbind(ttl = ttl_lvls[2], dt[G.type == 'mut', c('cohort', 'G', 'logHR_IGX')]); colnames(dt2)[3] <- 'x'
dt3 <- cbind(ttl = ttl_lvls[3], dt[G.type == 'loss', c('cohort', 'G', 'logHR_IGX')]); colnames(dt3)[3] <- 'x'
dt4 <- cbind(ttl = ttl_lvls[4], dt[G.type == 'gain', c('cohort', 'G', 'logHR_IGX')]); colnames(dt4)[3] <- 'x'
dt <- rbind(dt1, dt2, dt3, dt4)
dt$x <- factor(dt$x, levels = dt[, list(s=sum(abs(logHR_IGX))), by=x][order(s, decreasing = T)]$x)
dt$ttl <- factor(dt$ttl, levels = ttl_lvls)

plt <- ggplot()+theme_bw()+
    geom_bar(stat='identity', data=dt, aes(x=x, y=logHR_IGX, fill=cohort), color='black', linewidth=0.3)+
    geom_hline(yintercept = 0)+
    scale_fill_manual(values = cohort.colors, name='Cancer type')+
    facet_grid(~ttl, scales = 'free_x', space='free')+
    xlab('Genomic or immune features in IGXs')+
    ylab('Hazard ratio (log)')+
    theme(panel.grid = element_blank(), strip.background.x = element_blank(), 
          axis.line = element_line(linewidth = 0.4),
          panel.background = element_blank(),
          axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
          legend.key.size = unit(0.25, "cm"),     
          legend.text = element_text(size = 6),    
          legend.title = element_text(size = 6))

pdf(paste0(projDir, '/output/Fig-3A.pdf'), 10, 3.2)
plot(plt)
invisible(dev.off())
