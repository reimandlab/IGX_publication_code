# Plot the P-value and logHR of IGXs vs the best of G or I for each IGX

message('loading libraries')
library(data.table)
library(ggplot2)
library(patchwork)

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

# ***** use simple control model (clinical only) for this figure
dt <- do.call(rbind, lapply(1:length(catalog), function(i){
    ct <- catalog[[i]]
    hr.dt <- ct$hr.dt
    data.table(indx = i, 
               cohort = gsub('__', '.', ct$cohort), 
               nlog10.FDR_IGX_simple = -log10(ct$FDR_IGX_other_controls[3]), # *****
               nlog10.FDR_G = -log10(ct$FDR_G),
               nlog10.FDR_I = -log10(ct$FDR_I),
               logHR_IGX_simple = ct$logHR_IGX_in_simple_model_controlling_only_for_C,   # *****
               logHR_G = hr.dt[facet == 'G'][v == 'G']$logHR,
               logHR_I = hr.dt[facet == 'I'][v == 'I']$logHR)
}))

dt_nlog10.FDR <- dt[, list(y = nlog10.FDR_IGX_simple, x = max(nlog10.FDR_G, nlog10.FDR_I)), by=c('indx', 'cohort')]
dt_logHR <- dt[, list(y = logHR_IGX_simple, x = max(logHR_G, logHR_I)), by=c('indx', 'cohort')]

lim <- ceiling(max(unlist(dt_nlog10.FDR[, c('x', 'y')])))
plt1 <- ggplot()+theme_bw()+
    geom_abline(slope = 1, linetype = 'longdash', color='grey', linewidth=0.3)+
    geom_point(data=dt_nlog10.FDR, aes(x=x, y=y, fill=cohort), shape=21, stroke=0.2, size=1.75)+
    scale_fill_manual(values = cohort.colors, name='Cancer type')+
    guides(fill = guide_legend(override.aes = list(shape = 22, size=5, stroke=0.3)))+
    scale_x_continuous(limits = c(0, lim))+
    scale_y_continuous(limits = c(0, lim))+
    xlab('FDR, strongest of G or I (-log10)')+
    ylab('FDR of IGX (-log10)')+
    coord_fixed(ratio = 1)+
    theme(panel.grid = element_blank())

lim <- ceiling(20*max(abs(range(unlist(dt_logHR[, c('x', 'y')])))))/20
plt2 <- ggplot()+theme_bw()+
    geom_vline(xintercept = 0, linetype = 'longdash', color='grey', linewidth=0.3)+
    geom_hline(yintercept = 0, linetype = 'longdash', color='grey', linewidth=0.3)+
    geom_abline(slope = 1, linetype = 'longdash', color='grey', linewidth=0.3)+
    geom_point(data=dt_logHR, aes(x=x, y=y, fill=cohort), shape=21, stroke=0.2, size=1.75)+
    scale_fill_manual(values = cohort.colors, name='Cancer type')+
    guides(fill = guide_legend(override.aes = list(shape = 22, size=5, stroke=0.3)))+
    scale_x_continuous(limits = c(-lim, lim))+
    scale_y_continuous(limits = c(-lim, lim))+
    xlab('logHR, strongest of G or I')+
    ylab('logHR of IGX')+
    coord_fixed(ratio = 1)+
    theme(panel.grid = element_blank())

pdf(paste0(projDir, '/output/Fig-2BC.pdf'), 7.5/1.2, 4)
(plt1 | plt2) + plot_layout(guides = 'collect')
invisible(dev.off())
