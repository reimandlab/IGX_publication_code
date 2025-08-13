# Plot fraction of cancer samples in each IGX

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

dt <- do.call(rbind, lapply(1:length(catalog), function(indx){
    ct <- catalog[[indx]]
    tb <- table(ct$patient_level_dt[, c('I', 'G')])
    tb <- 100 * tb / sum(tb)
    G_lb <- gsub('\\(', '', gsub(')', '', ct$G_label))
    cohort <- gsub('__', '.', ct$cohort)
    cbind(cohort = cohort,
          IGX = paste0(cohort, ': ', G_lb, ' & ', ct$I_label), 
          prog = ifelse(ct$hr.dt[facet == 'IGX'][v == 'IGX']$logHR > 0, '-', '+'),
          data.table(tb))
}))
dt$G <- as.numeric(dt$G)
dt$I <- as.numeric(dt$I)
lvs <- c('IGX', 'G only', 'I only', 'neither')
dt$group <- factor(lvs[4-(dt$I+2*dt$G)], levels = rev(lvs))
igx_annot <- dt[, list(N=N[group == 'IGX']), by='IGX']
setorder(igx_annot, -N, IGX)
dt$IGX <- factor(dt$IGX, levels = igx_annot$IGX)

plt3 <- ggplot()+theme_bw()+
    geom_bar(stat='identity', data=dt, aes(x=IGX, y=N, fill=group), color='black', linewidth=0.3, width=1)+
    scale_fill_manual(values = c('#cbccce', '#86c6ee', '#b1d393', '#e17038'), name='Cancer\nsamples')+
    scale_x_discrete(expand = c(0, 0))+
    scale_y_continuous(expand = c(0, 0))+
    guides(fill = guide_legend(reverse = TRUE))+
    ylab('% cancer samples')+
    theme(plot.margin=margin(b=20),
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

plt1 <- ggplot()+theme_void()+
    geom_segment(data=data.table(x=as.numeric(dt$IGX)-1), aes(x=x-0.2, xend=x+0.2, y=0, yend=0), linewidth=0.3)+
    geom_segment(data=data.table(x=as.numeric(dt$IGX)-1)[x %in% (as.numeric(dt[prog == '+']$IGX)-1)], 
                 aes(x=x, xend=x, y=-1, yend=1), linewidth=0.3)+
    scale_x_continuous(expand = c(0, 0), limits = c(-0.5, length(catalog)-0.5))+
    scale_y_continuous(expand = c(0, 0))+
    theme(plot.margin = margin())

plt2 <- ggplot()+theme_void()+
    geom_tile(data=dt, aes(x=IGX, y=1, fill=cohort), color='black', linewidth=0.3)+
    scale_fill_manual(values=cohort.colors)+
    scale_x_discrete(expand = c(0, 0))+
    scale_y_continuous(expand = c(0, 0))+
    theme(plot.margin=margin(t=3, b=5), legend.position = 'n')

pdf(paste0(projDir, '/output/Fig-2D.pdf'), 7, 5.1)
(plt1 / plt2 / plt3)+plot_layout(heights = c(0.4, 0.95, 10))
invisible(dev.off())
