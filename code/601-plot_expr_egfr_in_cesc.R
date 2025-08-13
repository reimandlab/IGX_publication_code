# plot expression of EGFR vs the IGX "7p11.2 gain & CD8+ T low" in TCGA-CESC

message('loading libraries')
library(data.table)
library(ggplot2)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

# load the catalog
catalog <- readRDS(paste0(projDir, '/intermediate/catalog.rds'))

# focus on the catalog element for the IGX of inerest
ct <- catalog[[3]]
igx <- paste0('TCGA-', ct$cohort, '\n', gsub('\\(', '', gsub('\\)', '', ct$G_label)), ' & ', ct$I_label)

# double check if this is the catalog element we're interested in
stopifnot(igx == 'TCGA-CESC\n7p11.2 gain & CD8+ T low')

# load expression data (TPM)
tpm <- readRDS(paste0(projDir, '/input/expr_tpm_20_genes_in_tcga.rds'))
tpm <- data.table(tpm)
tpm <- t(data.frame(tpm[,-'gene_name'], row.names = tpm[['gene_name']], check.names = F))
tpm <- data.table(tpm, keep.rownames = 'patient_id')

# select gene
gene <- 'EGFR'

# merge igx data with expression data at patient level
tpm$expr <- tpm[[gene]]
dt <- merge(ct$patient_level_dt, tpm[, c('patient_id', 'expr')], by='patient_id')

# sanity check
stopifnot(0 == sum(is.na(dt)))

# make the 4 group levels
lvs <- c('IGX', 'G only', 'I only', 'neither')
dt$group <- lvs[4-(dt$I+2*dt$G)]#factor(, levels = rev(lvs))
t1 <- table(dt$group)
t2 <- paste0(names(t1), ' (', t1, ')')
names(t2) <- names(t1)
dt$group <- factor(t2[dt$group], levels = rev(t2[lvs]))

# two U tests to show the significance
P1 <- wilcox.test(dt[IGX==1 & G==1]$expr, dt[IGX==0 & G==1]$expr)$p.value
P2 <- wilcox.test(dt[IGX==1]$expr, dt[IGX==0]$expr)$p.value
ast <- Vectorize(function(x){
    if(x <= 0.001) return('***')
    if(x <= 0.01) return('**')
    if(x <= 0.05) return('*')
    return('ns')
})

# make the data.frames for placing the significance labels
sig.dt <- data.table(x=c(3.5, 2), 
                     y=c(1.15, 1.1)*max(log1p(dt$expr)), 
                     sig = c(ast(P1), ast(P2)))
line.dt <- rbind(data.table(x1=3, x2=4, y1=0.98*sig.dt[1]$y, y2=0.98*sig.dt[1]$y),
                 data.table(x1=2, x2=4, y1=0.98*sig.dt[2]$y, y2=0.98*sig.dt[2]$y),
                 data.table(x1=1, x2=3, y1=0.94*sig.dt[2]$y, y2=0.94*sig.dt[2]$y))
line.dt <- rbind(line.dt, 
                 data.table(x1=line.dt[1]$x1, x2=line.dt[1]$x1, y1=0.99*line.dt[1]$y1, y2=1.01*line.dt[1]$y1),
                 data.table(x1=line.dt[1]$x2, x2=line.dt[1]$x2, y1=0.99*line.dt[1]$y1, y2=1.01*line.dt[1]$y1),
                 data.table(x1=line.dt[2]$x2, x2=line.dt[2]$x2, y1=0.99*line.dt[2]$y1, y2=1.01*line.dt[2]$y1),
                 data.table(x1=line.dt[3]$x1, x2=line.dt[3]$x1, y1=0.99*line.dt[3]$y1, y2=1.01*line.dt[3]$y1),
                 data.table(x1=line.dt[3]$x2, x2=line.dt[3]$x2, y1=0.99*line.dt[3]$y1, y2=1.01*line.dt[3]$y1))
line.dt <- rbind(line.dt, data.table(x1=2, x2=2, y1=line.dt[2]$y1, y2=line.dt[3]$y1))

# plot 
plt <- ggplot()+theme_bw()+
    geom_boxplot(data=dt, aes(x=group, y=log1p(expr)), fill=NA, outliers = FALSE, linewidth=0.4)+
    geom_jitter(data=dt, aes(x=group, y=log1p(expr), fill=group), shape=21, width=0.1, stroke=0.3)+
    scale_fill_manual(values = c('#cbccce', '#86c6ee', '#b1d393', '#e17038'), name='Cancer\nsamples')+
    geom_text(data = sig.dt, aes(x=x, y=y, label=sig))+
    geom_segment(data=line.dt, aes(x=x1, xend=x2, y=y1, yend=y2), linewidth=0.3)+
    guides(fill = guide_legend(reverse = TRUE))+
    ggtitle(igx)+
    ylab(paste0(gene, ' expression (TPM, log1p)'))+
    theme(panel.grid = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())

pdf(paste0(projDir, '/output/Fig-S4.pdf'), 4, 4)
plot(plt)
invisible(dev.off())
