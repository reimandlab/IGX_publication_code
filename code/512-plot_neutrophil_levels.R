# Plot the Neutrophils levels in TCGA and METABRIC with respect to the IGX "11q13.1 loss & Neutrophils low".

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

message('loading libraries')
library(data.table)
library(ggplot2)

# remake the G and I features in tcga_brca_lumA for the IGX of interest
tcga_dt <- readRDS(paste0(projDir, '/intermediate/tcga_all_datasets.rds'))[['BRCA__LumA']][['cna']]$ft
tcga_dt$G <- as.numeric(tcga_dt$Del_11q13.1)-1
tcga_dt$I <- as.numeric(tcga_dt$Neutrophils <= median(tcga_dt$Neutrophils))

# the metabric dataset already has G and I features for this IGX
metabric_dt <- readRDS(paste0(projDir, '/intermediate/metabric_two_datasets.rds'))[['LumA']]

plt <- function(dt, imm, ttl){
    dt$imm <- 100*dt[[imm]]

    dt <- dt[,c('imm', 'I', 'G')]
    dt <- dt[order(imm, decreasing = T)]
    dt1 <- dt[I == 0 & G == 0]; n1 <- nrow(dt1)
    dt2 <- dt[I == 0 & G == 1]; n2 <- nrow(dt2)
    dt3 <- dt[I == 1 & G == 0]; n3 <- nrow(dt3)
    dt4 <- dt[I == 1 & G == 1]; n4 <- nrow(dt4)
    dt1 <- rbind(dt1, data.table(imm=0, I=NA, G=NA))
    dt <- rbind(dt1, dt2, dt3, dt4)
    dt$x <- factor(1:nrow(dt))

    ggplot()+theme_bw()+
        geom_bar(stat='identity', data=dt, aes(x=x, y=imm), fill='red', color=NA, width=0.7)+
        geom_vline(xintercept = nrow(dt1), linewidth=0.3)+
        geom_vline(xintercept = nrow(dt1)+nrow(dt2), linewidth=0.3)+
        geom_vline(xintercept = nrow(dt1)+nrow(dt2)+nrow(dt3), linewidth=0.3)+
        scale_y_continuous(expand = c(0,0))+
        scale_x_discrete(expand = c(0,1.2))+
        ylab(paste0(imm, ' ICF (%)'))+
        xlab('Tumor Samples')+
        annotate("text", x = nrow(dt1)/2, y = 0.9*max(dt$imm), label = paste0('G-\n(', n1, ')'))+
        annotate("text", x = nrow(dt1)+nrow(dt2)/2, y = 0.9*max(dt$imm), label = paste0('G+\n(', n2, ')'))+
        annotate("text", x = nrow(dt1)+nrow(dt2)+nrow(dt3)/2, y = 0.9*max(dt$imm), label = paste0('G-\n(', n3, ')'))+
        annotate("text", x = nrow(dt1)+nrow(dt2)+nrow(dt3)+nrow(dt4)/2, y = 0.9*max(dt$imm), label = paste0('G+\n(', n4, ')'))+
        annotate("text", x = (nrow(dt1)+nrow(dt2))/2, y = 1.07*max(dt$imm), label = paste0(imm, ', high'))+
        annotate("text", x = nrow(dt1)+nrow(dt2)+(nrow(dt3)+nrow(dt4))/2, y = 1.07*max(dt$imm), label = paste0(imm, ', low'))+
        coord_cartesian(ylim=c(0, 1.01*max(dt$imm)), clip = "off")+
        ggtitle(paste0(ttl, '\n'))+
        theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

pdf(paste0(projDir, '/output/Fig-5D.pdf'), 10, 3)
plt(tcga_dt, 'Neutrophils', 'TCGA')
plt(metabric_dt, 'Neutrophils', 'METABRIC')
invisible(dev.off())
