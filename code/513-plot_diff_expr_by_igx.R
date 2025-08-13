# Scatter plot of genes from the differential expression analyses of the IGX in the two datasets.

# point colors: merged P-value of the interaction term in differential expression analyses
# x/y axes: IGX interaction coefficient (log2FC)

# Apply directional P-value merging by ActivePathways (DMP)

message('loading libraries')
library(data.table)
library(ggplot2)
library(ggrepel)
library(ActivePathways)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

# load the list of Cancer Gene Census
CGC_genes <- fread(paste0(projDir, '/input/Cosmic_CancerGeneCensus_v102_GRCh37.tsv'))[['GENE_SYMBOL']]

# load the combined table of differential expression analyses
DEA <- readRDS(paste0(projDir, '/intermediate/combined_tcga_and_metabric_diff_expr.rds'))
DEA <- DEA[analysis == 'IGX']

# genes with P = 1 and log2FC == 0 are the ones with P = NA in the original diff expr results.
# remove these genes from the analysis here.
r <- DEA[P == 1 & log2FC == 0]$gene
DEA <- DEA[! gene %in% r]

# the table to hold the information for plot
dt <- dcast(DEA, 'gene ~ data', value.var = 'log2FC', fill = NA)

# matrix of log2FC directions
dir_mtx <- as.matrix(dt[,-'gene'], rownames = dt$gene)
dir_mtx <- sign(dir_mtx)

# matrix of P-values
p_mtx <- dcast(DEA, 'gene ~ data', value.var = 'P')
p_mtx <- as.matrix(p_mtx[,-'gene'], rownames = p_mtx$gene)

dt$P <- merge_p_values(scores = p_mtx, 
                       method = 'DPM', 
                       scores_direction = dir_mtx, 
                       constraints_vector = c(1, 1))
dt$FDR <- p.adjust(dt$P, method = 'BH')

message(paste0('Found ', nrow(dt[FDR <= 0.1]), ' Genes with merged FDR <= 0.1'))

scatter_plot <- function(dt, 
                         
                         x = 'TCGA',
                         y = 'METABRIC',
                         
                         xlab = 'Interaction Coefficient, TCGA',
                         ylab = 'Interaction Coefficient, METABRIC',
                       
                          xmin=NA, 
                          ymin=NA,
                          xmax=NA, 
                          ymax=NA,
                       
                          P_cap=NA, 
                       
                          P_thrsh_for_color, 
                          FDR_thrsh_for_shape,
                          FDR_thrsh_for_label_CGC_gene,
                          P_thrsh_for_label_gene, 
                         
                          random.seed=123){
    
    set.seed(random.seed)
    
    dt$x <- dt[[x]]
    dt$y <- dt[[y]]
    
    if(is.na(xmin)){ xmin <- min(dt$x) }; dt[x < xmin]$x <- xmin
    if(is.na(ymin)){ ymin <- min(dt$y) }; dt[y < ymin]$y <- ymin
    if(is.na(xmax)){ xmax <- max(dt$x) }; dt[x > xmax]$x <- xmax
    if(is.na(ymax)){ ymax <- max(dt$y) }; dt[y > ymax]$y <- ymax
    xmin <- xmin - (xmax-xmin)/50
    xmax <- xmax + (xmax-xmin)/50
    ymin <- ymin - (ymax-ymin)/50
    ymax <- ymax + (ymax-ymin)/50
    
    dt$P_for_color <- dt$P
    if(!is.na(P_cap)) dt[P < P_cap]$P_for_color <- P_cap
    dt[P > P_thrsh_for_color]$P_for_color <- NA
    
    dt$shape <- dt$FDR < FDR_thrsh_for_shape
    
    dt$label <- paste0("italic('", dt$gene, ifelse(dt$gene %in% CGC_genes, '*', ''), "')")
    dt[!((FDR < FDR_thrsh_for_label_CGC_gene & gene %in% CGC_genes) | P < P_thrsh_for_label_gene)]$label <- NA
    
    # order points so the stronger P's are printen later 
    dt <- dt[order(P, decreasing = T)]
    
    col_brks <- seq(2, -log10(P_cap), length.out = 4)
    col_labs <- col_brks
    col_labs[length(col_labs)] <- paste0(col_labs[length(col_labs)], '+')
    
    pnsz <- 2
    lbsz <- 3
    lnth <- 0.25
    ggplot()+
        geom_point(data=dt[is.na(label)], aes(x, y, fill=-log10(P_for_color), shape=shape), size=pnsz, color='lightgrey', stroke=0.2)+
        geom_point(data=dt[!is.na(label)], aes(x, y, fill=-log10(P_for_color), shape=shape), size=pnsz, stroke=lnth, color='black')+    
        scale_fill_distiller(palette='Spectral', 
                             na.value='lightgrey', 
                             breaks = col_brks, 
                             labels = col_labs, 
                             name = 'Merged P,\n-log10')+ 
        scale_shape_manual(values = c('TRUE'=24, 'FALSE'= 21), labels = c('TRUE'='TRUE', 'FALSE'='FALSE'), name='Significant,\nFDR < 0.05')+
        geom_text_repel(data=dt[!is.na(label)], aes(x, y, label=label), parse = TRUE, color = 'black', 
                        size = lbsz, max.overlaps = nrow(dt), segment.size=lnth, fontface=3)+ 
        guides(shape=guide_legend(reverse = T))+
        scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax))+
        scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax))+
        coord_fixed(ratio = (xmax-xmin)/(ymax-ymin))+
        xlab(xlab)+
        ylab(ylab)+
        theme_bw()+
        theme(panel.grid.minor = element_blank())
}

pdf(paste0(projDir, '/output/Fig-5E.pdf'), 7/1.1, 5/1.1)
scatter_plot(dt, 
             xmin=-5, 
             xmax=5,
             ymin=-2, 
             ymax=2,
             P_cap=1e-8, 
             P_thrsh_for_color=0.1, 
             FDR_thrsh_for_shape=0.05,
             FDR_thrsh_for_label_CGC_gene=0.1,
             P_thrsh_for_label_gene=1e-6)
invisible(dev.off())
