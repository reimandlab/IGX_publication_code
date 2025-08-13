# Plot expression of genes contributing to enrichment of each pathway shown in the enrichment map

message('loading libraries')
library(data.table)
library(stringr)
library(ggplot2)
library(ggplotify)
library(egg)
library(patchwork)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------


map_terms <- list(
'Cellular transport'=c( # 9
'GO:0006835',
'GO:0014047',
'GO:0015800',
'GO:0051050',
'GO:0051938',
'REAC:R-HSA-983712',
'GO:0015711',
'REAC:R-HSA-2672351',
'GO:0098657')
,
'Neurotransmitter release'=c( # 12
'GO:0099500',
'GO:0099003',
'GO:0006836',
'GO:0099643',
'GO:0099504',
'REAC:R-HSA-112310',
'GO:0007269',
'GO:0031629',
'GO:0048791',
'GO:0016079',
'REAC:R-HSA-212676',
'GO:0048499')
,
'Synaptic signaling'=c( # 7
'GO:0050806',
'GO:0007268',
'GO:0050804',
'GO:0098916',
'GO:0099536',
'GO:0099177',
'GO:0099537')
,
'Action potential'=c( # 4
'GO:0060079',
'REAC:R-HSA-112316',
'GO:0098900',
'GO:0042391')
,
'Metabolic process'=c( # 8
'GO:0120254',
'GO:0043436',
'GO:0032787',
'GO:0006082',
'REAC:R-HSA-196807',
'GO:0019369',
'GO:0019752',
'GO:0033559')
, 
'Secretion'=c( # 17
'GO:1903530',
'GO:0023061',
'GO:0017158',
'GO:0042886',
'GO:0046903',
'GO:0006887',
'GO:0017156',
'GO:0045055',
'GO:0051046',
'GO:0030073',
'GO:0015833',
'GO:0140029',
'GO:0017157',
'GO:0010817',
'GO:0032940',
'GO:0140352',
'GO:0002790')
,
'Immune response'=c( # 14
'GO:0002544',
'GO:0019730',
'GO:0002275',
'GO:0006953',
'GO:0002526',
'GO:0002675',
'REAC:R-HSA-6783783',
'GO:0002250',
'GO:0034341',
'GO:0071346',
'GO:0002523',
'GO:0002269',
'GO:0002673',
'REAC:R-HSA-1280215')
,
'Antigen presentation'=c( # 5
'GO:0002478',
'GO:0019884',
'GO:0019882',
'GO:0002495',
'GO:0002504')
,
'Response to stimulus'=c( # 5
'GO:1901652',
'GO:0034097',
'GO:0032102',
'GO:0071345',
'GO:0032496')
,
'Detoxification'=c( # 3 
'GO:0098869',
'GO:1990748',
'GO:0097237')
)

# there are 84 nodes in the enrichment map
stopifnot(84 == length(unique(unlist(map_terms))))

# load activepathways output
AP <- readRDS(paste0(projDir, '/intermediate/pathway_enrichment.rds'))

# add functional themes to the AP table
AP$functional_theme <- ''
for(theme in names(map_terms)) AP[term_id %in% map_terms[[theme]]]$functional_theme <- theme

# sanity check
stopifnot(all(unlist(map_terms) %in% AP$term_id))

# load the list of Cancer Gene Census
CGC_genes <- fread(paste0(projDir, '/input/Cosmic_CancerGeneCensus_v102_GRCh37.tsv'))[['GENE_SYMBOL']]

# load the list of OncoKB genes
OncoKB_genes <- fread(paste0(projDir, '/input/OncoKB_cancerGeneList.tsv'))[['Hugo Symbol']]

# load the diff. expr. resutls of both datasets
DEA <- readRDS(paste0(projDir, '/intermediate/combined_tcga_and_metabric_diff_expr.rds'))
DEA <- DEA[analysis == 'IGX']

# function for bubble plot
bubble_plt <- function(node, max_sig=6, min_log2FC=-2, max_log2FC=2, hide_legend=FALSE){
    dt <- DEA[gene %in% unlist(AP[term_id == node]$overlap)]
    dt$data <- factor(dt$data, levels = c('METABRIC', 'TCGA'))
 
    dt$sig <- -log10(dt$P)
    
    # dt$gene <- factor(dt$gene, levels = str_sort(unique(dt$gene), numeric = T))
    dt$gene <- factor(dt$gene, levels = dt[, mean(sig*log2FC), by='gene'][order(V1, decreasing = T)]$gene)
    
    dt[sig > max_sig]$sig <- max_sig
    dt[log2FC < min_log2FC]$log2FC <- min_log2FC
    dt[log2FC > max_log2FC]$log2FC <- max_log2FC
    
    gene_color <- rep('black', nrow(dt))
    names(gene_color) <- dt$gene
    gene_color[names(gene_color) %in% OncoKB_genes] <- '#08a633'
    gene_color[names(gene_color) %in% CGC_genes] <- 'red'
    
    gene_face <- rep('plain', nrow(dt))
    names(gene_face) <- dt$gene
    gene_face[names(gene_face) %in% c(CGC_genes, OncoKB_genes)] <- 'bold'
    
    # arrange the gene annotations in the order of appearing on the x-axis
    gene_color <- gene_color[levels(dt$gene)]
    gene_face <- gene_face[levels(dt$gene)]
    
    plt <- ggplot()+theme_bw()+
        geom_point(data=dt, aes(x=gene, y=data, size=sig, fill=log2FC), color='black', pch=21)+
        scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', 
                             limits=c(min_log2FC, max_log2FC), 
                             name ='Effect size (log2 FC)')+
        scale_size(limits = c(0, max_sig), 
                   breaks = c(1:max_sig), 
                   name = 'P-value (-log10)')+
        guides(size = guide_legend(reverse = T, nrow = 1, label.position = 'bottom', keywidth=0))+
        ggtitle(paste0('Theme: ', AP[term_id == node]$functional_theme,
                       '\nNode: ', AP[term_id == node]$term_name, '  -  ', node, 
                       '\nFDR = ', signif(AP[term_id == node]$adjusted_p_val, 2)))+
        suppressWarnings( # this is only because of "color" and "face" in axis.text.x. Will be updated soon.
            theme(axis.text.y = element_text(color = 'black'),
                  axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, color = gene_color, face = gene_face),
                  panel.grid = element_blank(),
                  axis.title = element_blank(), 
                  legend.position = ifelse(hide_legend, 'none', 'bottom'), 
                  legend.direction = "vertical",
                  legend.key.spacing.x = unit(0, 'lines'))
        )
    nx <- length(levels(dt$gene))
    ny <- length(levels(dt$data))
    as.ggplot(set_panel_size(plt, width  = unit(nx/2.1, 'cm'), height = unit(ny/2.1, 'cm')))
}



# for(theme in names(map_terms)){
#     for(node in AP[functional_theme == theme][order(adjusted_p_val)]$term_id){
#         plot(bubble_plt(node))
#     }
# }


# plot for the selected terms
message('plotting')
pdf(paste0(projDir, '/output/Fig-5G.pdf'), 31, 4)
for(node in c('REAC:R-HSA-6783783', 'GO:0019882', 'GO:0034341')){
    plot(bubble_plt(node))
}
invisible(dev.off())
