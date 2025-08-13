# Pathway enrichment analyses of the IGX in TCGA-BRCA-LumA and METABRIC-LumA.

# Use the combined table of the differential expression analyses.

# Exclude LM22 genes upstream.

# activepathways background is restricted to the genes present in this table. 


message('loading libraries')
library(data.table)
library(ActivePathways)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

message('Loading the combined table of differential expression analyses')

DEA <- readRDS(paste0(projDir, '/intermediate/combined_tcga_and_metabric_diff_expr.rds'))

# focus on diff expr analysis with respect to "IGX" effect
DEA <- DEA[analysis == 'IGX']

# exclude LM22 genes
DEA <- DEA[LM22 == FALSE]

# matrix of P-values
p_mtx <- dcast(DEA, 'gene ~ data', value.var = 'P')
p_mtx <- as.matrix(p_mtx[,-'gene'], rownames = p_mtx$gene)

# matrix of log2FC directions
dir_mtx <- dcast(DEA, 'gene ~ data', value.var = 'log2FC', fill = NA)
dir_mtx <- as.matrix(dir_mtx[,-'gene'], rownames = dir_mtx$gene)
dir_mtx <- sign(dir_mtx)

# Set the background for pathway enrichment analysis.
# Background genes are initially defined as all genes from the GMT, and then 
# restricted to the intersection with genes present in the diff expression analyses.
gmt <- read.GMT(paste0(projDir, '/input/hsapiens.REAC_GOBP.name.gmt'))
background <- makeBackground(gmt)
this_background <- intersect(background, rownames(p_mtx))

colors <- c('METABRIC' = rgb(255, 214, 54,  maxColorValue = 255), 
            'TCGA'     = rgb(0,   151, 71,  maxColorValue = 255),
            'combined' = rgb(202, 204, 205, maxColorValue = 255))

cytoscape_file_tag <- paste0(paste0(projDir, '/output/enrichment_map/'))
dir.create(cytoscape_file_tag)

message('Running ActivePathways')

res <- ActivePathways(scores = p_mtx, 
                      gmt = gmt, 
                      background = this_background, 
                      cutoff = 0.1, 
                      merge_method = 'DPM',
                      scores_direction = dir_mtx,
                      constraints_vector = c(1,1),
                      correction_method = 'fdr',
                      significant = 0.01,
                      geneset_filter = c(20, 1000),
                      cytoscape_file_tag = cytoscape_file_tag,
                      custom_colors = colors[c('METABRIC', 'TCGA')],
                      color_integrated_only = colors['combined'])

saveRDS(res, paste0(projDir, '/intermediate/pathway_enrichment.rds'))

