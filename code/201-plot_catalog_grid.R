# Plot the grid of the catalog

message('loading libraries')
library(data.table)
library(ggplot2)
library(stringr)

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

dt <- data.table(sapply(c('cohort', 'G_label', 'I_label'), function(x){ sapply(catalog, '[[', x) }))
dt$logHR_IGX <- sapply(catalog, function(ct){ ct$hr.dt[facet == 'IGX'][v == 'IGX']$logHR })
dt$cohort <- gsub('__', '.', dt$cohort)

dt$type <- as.character(NA)
tlvs <- c()
x <- 'lower-risk IGX; ICL high'; dt[logHR_IGX < 0 & grepl(' high$', dt$I_label)]$type <- x; tlvs <- c(tlvs, x)
x <- 'higher-risk IGX; ICL low'; dt[logHR_IGX > 0 & grepl(' low$', dt$I_label)]$type <- x; tlvs <- c(tlvs, x)
x <- 'higher-risk IGX; ICL high'; dt[logHR_IGX > 0 & grepl(' high$', dt$I_label)]$type <- x; tlvs <- c(tlvs, x)
x <- 'lower-risk IGX; ICL low'; dt[logHR_IGX < 0 & grepl(' low$', dt$I_label)]$type <- x; tlvs <- c(tlvs, x)
shapes <- c(24, 25, 21, 22)
names(shapes) <- tlvs
dt$type <- factor(dt$type, levels = tlvs)
stopifnot(all(!is.na(dt$type)))

dt$G <- dt$G_label
dt$G_string_rank <- str_rank(dt$G, numeric = TRUE)
dt$G.type <- sapply(strsplit(dt$G_label, ' '), '[', 2)
G.annot <- dt[, .N, by=c('G', 'G.type', 'G_string_rank')]
setorder(G.annot, -G.type, -N, G_string_rank)
dt$G <- factor(dt$G, levels = G.annot$G)

dt$I <- gsub(' high$', '', gsub(' low$', '', dt$I_label))
dt$I_string_rank <- str_rank(dt$I, numeric = TRUE)
I.annot <- dt[, .N, by=c('I', 'I_string_rank')]
setorder(I.annot, N, -I_string_rank)
dt$I <- factor(dt$I, levels = I.annot$I)

plt <- ggplot()+theme_bw()+
    geom_point(data=dt, aes(x=G, y=I, shape=type, fill=cohort), stroke=0.3, size=3)+
    scale_fill_manual(values = cohort.colors)+
    scale_shape_manual(values = shapes, name = 'IGX annotation')+
    guides(fill = guide_legend(override.aes = list(shape = 22, size=5, stroke=0.3)))+
    labs(fill = 'Cancer type')+
    xlab('\nGenomic features (G)')+
    ylab('Immune features (I)\n')+
    coord_fixed(ratio = 1)+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

pdf(paste0(projDir, '/output/Fig-2A.pdf'), 8, 7)
plot(plt)
invisible(dev.off())
