# write a TSV file containing the expression of LM22 genes for all the METABRIC samples.
# this file is suitable for being uploaded as input to CIBERSORTx (https://cibersortx.stanford.edu).

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

message('loading libraries')
library(data.table)

all_expr <- fread(paste0(projDir, '/input/brca_metabric/data_mrna_illumina_microarray.txt'))

lm22genes <- fread(paste0(projDir, '/input/LM22.txt'))[[1]]

# replace some of the genes with their aliases for compatibility with LM22
x2y <- list(
    c('PGGHG', 'ATHL1'), 
    c('SEPTIN8','SEPT8'),
    c('RADX','CXorf57'),
    c('ADGRE1','EMR1'),
    c('ADGRE3', 'EMR3'),
    c('FCMR','FAIM3'),
    c('GASK1B','FAM198B'),
    c('INKA2','FAM212B'),
    c('RIPOR2','FAM65B'),
    c('KLF3-AS1','FLJ13197'),
    c('H2AC8', 'HIST1H2AE'),
    c('H2BC8','HIST1H2BG'),
    c('RUBCNL','KIAA0226L'),
    c('KIRREL1','KIRREL'),
    c('MARCHF3', 'MARCH3'),
    c('NPIPL2','NPIPB15')
)
lm22expr <- all_expr
for(xy in x2y) lm22expr[Hugo_Symbol == xy[1]]$Hugo_Symbol <- xy[2]
lm22expr <- lm22expr[Hugo_Symbol %in% lm22genes]

dup_genes <- lm22expr[duplicated(Hugo_Symbol)]$Hugo_Symbol
lm22expr <- lm22expr[! Hugo_Symbol %in% dup_genes]
cat('genes present:', nrow(lm22expr), '/', length(lm22genes), '\n')

lm22expr <- data.frame(lm22expr[,-c('Hugo_Symbol','Entrez_Gene_Id')], row.names = lm22expr$Hugo_Symbol)
n <- nrow(lm22expr)
lm22expr <- lm22expr[rowSums(is.na(lm22expr)) == 0,]
cat('rm rows containing NA:', n-nrow(lm22expr), '\n'); rm(n)
cat('dim lm22expr:', dim(lm22expr), '\n')

fwrite(lm22expr, paste0(projDir, '/intermediate/metabric_lm22expr.tsv'), sep='\t', quote = F, row.names = T)
