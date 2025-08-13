# load gene-level CNA data in METABRIC.

# subset to genes overlapping with "11q13.1 and 4q35.2 losses" obtained from GISTIC analysis in TCGA-BRCA.

# obtain "11q13.1 loss" and "4q35.2 loss" features as follows:
# the "Genomic" feature (loss event in this case) is determined by whether >75% of genes are affected by the CNA event in each sample.

message('loading libraries')
library(data.table)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

message('loading METABRIC CNA data')
cna <- fread(paste0(projDir, '/input/brca_metabric/data_cna.txt'))
cna <- cna[, -'Entrez_Gene_Id']
dup_genes <- cna[duplicated(Hugo_Symbol)]$Hugo_Symbol
message('removed ', length(dup_genes), ' duplicated genes from cna data: { ', 
        paste(dup_genes, collapse=', '), ' }\n')
cna <- cna[! Hugo_Symbol %in% dup_genes]

cna <- lapply(c('11q13.1', '4q35.2'), function(cytoband){

    # obtain genes overlapping with the loss of this cytoband from GISTIC analysis of TCGA-BRCA
    dr <- paste0(projDir, '/input/tcga_cna_gistic/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0')
    f <- paste0(dr, '/del_genes.conf_99.txt')
    select_genes <- fread(f)
    select_genes <- select_genes[, which(select_genes[1] == cytoband), with=F]
    select_genes <- setdiff(unlist(select_genes[-(1:4)]), '')
    message(paste0('Focus on the ', length(select_genes), 
                   ' genes overlapping with ', cytoband, ' loss obtained from TCGA GISTIC analysis: { '), 
            paste(select_genes, collapse = ', '), ' }\n')

    # subset to the selected genes
    dt <- cna[Hugo_Symbol %in% select_genes]
    missing <- setdiff(select_genes, dt$Hugo_Symbol)
    message(length(missing), ' of them are missing in the METABRIC CNA data: { ', 
            paste(missing, collapse=', '), ' }\n')

    # place genes in row.names
    dt <- data.frame(dt[, -'Hugo_Symbol'], row.names = dt$Hugo_Symbol, check.names = F)

    # make the table binary where loss events (-1, or -2) are set to 1 and the rest are set to 0
    dt[dt > 0] <- 0
    dt[dt < 0] <- 1

    # transpose: patients as rows (listed in the first column) and genes in the rest of columns
    dt <- data.table(t(dt), keep.rownames = 'PATIENT_ID')

    # set the ultimate CNA feature based on the 75% cutoff
    G <- paste0('Del_', cytoband)
    dt[[G]] <- as.numeric(rowMeans(dt[,-'PATIENT_ID']) >= 0.75)
    
    # focus on the two columns: 'patients' and the 'CNA feature'
    return(dt[, c('PATIENT_ID', G), with=F])
})

# merge to have both CNA features in the same data.table
cna <- do.call(merge, cna)

f <- 'metabric_cna_features.rds'
message('saved ', f)
saveRDS(cna, paste0(projDir, '/intermediate/', f))
