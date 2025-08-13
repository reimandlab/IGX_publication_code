# Make two datasets from METABRIC:
# 1) for analysis of the IGX "11q13.1 loss & Neutrophils low" in LumA subtype (subset to discovery set of metabric)
# 2) for analysis of the IGX "4q35.2 loss & Rest. Dendritic low" in LumB subtype (subset to discovery set of metabric)


# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

message('loading libraries')
library(data.table)

# load the clinical, cna, and immune features built for METABRIC
clin <- readRDS(paste0(projDir, '/intermediate/metabric_clinical_features.rds'))
cna <- readRDS(paste0(projDir, '/intermediate/metabric_cna_features.rds'))
imm <- readRDS(paste0(projDir, '/intermediate/metabric_immune_features.rds'))

# focus on the "discovery" subset of METABRIC.
clin <- clin[case == 'discovery']

# select RFS (relapse-free survival) for type survival data
clin$endpoint.used <- 'RFS'
clin$time <- clin$RFS_MONTHS
clin$status <- clin$RFS_STATUS

# check time / status
stopifnot(all(!is.na(clin$time)) & all(!is.na(clin$status)))
stopifnot(all(clin$time > 0) & all(clin$status %in% c(0, 1)))

# right-censor any row with time > 10 years (10x12 months)
mxtm <- 10 * 12
indx <- clin$time > mxtm
if (length(indx) > 0) {
    clin[indx]$time <- mxtm
    clin[indx]$status <- 0
}

message('making the two datasets for LumA and LumB')

dt <- sapply(c('LumA', 'LumB'), function(subtype){
    
    # IGX "11q13.1 loss & Neutrophils low" in LumA
    if(subtype == 'LumA'){
        cna.event <- 'Del_11q13.1'
        imm.cell <- 'Neutrophils'
    }
    
    # IGX "4q35.2 loss & Rest. Dendritic low" in LumB
    if(subtype == 'LumB'){
        cna.event <- 'Del_4q35.2'
        imm.cell <- 'Dendritic_cells_resting'
    }
    
    # focus on this subtype
    clin <- clin[CLAUDIN_SUBTYPE == subtype]
    
    # set the "genomic" feature from cna
    cna$G <- cna[[cna.event]]
    
    # the immune feature in this case is defined by the Immune cell level being "less than or equal the median in the cohort".
    imm$I <- as.numeric(imm[[imm.cell]] <= median(imm[[imm.cell]]))

    # the baesline covariates
    C <- c('age','grade','stage')

    # merge clinical, G, and I features
    dt <- Reduce(function(x, y){ merge(x, y, by = 'PATIENT_ID') }, 
                 list(clin[, c('PATIENT_ID', 'endpoint.used', 'time', 'status', C), with=F], 
                      cna[, c('PATIENT_ID', 'G')], 
                      imm[, c('PATIENT_ID', imm.cell, 'I'), with=F]))

    # remove rows with NA stage or grade and drop their missing levels
    dt <- dt[!is.na(stage) & !is.na(grade)]
    dt$stage <- droplevels(dt$stage)
    dt$grade <- droplevels(dt$grade)

    # make the levels of the "stage" variable consistent with TCGA
    lvs <- levels(dt$stage)
    dt$stage <- as.character(dt$stage)
    dt$stage <- gsub('_1', '_I', gsub('_2', '_II', gsub('_3', '_III', gsub('_4', '_IV', dt$stage))))
    lvs <- gsub('_1', '_I', gsub('_2', '_II', gsub('_3', '_III', gsub('_4', '_IV', lvs))))
    dt$stage <- factor(dt$stage, levels = lvs)
    # check if levels of the "grade" variable is already consistent with TCGA
    stopifnot(all(dt$grade %in% paste0('G', 1:4)))
    
    # no remaining NA
    stopifnot(0 == sum(is.na(dt)))
    
    return(dt)
}, simplify=FALSE, USE.NAMES=TRUE)

f <- 'metabric_two_datasets.rds'
message('saved ', f)
saveRDS(dt, paste0(projDir, '/intermediate/', f))
