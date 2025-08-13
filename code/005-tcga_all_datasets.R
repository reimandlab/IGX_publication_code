message('loading libraries')
library(data.table)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

clin <- readRDS(paste0(projDir, '/intermediate/tcga_clinical_features.rds'))
imm <- readRDS(paste0(projDir, '/intermediate/tcga_immune_features.rds'))
mut <- readRDS(paste0(projDir, '/intermediate/tcga_mutation_features.rds'))
cna <- readRDS(paste0(projDir, '/intermediate/tcga_cna_features.rds'))

clin <- clin[cohort != 'LAML']

all.dt <- Reduce(function(x, y){ merge(x, y, by='patient_id', all = T) }, list(clin, imm, mut, cna))

# default columns from clinical data
C0 <- c('patient_id', 'cohort', 'endpoint.used', 'time', 'status')

dt <- sapply(unique(clin$cohort), function(this_cohort){
    sapply(c('mut', 'cna'), function(G){
        # subset to this cohort
        dt <- all.dt[cohort == this_cohort]
        # remove columns being entirely NA
        dt <- dt[, nrow(dt) != colSums(is.na(dt)), with=F]
        # obtain available features
        I <- intersect(colnames(dt), colnames(imm[,-'patient_id']))
        C <- intersect(colnames(dt), c('age', 'sex', 'stage', 'grade'))
        if(identical(G, 'mut')) G <- intersect(colnames(dt), colnames(mut[,-'patient_id']))
        if(identical(G, 'cna')) G <- intersect(colnames(dt), colnames(cna[,-'patient_id']))
        # focus on relevant features
        dt <- dt[, c(C0, C, I, G), with=F]
        # now remove rows containing NA
        dt <- dt[0 == rowSums(is.na(dt))]
        # among <C>s that are factor, drop absent levels and remove any with single level
        for(x in C){
            if(class(dt[[x]]) == 'factor'){
                dt[[x]] <- droplevels(dt[[x]])
                if(length(levels(dt[[x]])) == 1){
                    C <- setdiff(C, x)
                    dt <- dt[, -x, with=F]
                }
            }
        }
        # convert <G>s into factors and remove the sparse ones
        for(x in G){
            dt[[x]] <- factor(ifelse(dt[[x]]==0, 'no', 'yes'), c('no', 'yes'))
            if(min(100*table(dt[[x]])/nrow(dt)) < 5){
                G <- setdiff(G, x)
            }
        }
        # include only if the data is not sparse and at least one <G> is available
        if(nrow(dt) >= 100 & 100*nrow(dt[status == 1])/nrow(dt) >= 5 & length(G) > 0){
            # focus on available features
            dt <- dt[, c(C0, C, I, G), with=F]            
            return(list(ft=dt, C=C, I=I, G=G))
        }
    }, simplify = F)
}, simplify = F)

n <- sapply(dt, function(x){ sum(sapply(x, is.null)) })
stopifnot(all(n %in% c(0, 2)))
dt <- dt[names(which(n == 0))]

saveRDS(dt, paste0(projDir, '/intermediate/tcga_all_datasets.rds'))







# # Do not run, this was as internal check

# x <- readRDS('../../src_B/current/all_datasets.rds')
# stopifnot(setequal(names(x), names(dt)))
# stopifnot(all(sapply(names(dt), function(cohort){
#     all(sapply(1:2, function(i){
#         A <- x[[cohort]][[i]]$ft
#         B <- dt[[cohort]][[i]]; setkey(B, NULL)
#         setnames(A, old = c('cancer', 'patient_barcode', 'endpoint'), 
#                  new = c('cohort', 'patient_id', 'endpoint.used'))
#         # *********************************************
#         if(cohort == 'BRCA__Basal') A <- A[,-'sex']
#         colnames(A) <- gsub('_A$', '-A', gsub('_B$', '-B', colnames(A)))
#         for(cl in colnames(A)){
#             if(class(A[[cl]]) == 'factor'){
#                 if('_altered' %in% levels(A[[cl]])){
#                     A[[cl]] <- plyr::mapvalues(A[[cl]], from=levels(A[[cl]]), to=c('no', 'yes'))
#                 }
#             }
#         }
#         # *********************************************
#         stopifnot(all(colnames(B) %in% colnames(A)))
#         A <- A[, colnames(B), with=F]
#         all.equal(B, A)
#     }))
# })))
# t(sapply(names(dt), function(cohort){
#     n1 <- sapply(x[[cohort]], function(y){ ncol(y$ft) })
#     names(n1) <- c('mut', 'cna')
#     n2 <- sapply(dt[[cohort]], ncol)
#     n1-n2
# }))
