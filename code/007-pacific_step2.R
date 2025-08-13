# Run PACIFIC-Step2

# check the iterations accumulated so far, then take thier difference with <needed_iters>.
# Runs only if enough iterations for all cohorts are already stored in intermediate/tcga_pacific_runs.

####################################################
# total number of Step1 iterations ultimately needed
needed_iters <- 100000
####################################################

message('loading libraries')
library(data.table)
library(PACIFIC)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

all_datasets <- readRDS(paste0(projDir, '/intermediate/tcga_all_datasets.rds'))


N <- data.table(expand.grid(V1=names(all_datasets), V2=c('mut', 'cna'), V3=0, stringsAsFactors = F))
f <- list.files(paste0(projDir, '/intermediate/tcga_pacific_runs'), pattern = 'iter', recursive = T)
N <- rbind(N, data.table(do.call(rbind, strsplit(f, '[\\/|-]'))[,c(1, 2, 6)]))
N <- N[, list(diff = needed_iters - sum(as.numeric(V3))), by=c('V1', 'V2')][diff > 0]
if(nrow(N) > 0){
    message(paste0('Not enough iterations accumulated so far. Run the previous script as follows:'))
    message(paste(paste0(N$diff, ' iterations for ', N$V1, ' (', N$V2, ')'), collapse = '\n'))
    quit(save = "no", status = 0)
}

for(cohort in names(all_datasets)){
    for(G.type in c('mut', 'cna')){
        message(paste0('\nFor ', cohort, ' ', G.type))
        C <- all_datasets[[cohort]][[G.type]]$C
        output_dir <- paste0(projDir, '/intermediate/tcga_pacific_runs/', cohort, '/', G.type)
        res <- PACIFIC_survival_step2(output_dir = output_dir, anova_baseline = C)
        if(is.null(res)) saveRDS(NULL, paste0(output_dir, '/aggregated-results.rds'))
    }
}



