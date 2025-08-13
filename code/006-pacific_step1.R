# Run PACIFIC-Step1

# First argument is mandatory: the number of iterations to run.

# Second and third arguments are optional: the specific cohort and G.type to run for.
# If not passed, the procedure will be called for all cohorts and G.types.

# This script can be called repeatedly for as many times as needed. 
# The outputs are automatically stored in the "intermediate" directory.
# The goal is to accumulate enough number of iterations predefined by <needed_iters>.

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

# ---------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
iters <- suppressWarnings(as.numeric(args[1]))
if(is.na(iters)){
    message('Enter the number of iterations')
    quit(save = "no", status = 0)
}
select_cohorts <- args[2]
if(is.na(select_cohorts)){
    select_cohorts <- names(all_datasets)
} else {
    if(! select_cohorts %in% names(all_datasets)){
        message(paste0('Invalid input cohort: ', select_cohorts))
        quit(save = "no", status = 0)
    }
}
select_G.types <- args[3]
if(is.na(select_G.types)){
    select_G.types <- c('mut', 'cna')
} else {
    if(! select_G.types %in% c('mut', 'cna')){
        message(paste0('Invalid input G.type: ', select_G.types))
        quit(save = "no", status = 0)
    }
}
# ---------------------------------------------

for(cohort in select_cohorts){
    for(G.type in select_G.types){
        out <- paste0(projDir, '/intermediate/tcga_pacific_runs/', cohort, '/', G.type)
        message(paste0('\nFor ', cohort, ' ', G.type))
        # obtain the number of iterations accumulated so far
        f <- list.files(out, pattern = 'iter')
        N <- ifelse(length(f)==0, 0, sum(as.numeric(sapply(strsplit(f, '-'), '[', 4))))
        # if more than needed, do nothing, otherwise, run
        if(N >= needed_iters){
            message(paste0(N, ' iterations have been accumulated already. No need for more.'))
        } else {
            this_iters <- min(iters, needed_iters - N)
            message(paste0('will run ', this_iters, ' iterations'))
            dtst <- all_datasets[[cohort]][[G.type]]
            PACIFIC_survival_step1(data = dtst$ft,
                                   baseline = dtst$C,
                                   feat1 = dtst$I,
                                   feat2 = dtst$G,
                                   num_iterations = this_iters,
                                   output_dir = out)
            message(paste0('Saved intermediate files in ', out))
        }
    }
}