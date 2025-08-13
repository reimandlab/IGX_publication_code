# make the supplementary table for the catalog of 34 IGXs

message('loading libraries')
library(data.table)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

# load the catalog
catalog <- readRDS(paste0(projDir, '/intermediate/catalog.rds'))

# make the formatted supplementary table for the catalog
dt <- do.call(rbind, lapply(1:length(catalog), function(indx){
    ct <- catalog[[indx]]
    
    logHR_info_IGX <- signif(unlist(ct$hr.dt[facet == 'IGX'][v == 'IGX'][, c('logHR', 'logHR_lower', 'logHR_upper')]), 2)
    logHR_info_G <- signif(unlist(ct$hr.dt[facet == 'G'][v == 'G'][, c('logHR', 'logHR_lower', 'logHR_upper')]), 2)
    logHR_info_I <- signif(unlist(ct$hr.dt[facet == 'I'][v == 'I'][, c('logHR', 'logHR_lower', 'logHR_upper')]), 2)
    
    HR_info_IGX <- signif(unlist(ct$hr.dt[facet == 'IGX'][v == 'IGX'][, c('HR', 'HR_lower', 'HR_upper')]), 2)
    HR_info_G <- signif(unlist(ct$hr.dt[facet == 'G'][v == 'G'][, c('HR', 'HR_lower', 'HR_upper')]), 2)
    HR_info_I <- signif(unlist(ct$hr.dt[facet == 'I'][v == 'I'][, c('HR', 'HR_lower', 'HR_upper')]), 2) 
    
    dt <- ct$patient_level_dt
    
    r <- data.table(`IGX index in supplementary` = indx)

    r[['Cancer type']] <- gsub('__', '-', ct$cohort)
    r[['Genomic feature']] <- ct$G_label
    r[['Immune feature']] <- ct$I_label
    r[['IGX feature selection frequency (%)']] <- signif(ct$EN_score, 2)
    r[['IGX prognosis']] <- ifelse(logHR_info_IGX[1] > 0, 'High risk', 'Low risk')
    
    r[['Number (percent) of samples with IGX feature present']] <- paste0(sum(dt$IGX), ' (', signif(100*sum(dt$IGX)/nrow(dt), 2), ')')
    r[['Number (percent) of samples with genomic feature present']] <- paste0(sum(dt$G), ' (', signif(100*sum(dt$G)/nrow(dt), 2), ')')
    r[['Number (percent) of samples with immune feature present']] <- paste0(sum(dt$I), ' (', signif(100*sum(dt$I)/nrow(dt), 2), ')')
    
    
    
     
    r[['logHR (logHR conf. interval) of IGX feature controlled by baseline clinical variables and both genomic and immune features']] <- 
        paste0(logHR_info_IGX[1], ' (', logHR_info_IGX[2], ', ', logHR_info_IGX[3], ')')
    
    r[['logHR of IGX feature controlled by baseline clinical variables only']] <- 
        signif(ct$logHR_IGX_in_simple_model_controlling_only_for_C, 2)
    
    r[['logHR (logHR conf. interval) of genomic feature controlled by baseline clinical variables']] <-
        paste0(logHR_info_G[1], ' (', logHR_info_G[2], ', ', logHR_info_G[3], ')')
    r[['logHR (logHR conf. interval) of immune feature controlled by baseline clinical variables']] <-
        paste0(logHR_info_I[1], ' (', logHR_info_I[2], ', ', logHR_info_I[3], ')')
    
    
    
    r[['HR (HR conf. interval) of IGX feature controlled by baseline clinical variables and both genomic and immune features']] <- 
        paste0(HR_info_IGX[1], ' (', HR_info_IGX[2], ', ', HR_info_IGX[3], ')')    
    r[['HR (HR conf. interval) of genomic feature controlled by baseline clinical variables']] <-
        paste0(HR_info_G[1], ' (', HR_info_G[2], ', ', HR_info_G[3], ')')
    r[['HR (HR conf. interval) of immune feature controlled by baseline clinical variables']] <-
        paste0(HR_info_I[1], ' (', HR_info_I[2], ', ', HR_info_I[3], ')')
    
    
    
    r[['P-value (adjusted P-value) of ANOVA test for IGX feature controlled by baseline clinical variables and both genomic and immune features']] <- 
        paste0(signif(ct$P_IGX, 2), ' (', signif(ct$FDR_IGX, 2), ')')
    r[['P-value (adjusted P-value) of ANOVA test for IGX feature controlled by baseline clinical variables and genomic feature']] <-
        paste0(signif(ct$P_IGX_other_controls[1], 2), ' (', signif(ct$FDR_IGX_other_controls[1], 2), ')')
    r[['P-value (adjusted P-value) of ANOVA test for IGX feature controlled by baseline clinical variables and immune feature']] <- 
        paste0(signif(ct$P_IGX_other_controls[2], 2), ' (', signif(ct$FDR_IGX_other_controls[2], 2), ')')
    r[['P-value (adjusted P-value) of ANOVA test for IGX feature controlled by baseline clinical variables']] <-
        paste0(signif(ct$P_IGX_other_controls[3], 2), ' (', signif(ct$FDR_IGX_other_controls[3], 2), ')')
    
    r[['P-value (adjusted P-value) of ANOVA test for genomic feature controlled by baseline clinical variables']] <-
        paste0(signif(ct$P_G, 2), ' (', signif(ct$FDR_G, 2), ')')
    r[['P-value (adjusted P-value) of ANOVA test for immune feature controlled by baseline clinical variables']] <-
        paste0(signif(ct$P_I, 2), ' (', signif(ct$FDR_I, 2), ')')
    
    return(r)
}))


# this is a permutation of the rows to match to the supplementary table in the publication
x <- c(1, 2, 3, 6, 12, 7, 11, 9, 10, 4, 8, 5, 13, 16, 15, 14, 17, 18, 19, 20, 21, 22, 23, 24, 26, 25, 28, 27, 30, 29, 31, 34, 33, 32)
# sanity check
stopifnot(setequal(x, 1:nrow(dt)))
# reorder the rows
dt <- dt[x] 
# re-index the rows
dt[[1]] <- 1:nrow(dt)


fwrite(dt, paste0(projDir, '/output/Table-S2.tsv'), sep='\t')
