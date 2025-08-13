# plot the KM-HR panels for the two BRCA IGXs applied to METABRIC

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

# import the functions for plotting
source(paste0(projDir, '/code/010-function_km_plot.R'))

# function to extract P-value from ANOVA of two nested CoxPH models
P_anova_nested_coxph <- function(m1, m2){
    anv <- anova(m1, m2)
    n <- grep('Pr', colnames(anv))
    stopifnot(1 == length(n))
    return(anv[2,n])
}

# Given the required data in the input catalog element, add the results of 
# survival analysis to this element. This is tailored to metabric dataset.
add_surv_results_to_metabric_catalog_element <- function(ct){
    
    ct$endpoint.used <- 'RFS'
    
    ct$C <- c('age', 'grade', 'stage')
    C <- ct$C
    
    ft <- ct$patient_level_dt
    ft$IGX <- ft$G * ft$I
    
    ct$patient_level_dt <- ft

    # CoxPH fit
    CoxPH <- function(x){
        frm <- paste('Surv(time, status) ~', paste(x, collapse = '+'))
        frm <- as.formula(frm)
        suppressWarnings(coxph(frm, data=ft))
    }
    # multivariate CoxPH models
    ct$CoxPH.C.IGX <- CoxPH(c(C, 'IGX'))
    ct$CoxPH.C.G.I.IGX <- CoxPH(c(C, 'I', 'G', 'IGX'))
    ct$CoxPH.C.G.I <- CoxPH(c(C, 'I', 'G'))
    ct$CoxPH.C.G <- CoxPH(c(C, 'G'))
    ct$CoxPH.C.I <- CoxPH(c(C, 'I'))
    ct$CoxPH.C <- CoxPH(C)
    
    # For this specific case, use the simpler anova for P-value of the IGX 
    # that is controlling for clinical covariates only:
    ct$P_IGX <- P_anova_nested_coxph(ct$CoxPH.C.IGX, ct$CoxPH.C) # * * *
    
    ct$P_G <- P_anova_nested_coxph(ct$CoxPH.C.G, ct$CoxPH.C)
    ct$P_I <- P_anova_nested_coxph(ct$CoxPH.C.I, ct$CoxPH.C)
    
    ct$hr.dt <- data.table(facet = c('IGX', 'G', 'I'), 
                           v = c('IGX', 'G', 'I'), 
                           logHR = c(summary(ct$CoxPH.C.IGX)$coefficients['IGX', 'coef'],
                                     summary(ct$CoxPH.C.G)$coefficients['G', 'coef'],
                                     summary(ct$CoxPH.C.I)$coefficients['I', 'coef']))
    
    # Kaplan-Meier estimator for univariate models
    KMFit <- function(v){
        frm <- paste('Surv(time, status) ~', v)
        frm <- as.formula(frm)
        suppressWarnings(survfit(frm, data=ft))
    }
    ct$KMFit.IGX <- KMFit('IGX')
    ct$KMFit.G <- KMFit('G')
    ct$KMFit.I <- KMFit('I')
    return(ct)
}

# load the two metabric datasets
metabric_two_datasets <- readRDS(paste0(projDir, '/intermediate/metabric_two_datasets.rds'))

# make catalog element suitable for KM-HR plot
ct <- list()
ct$cohort <- 'LumA (discovery set)'
ct$patient_level_dt <- metabric_two_datasets[['LumA']]
ct$G_label <- '11q13.1 loss'
ct$I_label <- 'Neutrophils low'
ct <- add_surv_results_to_metabric_catalog_element(ct)
plt1 <- get_km_plots_for_a_catalog_element(ct, 'METABRIC', 'P')

# make catalog element suitable for KM-HR plot
ct <- list()
ct$cohort <- 'LumB (discovery set)'
ct$patient_level_dt <- metabric_two_datasets[['LumB']]
ct$G_label <- '4q35.2 loss'
ct$I_label <- 'Rest. Dendritic low'
ct <- add_surv_results_to_metabric_catalog_element(ct)
plt2 <- get_km_plots_for_a_catalog_element(ct, 'METABRIC', 'P')

pdf(paste0(projDir, '/output/Fig-S6.pdf'), 8.5, 2.9)
plot(plt1)
plot(plt2)
invisible(dev.off())
