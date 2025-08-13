# make the catalog of IGXs based on PACIFIC results (step 1 & 2).

# the catalog object is a two-layer nested list: 
# -- in the first layer it points to each IGX in the catalog.
# -- in the second layer it holds, in list format, detailed information about the corresponding IGX, including:
# ----- cohort name and the end point used.
# ----- G/I components and their labels.
# ----- Feature selection score of the IGX.
# ----- logHR, P-value, and FDR associated with G, I, and IGX feeatures.
# ----- Objects of CoxPH and KM-fit models associated with the IGX, G, and I features.
# ----- annotation of all patients of the cohort with the G and I features, and clinical data.

message('loading libraries')
library(data.table)
library(survival)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

# function to format the names of the immune cells (I features)
format_I <- function(I){
    imm_names <- c('Neutrophils' = 'Neutrophils',
                   'Monocytes' = 'Monocytes',
                   'Eosinophils' = 'Eosinophils',
                   'T_cells_gamma_delta' = 'Gamma Delta T',
                   'T_cells_CD4_naive' = 'Naive CD4+ T',
                   'Dendritic_cells_resting' = 'Rest. Dendritic',
                   'T_cells_regulatory_Tregs' = 'Regulatory T',
                   'T_cells_follicular_helper' = 'Follicular Helper T',
                   'Dendritic_cells_activated' = 'Act. Dendritic',
                   'Macrophages_M0' = 'M0 Macrophages',
                   'Macrophages_M1' = 'M1 Macrophages',
                   'Macrophages_M2' = 'M2 Macrophages',
                   'T_cells_CD4_memory_resting' = 'Rest. Mem. CD4+ T',
                   'T_cells_CD4_memory_activated' = 'Act. Mem. CD4+ T',
                   'NK_cells_resting' = 'Rest. NK',
                   'NK_cells_activated' = 'Act. NK',
                   'Mast_cells_resting' = 'Rest. Mast',
                   'T_cells_CD8' = 'CD8+ T',
                   'Mast_cells_activated' = 'Act. Mast',
                   'B_cells_memory' = 'Mem. B',
                   'B_cells_naive' = 'Naive B',
                   'Plasma_cells' = 'Plasma cells')
    return(imm_names[I])
}

# function to format the names of the genomic events (G features)
format_G <- function(G){
    paste0(gsub('Amp_', '', gsub('Del_', '', G)), 
           ifelse(grepl('Amp_', G), ' gain', ifelse(grepl('Del_', G), ' loss', ' mut')))
}

# function to extract P-value from ANOVA of two nested CoxPH models
P_anova_nested_coxph <- function(m1, m2){
    anv <- anova(m1, m2)
    n <- grep('Pr', colnames(anv))
    stopifnot(1 == length(n))
    return(anv[2,n])
}

message('merge the aggregated results from all cohorts')

# load all datasets
all_datasets <- readRDS(paste0(projDir, '/intermediate/tcga_all_datasets.rds'))

tb <- list() # the merged table
i <- 0; n <- 2*length(all_datasets)
for(cohort in names(all_datasets)){
    for(G.type in c('mut', 'cna')){
        i <- i+1; cat(paste0('\rProgress: ', floor(100*i/n), '%    ')); flush.console()
        out_dir <- paste0(projDir, '/intermediate/tcga_pacific_runs/', cohort, '/', G.type)
        res <- readRDS(paste0(out_dir, '/aggregated-results.rds'))
        if(!is.null(res)) tb <- c(tb, list(cbind(cohort=cohort, G.type=G.type, res$top_interactions)))
    }
}
cat('\n')
tb <- as.data.table(do.call(rbind, tb))
message(nrow(tb), ' candidate IGXs (elastic net selection frequency > 50%)')

# apply FDR correction to each P-value column. 
igx_fdrs <- list()
for(cl in grep('_P', colnames(tb), value=T)){
    fdr <- p.adjust(tb[[cl]], method = 'BH')  
    tb[[gsub('_P', '_FDR', cl)]] <- fdr
    if(grepl('^intr_', cl)) igx_fdrs <- c(igx_fdrs, list(fdr))
}

# keep only rows where all FDRs associated with IGXs are significant (<= 0.05)
tb <- tb[0 == rowSums(do.call(cbind, igx_fdrs) > 0.05)]
message(nrow(tb), ' IGXs in final catalog (FDR < 0.05 for the 4 IGX anova P values)')

# order rows by cohort and feature selection frequency
setorder(tb, cohort, -EN_score)

message('make the catalog')

# turn the table into a nested list
catalog <- lapply(split(tb, seq(nrow(tb))), as.list)

# for each element of catalog, add furher information
for(indx in 1:length(catalog)){    
    cat(paste0('\rProgress: ', floor(100*indx/length(catalog)), '%    ')); flush.console()
    
    ct <- catalog[[indx]]
    
    # load the dataset specific for this IGX
    dtst <- all_datasets[[ct$cohort]][[ct$G.type]]
    ft <- dtst$ft
    C <- dtst$C
    I <- strsplit(ct$intr, '\\*')[[1]][1]
    G <- strsplit(ct$intr, '\\*')[[1]][2]
    ft$I <- as.numeric(ft[[I]] > median(ft[[I]]))
    if(grepl('lower', ct$feat1_level)) ft$I <- 1-ft$I
    ft$G <- as.numeric(ft[[G]])-1
    ft$IGX <- ft$I * ft$G
    
    # ----------------------------------------------------------------------------------------------------
    # pick P values and FDRs from PACIFIC restuls
    ct$P_IGX <- ct$intr_P_C3 # igx vs G+I+C
    ct$P_IGX_other_controls <- c(ct$intr_P_C2, ct$intr_P_C1, ct$intr_P) # igx vs G+C; igx vs I+C; igx vs C
    ct$P_G <- ct$feat2_P
    ct$P_I <- ct$feat1_P
    ct$FDR_IGX <- ct$intr_FDR_C3 # igx vs G+I+C
    ct$FDR_IGX_other_controls <- c(ct$intr_FDR_C2, ct$intr_FDR_C1, ct$intr_FDR) # igx vs G+C; igx vs I+C; igx vs C
    ct$FDR_G <- ct$feat2_FDR
    ct$FDR_I <- ct$feat1_FDR
    
    # CoxPH fit
    CoxPH <- function(x){
        frm <- paste('Surv(time, status) ~', paste(x, collapse = '+'))
        frm <- as.formula(frm)
        suppressWarnings(coxph(frm, data=ft))
    }
    
    # multivariate CoxPH models
    ct$CoxPH.C.G.I.IGX <- CoxPH(c(C, 'I', 'G', 'IGX'))
    ct$CoxPH.C.G.IGX <- CoxPH(c(C, 'G', 'IGX'))
    ct$CoxPH.C.I.IGX <- CoxPH(c(C, 'I', 'IGX'))
    ct$CoxPH.C.IGX <- CoxPH(c(C, 'IGX'))
    ct$CoxPH.C.G.I <- CoxPH(c(C, 'I', 'G'))
    ct$CoxPH.C.G <- CoxPH(c(C, 'G'))
    ct$CoxPH.C.I <- CoxPH(c(C, 'I'))
    ct$CoxPH.C <- CoxPH(C)
    
    # sanity check: re-calculate the P-values and check if consistent with those taken from PACIFIC results
    stopifnot(ct$P_IGX == P_anova_nested_coxph(ct$CoxPH.C.G.I.IGX, ct$CoxPH.C.G.I))
    stopifnot(identical(ct$P_IGX_other_controls, c(P_anova_nested_coxph(ct$CoxPH.C.G.IGX, ct$CoxPH.C.G), 
                                                   P_anova_nested_coxph(ct$CoxPH.C.I.IGX, ct$CoxPH.C.I),
                                                   P_anova_nested_coxph(ct$CoxPH.C.IGX, ct$CoxPH.C))))
    stopifnot(ct$P_G == P_anova_nested_coxph(ct$CoxPH.C.G, ct$CoxPH.C))
    stopifnot(ct$P_I == P_anova_nested_coxph(ct$CoxPH.C.I, ct$CoxPH.C))
    # ----------------------------------------------------------------------------------------------------
    
    
    
    # ----------------------------------------------------------------------------------------------------
    # extract HR values and their confidence intervals for all variables
    hr.dt <- rbind(cbind(facet = 'IGX', data.table(do.call(cbind, summary(ct$CoxPH.C.G.I.IGX)[c('coefficients', 'conf.int')])[, c('exp(coef)', 'lower .95', 'upper .95', 'coef')], keep.rownames = 'v')),
                   cbind(facet = 'G', data.table(do.call(cbind, summary(ct$CoxPH.C.G)[c('coefficients', 'conf.int')])[, c('exp(coef)', 'lower .95', 'upper .95', 'coef')], keep.rownames = 'v')),
                   cbind(facet = 'I', data.table(do.call(cbind, summary(ct$CoxPH.C.I)[c('coefficients', 'conf.int')])[, c('exp(coef)', 'lower .95', 'upper .95', 'coef')], keep.rownames = 'v')))
    hr.dt$v <- gsub('^grade', 'Grade ', 
                    gsub('_Grade', '', 
                         gsub('^stage', '', 
                              gsub('Stage_', 'Stage ', 
                                   gsub('sexMALE', 'Sex M', 
                                        gsub('^age$', 'Age', hr.dt$v))))))
    colnames(hr.dt)[-(1:2)] <- c('HR', 'HR_lower', 'HR_upper', 'logHR')
    hr.dt$logHR_lower <- log(hr.dt$HR_lower)
    hr.dt$logHR_upper <- log(hr.dt$HR_upper)
    ct$hr.dt <- hr.dt
    # ----------------------------------------------------------------------------------------------------
    
    # ----------------------------------------------------------------------------------------------------
    ct$logHR_IGX_in_simple_model_controlling_only_for_C <- summary(ct$CoxPH.C.IGX)[['coefficients']]['IGX', 'coef']
    # ----------------------------------------------------------------------------------------------------
    
    # Kaplan-Meier estimator for univariate models
    KMFit <- function(v){
        frm <- paste('Surv(time, status) ~', v)
        frm <- as.formula(frm)
        suppressWarnings(survfit(frm, data=ft))
    }
    ct$KMFit.IGX <- KMFit('IGX')
    ct$KMFit.G <- KMFit('G')
    ct$KMFit.I <- KMFit('I')
    
    # additional info
    ct$G_var <- G
    ct$G_var_level <- ct$feat2_level
    ct$G_label <- format_G(G)
    
    ct$I_var <- I
    ct$I_var_level <- ct$feat1_level
    ct$I_label <- paste0(format_I(I), ifelse(grepl('lower', ct$feat1_level), ' low', ' high'))
    
    endpoint <- unique(ft$endpoint.used)
    stopifnot(1 == length(endpoint))
    ct$endpoint.used <- endpoint
    ct$patient_level_dt <- ft[, c('patient_id', 'time', 'status', C, 'I', 'G', 'IGX'), with=F]
    
    
    
    # pick the final facets
    ct <- ct[c(
    'cohort',
    'G.type',
    'EN_score',
    'endpoint.used',
    'patient_level_dt',
    'P_IGX',
    'P_IGX_other_controls',
    'P_G',
    'P_I',
    'FDR_IGX',
    'FDR_IGX_other_controls',
    'FDR_G',
    'FDR_I',
    'CoxPH.C.G.I.IGX',
    'CoxPH.C.G.IGX',
    'CoxPH.C.I.IGX',
    'CoxPH.C.IGX',
    'CoxPH.C.G.I',
    'CoxPH.C.G',
    'CoxPH.C.I',
    'CoxPH.C',
    'hr.dt',
    'logHR_IGX_in_simple_model_controlling_only_for_C',
    'KMFit.IGX',
    'KMFit.G',
    'KMFit.I',
    'G_var',
    'G_var_level',
    'G_label',
    'I_var',
    'I_var_level',
    'I_label')]
    
    # update the catalog
    catalog[[indx]] <- ct
}

cat('\n')

message('save the catalog')
saveRDS(catalog, paste0(projDir, '/intermediate/catalog.rds'))

