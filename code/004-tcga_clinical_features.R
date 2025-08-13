message('loading libraries')
library(data.table)
library(readxl)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

pth <- paste0(projDir, '/input/TCGA-CDR-SupplementalTableS1.xlsx')
clin <- data.table(read_xlsx(pth, sheet = 'TCGA-CDR', guess_max = 5000, .name_repair = 'minimal'))

clin$cohort <- clin$type

# the default endpoint used is OS
clin$endpoint.used <- 'OS'
clin$time <- clin$OS.time
clin$status <- clin$OS

# the endpoint used for certain studies is PFI (instead of OS)
i <- clin$type %in% c('BRCA', 'DLBC', 'LGG', 'PCPG', 'PRAD', 'READ', 'TGCT', 'THCA', 'THYM')
clin[i]$endpoint.used <- 'PFI'
clin[i]$time <- clin[i]$PFI.time
clin[i]$status <- clin[i]$PFI

select_cols <- c('patient_id' = 'bcr_patient_barcode',
                 'study' = 'type',
                 'cohort' = 'cohort',
                 'OS.time' = 'OS.time',
                 'OS.status' = 'OS',
                 'PFI.time' = 'PFI.time',
                 'PFI.status' = 'PFI',
                 'endpoint.used' = 'endpoint.used',
                 'time' = 'time',
                 'status' = 'status',
                 'age' = 'age_at_initial_pathologic_diagnosis',
                 'sex' = 'gender',
                 'stage' = 'ajcc_pathologic_tumor_stage',
                 'grade' = 'histological_grade')
setnames(clin, old = select_cols, new = names(select_cols))
clin <- clin[, names(select_cols), with=F]

# time: 0 --> NA
clin[time == 0]$time <- NA
# time: cap at 10*365 and right censor if beyond
clin[time > 10*365]$status <- 0
clin[time > 10*365]$time <- 10*365

# correct sex levels
clin$sex <- factor(clin$sex, levels = c('FEMALE', 'MALE'))

# correct stage levels
clin$stage <- gsub(' ', '_', clin$stage)
S <- c('IS', paste0('Stage_', c('I', 'II', 'III', 'IV')))
E <- c('A', 'B', 'C')
for(s in S){ for(e in E){ clin[stage == paste0(s, e)]$stage <- s } }
clin$stage <- factor(clin$stage, levels = S)

# correct grade levels
clin$grade <- gsub(' ', '_', clin$grade)
clin$grade <- factor(clin$grade, levels = c('Low_Grade', paste0('G', 1:4), 'High_Grade'))

# split LGG cohort into IDH.mut and IDH.wt subtypes
pth <- paste0(projDir, '/intermediate/tcga_mutation_features.rds')
idh <- readRDS(pth)[, list(IDH = ifelse(IDH1==1 | IDH2==1, 'IDH.mut', 'IDH.wt')), by='patient_id']
clin <- merge(clin, idh, by='patient_id', all = T)
clin_org <- clin[study != 'LGG']
clin_lgg <- clin[study == 'LGG' & !is.na(IDH)]
clin_lgg[, cohort := paste0(cohort, '__', IDH)]
clin <- rbind(clin_org, clin_lgg)[,-'IDH']

# split BRCA cohort into PAM50 subtypes
pth <- paste0(projDir, '/input/PanCancerAtlas_subtypes.rds')
pam50 <- data.table(readRDS(pth))[cancer.type == 'BRCA', c('pan.samplesID', 'Subtype_mRNA')]
colnames(pam50) <- c('patient_id', 'PAM50')
pam50 <- pam50['01' == substr(pam50$patient_id, 14, 15)]
pam50$patient_id <- substr(pam50$patient_id, 1, 12)
stopifnot(all(!duplicated(pam50$patient_id)))
stopifnot(all(!is.na(pam50$PAM50)))
clin <- merge(clin, pam50, by='patient_id', all = T)
clin_org <- clin[study != 'BRCA']
clin_brca <- clin[study == 'BRCA' & !is.na(PAM50)]
clin_brca[, cohort := paste0(cohort, '__', PAM50)]
clin <- rbind(clin_org, clin_brca)[,-'PAM50']

saveRDS(clin, paste0(projDir, '/intermediate/tcga_clinical_features.rds'))
