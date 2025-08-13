# prepare the patient-level clinical features of METABRIC

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

message('loading libraries')
library(data.table)
library(stringr)

# obtain the list of samples in the two subsets of METABRIC, i.e. validation and discovery sets

f <- paste0(projDir, '/input/NIHMS45243-supplement-7/table_S43/IntClustFeatures_validationdataset995.txt')
val_cases <- as.character(read.table(f, nrows = 1))

f <- paste0(projDir, '/input/NIHMS45243-supplement-7/table_S43/IntClustFeatures_discoverydataset997.txt')
dis_cases <- as.character(read.table(f, nrows = 1))


# clean the clinical data of METABRIC

clin_patient <- fread(paste0(projDir, '/input/brca_metabric/data_clinical_patient.txt'))
colnames(clin_patient) <- as.character(clin_patient[4])
clin_patient <- clin_patient[-c(1:4)]

stopifnot(0 == nrow(clin_patient[duplicated(PATIENT_ID)]))

cat('missing val_cases:', length(setdiff(val_cases, clin_patient$PATIENT_ID)), '\n')
cat('missing dis_cases:', length(setdiff(dis_cases, clin_patient$PATIENT_ID)), '\n')

TFoutOfN <- function(b) paste0(length(which(b)), '/', length(b))
cat('val & dis cases, pam50:', TFoutOfN(clin_patient[  PATIENT_ID %in% c(val_cases, dis_cases)]$CLAUDIN_SUBTYPE != ''), '\n')
cat('other cases, pam50:', TFoutOfN(clin_patient[! PATIENT_ID %in% c(val_cases, dis_cases)]$CLAUDIN_SUBTYPE != ''), '\n')

clin_patient <- clin_patient[CLAUDIN_SUBTYPE %in% c('LumA', 'LumB')]

# -----------------------------------------------------------------------------------
clin_patient$RFS_MONTHS <- as.numeric(clin_patient$RFS_MONTHS)
clin_patient$OS_MONTHS <- as.numeric(clin_patient$OS_MONTHS)

stopifnot(setequal(setdiff(clin_patient$RFS_STATUS, NA), c('0:Not Recurred','1:Recurred')))
stopifnot(setequal(setdiff(clin_patient$OS_STATUS, NA), c('0:LIVING','1:DECEASED')))

clin_patient$RFS_STATUS <- ifelse(clin_patient$RFS_STATUS == '0:Not Recurred', 0, 1)
clin_patient$OS_STATUS <- ifelse(clin_patient$OS_STATUS == '0:LIVING', 0, 1)
stopifnot(all(c('Died of Disease', 'Living') %in% clin_patient$VITAL_STATUS))
clin_patient[ ! VITAL_STATUS %in% c('Died of Disease', 'Living')]$OS_STATUS <- NA
# -----------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------
clin_sample <- fread(paste0(projDir, '/input/brca_metabric/data_clinical_sample.txt'))
colnames(clin_sample) <- as.character(clin_sample[4])
clin_sample <- clin_sample[-c(1:4)]
stopifnot(all(clin_sample$SAMPLE_TYPE == 'Primary'))
stopifnot(all(!duplicated(clin_sample$PATIENT_ID)))
# -----------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------
stopifnot(all(clin_patient$PATIENT_ID %in% clin_sample$PATIENT_ID))
stopifnot(identical(clin_sample$PATIENT_ID, clin_sample$SAMPLE_ID))

clin_patient$stage <- clin_sample[match(clin_patient$PATIENT_ID, PATIENT_ID)]$TUMOR_STAGE
clin_patient$grade <- clin_sample[match(clin_patient$PATIENT_ID, PATIENT_ID)]$GRADE

clin_patient[!is.na(stage)]$stage <- paste0('Stage_', clin_patient[!is.na(stage)]$stage)
clin_patient[!is.na(grade)]$grade <- paste0('G', clin_patient[!is.na(grade)]$grade)

clin_patient$stage <- factor(clin_patient$stage, levels = str_sort(unique(clin_patient[!is.na(stage)]$stage), numeric = T))
clin_patient$grade <- factor(clin_patient$grade, levels = str_sort(unique(clin_patient[!is.na(grade)]$grade), numeric = T))

clin_patient$age <- as.numeric(clin_patient$AGE_AT_DIAGNOSIS)
# -----------------------------------------------------------------------------------

clin_patient$case <- ''
clin_patient[PATIENT_ID %in% dis_cases]$case <- 'discovery'
clin_patient[PATIENT_ID %in% val_cases]$case <- 'validation'
stopifnot(all(clin_patient$case != ''))
stopifnot(all(clin_patient$PATIENT_ID %in% c(val_cases, dis_cases))) # equivalent as above

message('saved metabric_clinical_features.rds')
saveRDS(clin_patient, paste0(projDir, '/intermediate/metabric_clinical_features.rds'))
