# Plot the average values of immune features per each IGX group.

message('loading libraries')
library(data.table)
library(ggplot2)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

cohort.colors <- 
    c('BLCA'='#E9A22E', 
      'BRCA.LumA'='#C25A84', 
      'BRCA.LumB'='#F176A7', 
      'CESC'='#73C5C5', 
      'COAD'='#002CCF', 
      'GBM'='#000000', 
      'HNSC'='#7D2322', 
      'KIRC'='#F9431D', 
      'KIRP'='#EF7727', 
      'LGG'='#6F6F6F', 
      'LIHC'='#145713', 
      'LUAD'='#92905A', 
      'LUSC'='#EBE45C', 
      'OV'='#1F9293', 
      'PAAD'='#B05920', 
      'READ'='#29B376', 
      'SARC'='#B7C98D', 
      'STAD'='#B9ECFD', 
      'UCEC'='#8527B7')

catalog <- readRDS(paste0(projDir, '/intermediate/catalog.rds'))
all_datasets <- readRDS(paste0(projDir, '/intermediate/tcga_all_datasets.rds'))

dt <- do.call(rbind, lapply(1:length(catalog), function(indx){
    
    ct <- catalog[[indx]]
    
    logHR_IGX <- ct$hr.dt[facet == 'IGX'][v == 'IGX']$logHR
    
    ft <- all_datasets[[ct$cohort]][[ct$G.type]]$ft
    I0_group <- ct$patient_level_dt[I == 0]$patient_id
    I1_group <- ct$patient_level_dt[I == 1]$patient_id
    I0_avg <- mean(ft[patient_id %in% I0_group][[ct$I_var]])
    I1_avg <- mean(ft[patient_id %in% I1_group][[ct$I_var]])
    I_dir <- ifelse(grepl(' high$', ct$I_label), 'high', 'low')
    I <- gsub(' high$', '', gsub(' low$', '', ct$I_label))
    if(logHR_IGX < 0){
        temp <- I0_avg
        I0_avg <- I1_avg
        I1_avg <- temp
    }
    data.table(cohort = ct$cohort, I = I,
               low_risk_y = I0_avg, high_risk_y = I1_avg, 
               prog = ifelse(logHR_IGX > 0, '-', '+'))
}))

dt$cohort <- gsub('__', '.', dt$cohort)

I_annot <- dt[, list(.N, v = max(c(low_risk_y, high_risk_y))), by='I']
setorder(I_annot, -N, -v)

# dt$I <- factor(dt$I, levels = I_annot$I)

# hard-code the order of the Immune features in this panel to match to the figure in the publication 
hardcoded_levels <- c('Naive B', 'Act. Mem. CD4+ T', 'CD8+ T', 'Regulatory T', 'Monocytes', 
                      'Rest. Dendritic', 'Neutrophils', 'Mem. B', 'Act. Dendritic', 'Rest. Mem. CD4+ T', 
                      'Gamma Delta T', 'Rest. NK', 'Act. NK', 'M0 Macrophages', 'Eosinophils')
# sanity check
stopifnot(setequal(hardcoded_levels, unique(dt$I)))
# factorize to with the given levels
dt$I <- factor(dt$I, levels = hardcoded_levels)

dt$high_risk_y <- 100*dt$high_risk_y # make percent
dt$low_risk_y <- 100*dt$low_risk_y # make percent


plt <- ggplot()+theme_bw()+
    geom_segment(data=dt, aes(x=0, xend=1, y=low_risk_y, yend=high_risk_y, color=cohort), linewidth=1)+

    geom_point(data=dt[prog == '-'], aes(x=0, y=low_risk_y, shape='neg'), size=1)+
    geom_point(data=dt[prog == '-'], aes(x=1, y=high_risk_y, shape='pos'), size=3)+
    
    geom_point(data=dt[prog == '+'], aes(x=0, y=low_risk_y, shape='pos'), size=3)+
    geom_point(data=dt[prog == '+'], aes(x=1, y=high_risk_y, shape='neg'), size=1)+

    scale_shape_manual(values = c('neg' = 16, 'pos' = 18),
                       labels = c('neg' = 'IGX-negative', 'pos' = 'IGX-positive'), 
                       name = 'Tumor samples')+

    scale_x_continuous(breaks = c(0, 1), labels = c('Lower\nRisk', 'Higher\nRisk'))+
    
    xlab('Tumor samples grouped by IGX status and corresponding risk')+
    ylab('Average ICL in each risk group (%)')+

    facet_wrap(~I, scales = 'free_y', ncol = 5)+
    scale_color_manual(values = cohort.colors, name='Cancer type')+
    theme(panel.grid = element_blank(), strip.background.x = element_blank())

pdf(paste0(projDir, '/output/Fig-3B.pdf'), 10, 5)
plot(plt)
invisible(dev.off())
