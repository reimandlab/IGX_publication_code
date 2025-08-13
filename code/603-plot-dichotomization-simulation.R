# this script depends on 10000 iterations of PACIFIC for all cohort for the purpose of Figure S3.

# it cannot readily run as it depends on a separate PACIFIC run.

# but it's a simple procedure: it takes the median values of immune cell values used as 
# threshold for dichotomization over all the iterations and it plots their 
# distribution in the context of each of 34 IGXs

library(data.table)
library(ggplot2)
library(patchwork)
set_plot_size <- function(w, h) options(repr.plot.width = w, repr.plot.height = h)
pth <- function(p){ paste0('/.mounts/labs/reimandlab/private/users/mbayati/prognostic_immunogenomic_markers/', p)}
imm2abb <- Vectorize(function(imm){
    dt <- as.data.table(t(simplify2array(list(
        c('Neutrophils','Neutrophils'),
        c('Monocytes','Monocytes'),
        c('Eosinophils','Eosinophils'),
        c('T_cells_gamma_delta','Gamma Delta T'),
        c('T_cells_CD4_naive','Naive CD4+ T'),
        c('Dendritic_cells_resting','Rest. Dendritic'),
        c('T_cells_regulatory_Tregs','Regulatory T'),
        c('T_cells_follicular_helper','Follicular Helper T'),
        c('Dendritic_cells_activated','Act. Dendritic'),
        c('Macrophages_M0', 'M0 Macrophages'),
        c('Macrophages_M1', 'M1 Macrophages'),
        c('Macrophages_M2', 'M2 Macrophages'),
        c('T_cells_CD4_memory_resting','Rest. Mem. CD4+ T'),
        c('T_cells_CD4_memory_activated','Act. Mem. CD4+ T'),
        c('NK_cells_resting','Rest. NK'),
        c('NK_cells_activated','Act. NK'),
        c('Mast_cells_resting','Rest. Mast'),
        c('T_cells_CD8','CD8+ T'),
        c('Mast_cells_activated','Act. Mast'),
        c('B_cells_memory','Mem. B'),
        c('B_cells_naive','Naive B'),
        c('Plasma_cells','Plasma cells')
    )))) 
    return(dt[V1 == imm]$V2)
})

all_datasets <- readRDS(pth('src_B/current/all_datasets.rds'))
ctlg <- readRDS(pth('src_B/current/catalog.rds'))

res <- list()
for(cn in unique(ctlg$cancer)){
    if(! cn %in% names(res)) res[[cn]] <- list()
    for(at in unique(ctlg[cancer == cn]$alter)){
        f <- paste0(pth('src_B/current/out-with-med-analysis'), '/', cn, '/', at)
        res[[cn]][[at]] <- unlist(lapply(list.files(f, full.names = T), readRDS), recursive = FALSE)
    }
}

dt <- do.call(rbind, apply(ctlg, 1, function(r){
    cbind(id = paste0(gsub('__', '.', r[['cancer']]), ':\n', 
                      gsub('Amp_', '', gsub('Del_', '', r[['alt']])), 
                      ifelse(grepl('^Amp_', r[['alt']]), ' gain', ifelse(grepl('^Del_', r[['alt']]), ' loss', ' mut')), ' &\n',
                      imm2abb(r[['imm']]), ifelse(grepl('lower', r[['imm_lev']]), ' low', ' high')),
          rbind(data.table(type = 'value in cohort', value = 100 * all_datasets[[r[['cancer']]]][[r[['alter']]]][['ft']][[r[['imm']]]]), 
                data.table(type = 'threshold in iter', value = 100 * sapply(lapply(res[[r[['cancer']]]][[r[['alter']]]], '[[', 'med_imms'), 
                                                                            function(m){ m[[r[['imm']]]] }))))
}))
dt2 <- dt[type == 'threshold in iter', list(l = quantile(value, 0.025), m = median(value), u = quantile(value, 0.975)), by='id']
plts <- lapply(dt[type == 'value in cohort', list(m=median(value)), by='id'][order(m, decreasing = T)]$id, function(i){
    ggplot()+theme_bw()+
        geom_violin(data = dt[type == 'value in cohort'][id == i], aes(x=id, y=value), scale = 'width', linewidth=0.25, fill='lightgray')+
        geom_hline(yintercept = median(dt[type == 'value in cohort'][id == i]$value), color='red', linewidth=0.25)+
        geom_rect(data=dt2[id == i], aes(xmin = 0.75, xmax=1.25, ymin=l, ymax=u), color='blue', fill=NA, linewidth=0.25)+
        geom_segment(data=dt2[id == i], aes(x = 0.75, xend=1.25, y=m, yend=m), color='blue', linewidth=0.25)+
        scale_y_continuous(limits = c(-0.02*max(dt[type == 'value in cohort'][id == i]$value), NA),
                           expand = c(0, 0), name = 'Immune cell levels in the cohort')+
        xlab('Cancer : IGX')+
        theme(plot.margin = margin(0,0,0,0),
              axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
              panel.grid = element_blank(), 
              panel.border = element_blank(), 
              axis.line = element_line())
})
pdf('median dichotomization of immune features.pdf', 10, 12)
wrap_plots(plts, nrow = 3, axis_titles = "collect")
dev.off()
