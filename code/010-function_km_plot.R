message('loading libraries')
library(data.table)
library(ggplot2)
library(survival)
library(patchwork)


# function to turn "survfit" object into the data.frame suitable for KM plot

get_km_df <- function(fit, # requires "survfit" object
                      custom_levels, 
                      max_time){
    stopifnot(length(fit$strata) == length(custom_levels))
    df <- rbind(
        do.call(rbind, lapply(names(fit$strata), function(lv){
            strata <- summary(fit)$strata
            time <- summary(fit)$time[strata == lv]
            surv <- summary(fit)$surv[strata == lv]
            n <- sum(strata == lv)
            data.frame(x = c(0, time, time),
                       xend = c(time, max_time, time),
                       y = c(1, surv, 1, surv[-n]),
                       yend = c(1, surv, surv),
                       strata = lv,
                       linewidth = 2)
        })),
        subset(data.frame(x = fit$time,
                          xend = fit$time,
                          y = pmax(0, fit$surv - 0.02),
                          yend = pmin(1, fit$surv + 0.02),
                          strata = rep(names(fit$strata), fit$strata),
                          linewidth = 1,
                          n.censor = fit$n.censor), n.censor > 0, select=-n.censor))
    df$strata <- factor(custom_levels[match(df$strata, names(fit$strata))], levels = custom_levels)
    return(df)
}



# function to make KM and HR plots for an element of the catalog.
# the previous script explains what we mean by "and element of the catalog".

get_km_plots_for_a_catalog_element <- function(ct, consortium='TCGA', sig='P'){ # sig: P or FDR
    
    # make the data.frame for KM plots
    
    max_time <- max(ct$patient_level_dt$time)
    km.df <- rbind(cbind(facet = 'I', get_km_df(ct$KMFit.I, c('Absent', 'G or I present'), max_time)),
                   cbind(facet = 'G', get_km_df(ct$KMFit.G, c('Absent', 'G or I present'), max_time)),
                   cbind(facet = 'IGX', get_km_df(ct$KMFit.IGX, c('Absent', 'IGX present'), max_time)))
    km.df$facet <- factor(km.df$facet, levels = rev(unique(km.df$facet)))
    
    sig_annot <- data.table(facet=c('I', 'G', 'IGX'), sig=c(ct[[paste0(sig, '_I')]], 
                                                            ct[[paste0(sig, '_G')]],    
                                                            ct[[paste0(sig, '_IGX')]]))
    sig_annot$facet <- factor(sig_annot$facet, levels = levels(km.df$facet))
    
    loghr_annot <- data.table(facet=c('I', 'G', 'IGX'), loghr=c(ct$hr.dt[facet == 'I'][v == 'I']$logHR, 
                                                                ct$hr.dt[facet == 'G'][v == 'G']$logHR, 
                                                                ct$hr.dt[facet == 'IGX'][v == 'IGX']$logHR))
    loghr_annot$facet <- factor(loghr_annot$facet, levels = levels(km.df$facet))
    
    n_annot <- data.table(t(apply(ct$patient_level_dt[, c('I', 'G', 'IGX')], 2, table)), keep.rownames = 'facet')
    colnames(n_annot)[-1] <- paste0('n', colnames(n_annot)[-1])
    n_annot$facet <- factor(n_annot$facet, levels = levels(km.df$facet))
    n_annot <- melt(n_annot, 'facet')
    n_annot$x <- 0.75*max_time
    n_annot$y <- 0.92
    n_annot[variable == 'n1']$y <- 1
    n_annot$color <- 'Absent'
    n_annot[variable == 'n1']$color <- 'G or I present'
    n_annot[variable == 'n1'][facet == 'IGX']$color <- 'IGX present'

    # KM plots
    sig.lb <- sig
    km.plt <- ggplot()+theme_bw()+
        geom_segment(data = km.df, aes(x=x, xend=xend, y=y, yend=yend, color=strata, linewidth=linewidth))+
        scale_linewidth_continuous(range = c(0.3, 0.4), guide = 'none')+
        geom_text(data=sig_annot, aes(x=0, y=0.13, label = paste0(sig.lb, ' = ', signif(sig, 2))), hjust=0, size=3)+
        geom_text(data=loghr_annot, aes(x=0, y=0.05, label = paste0('logHR = ', signif(loghr, 2))), hjust=0, size=3)+
        geom_text(data=n_annot, aes(x=x, y=y, label=paste0('n = ', value), color=color), size=3, hjust=0, show.legend = FALSE)+
        scale_y_continuous(limits = c(0, 1))+
        scale_x_continuous(limits = c(0, max_time))+
        facet_wrap(~facet)+
        scale_color_manual(values = c('#8C8C8C', 'black', 'red'), name=NULL)+
        guides(color = guide_legend(reverse = TRUE))+
        ylab('Survival Probability')+
        xlab('Time')+
        ggtitle(paste0(consortium, ' - ', 
                       gsub('__', ' ', ct$cohort), ': ', 
                       ct$G_label, ' & ', 
                       ct$I_label, 
                       ' (', ct$endpoint.used,')'))+
        theme(panel.grid = element_blank(),
              strip.text.x = element_blank(),
              strip.background.x = element_blank())
    
    
    
    return(km.plt)
}
