# Figure 4:

# Associations of IGXs with genomic and immunogenic characteristics. 
# IGX-positive samples were compared with IGX-negative samples across 
# a panel of genomic and TME characteristics (FDR < 0.05, Mann–Whitney U test).

# High-confidence associations were identified by comparing the 
# IGX-positive group with the other three groups (P < 0.05 from Mann-Whitney U tests


message('loading libraries')
library(data.table)
library(ggplot2)
library(patchwork)
library(readxl)
library(coin)

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------


# load TPM values of 6 genes: CTLA4, LAG3, CD274, PDCD1, GZMA, PRF1, MEN1
tpms <- readRDS(paste0(projDir, '/input/expr_tpm_20_genes_in_tcga.rds'))

# give each gene a title
tpms[[1]] <- paste0(tpms[[1]], ' expression (TPM)')

# reshape the table so the first column is patient id and the other columns are the genes
tpms <- data.frame(tpms[,-1], row.names = tpms[[1]], check.names = F)
tpms <- data.table(t(tpms), keep.rownames = 'patient_id')

# sanity check
stopifnot(0 == sum(is.na(tpms)))

# laod the data from "The Immune Landscape of Cancer"
panimmune_features <- data.table(read_xlsx(paste0(projDir, '/input/NIHMS958212-supplement-2.xlsx'), 
                                           sheet=1, .name_repair='minimal'))

# remove duplicated columns from "panimmune_features". They are not used here anyway.
dup_cols <- names(which(1 != table(colnames(panimmune_features))))
panimmune_features <- panimmune_features[, -dup_cols, with=F]

# merge the two tables
feat.dt <- merge(tpms, panimmune_features, by.x = 'patient_id', by.y = 'TCGA Participant Barcode')

# selected features
features <- c('Leukocyte Fraction',
              'Stromal Fraction',
              'Intratumor Heterogeneity',
              'Proliferation',
              'Wound Healing',
              'Macrophage Regulation',
              'Lymphocyte Infiltration Signature Score',
              'IFN-gamma Response',
              'TGF-beta Response',
              'SNV Neoantigens',
              'Aneuploidy Score',
              'Homologous Recombination Defects',
              'Silent Mutation Rate',
              'Nonsilent Mutation Rate',
              'CTA Score',
              'Th1 Cells',
              'Th2 Cells',
              'Th17 Cells',
              'CTLA4 expression (TPM)',
              'LAG3 expression (TPM)',
              'CD274 expression (TPM)',
              'PDCD1 expression (TPM)',
              'GZMA expression (TPM)',
              'PRF1 expression (TPM)')

stopifnot(all(features %in% colnames(feat.dt)))

# make all of them numeric (this also turns the non-numeric values into NA)
for(f in features) feat.dt[[f]] <- suppressWarnings(as.numeric(feat.dt[[f]]))

# make three new features:
# 'Th1:Th2 cell ratio' = 'Th1 Cells' / 'Th2 Cells'. Still retain those two in 'features'.
new_feat <- 'Th1:Th2 cell ratio'
feat.dt[[new_feat]] <- feat.dt[['Th1 Cells']] / feat.dt[['Th2 Cells']]
features <- c(new_feat, features)
# 'Mutation Rate' is the sum of 'Silent and Nonsilent Mutation Rates'. Then remove those two from 'features'.
new_feat <- 'Mutation Rate'
x <- c('Silent Mutation Rate', 'Nonsilent Mutation Rate')
feat.dt[[new_feat]] <- rowSums(feat.dt[, x, with=F])
features <- c(new_feat, setdiff(features, x))
# 'Cytolytic Activity' = geometric mean of 'GZMA' and 'PRF1' tpm values. Then remove those two from 'features'.
# add a 0.01 to the tpm values prior to applying geometric mean
new_feat <- 'Cytolytic Activity'
x <- c('GZMA expression (TPM)', 'PRF1 expression (TPM)')
feat.dt[[new_feat]] <- sqrt((0.01+feat.dt[[x[1]]])*(0.01+feat.dt[[x[2]]])) # exp(rowMeans(log(0.01 + feat.dt[, x, with=F])))
features <- c(new_feat, setdiff(features, x))

# load the catalog of IGXs
catalog <- readRDS(paste0(projDir, '/intermediate/catalog.rds'))

message('Running tests')

Utest <- function(y, x, dt){
    U.test <- coin::wilcox_test(formula = as.formula(paste(y, '~', x)), data = dt)
    P <- as.numeric(coin::pvalue(U.test)) # P-value of the test
    Z <- U.test@statistic@teststatistic
    N <- nrow(dt)
    R <- Z / sqrt(N) # effect size of the test
    return(list(P=P, R=R))
}

# loop over all feature-IGX pair and test for (binary) association using "Wilcoxon rank sum test" ("Mann-Whitney U test")
res <- do.call(rbind, lapply(features, function(f){
    feat.dt$f <- feat.dt[[f]]
    do.call(rbind, lapply(1:length(catalog), function(ct.indx){
        ct <- catalog[[ct.indx]]
        dt <- ct$patient_level_dt
        dt <- merge(dt[, c('patient_id', 'I', 'G', 'IGX')], feat.dt[, c('patient_id', 'f')], by='patient_id')
        # remove rows with missing data (i.e. NAs in this specific feature)
        dt <- dt[!is.na(f)]
        # to be consistent with coin::wilcox_test(), factorize with positive group being the first level 
        dt$I <- factor(dt$I, levels = c(1, 0))
        dt$G <- factor(dt$G, levels = c(1, 0))
        dt$IGX <- factor(dt$IGX, levels = c(1, 0))
        # make the output structure
        res <- data.table(ct.indx = ct.indx,
                          IGX = paste0(gsub('__', '.', ct$cohort), ': ', ct$G_label, ' & ', ct$I_label), 
                          feature = f, 
                          P = as.numeric(NA),
                          R = as.numeric(NA),
                          P.G = as.numeric(NA),
                          P.I = as.numeric(NA),
                          P.o = as.numeric(NA))
        # Apply test for this pair only if there is enough data - at least 10 samples in each IGX level.
        if(min(table(dt$IGX)) >= 10){
            ut <- Utest('f', 'IGX', dt)
            res$P <- ut$P
            res$R <- ut$R
            if(ut$P <= 0.05){
                res$P.G <- Utest('f', 'IGX', dt[G == 1])$P
                res$P.I <- Utest('f', 'IGX', dt[I == 1])$P
                res$P.o <- Utest('f', 'IGX', dt[G == I])$P
            }
        }
        logHR_IGX <- ct$hr.dt[facet == 'IGX'][v == 'IGX']$logHR
        res$prognosis <- ifelse(logHR_IGX > 0, 'Worse prognosis', 'Better prognosis')
        return(res)
    })) 
}))

# FDR correction over all binary tests
res$FDR <- p.adjust(res$P, method = 'BH')

# sanity check
stopifnot(nrow(res) == length(features) * length(catalog))

message('Make the scatter plot')

# apparently some FDRs were 0, make them equal the least non-zero FDR divided by 100
res[FDR == 0]$FDR <- res[FDR != 0][order(FDR)][1]$FDR/100

# annotate each "IGX" with the {number of significant hits}, {mean of -log10 of significant FDRs}, and {number of performed tests}
igx_annot <- res[, list(N.sig = sum(FDR <= 0.05, na.rm=TRUE), 
                        mean_score = mean(-log10(FDR[FDR <= 0.05]), na.rm=TRUE), 
                        N.tests = sum(!is.na(FDR))), by='IGX']
igx_annot[is.nan(mean_score)]$mean_score <- 0
# sort rows by the metrics in order
setorder(igx_annot, -N.sig, -mean_score, -N.tests, +IGX)
# factorize the "IGX" column with the sorted levels
res$IGX <- factor(res$IGX, levels = igx_annot$IGX)


# annotate each "feature" with the {number of significant hits}, {mean of -log10 of significant FDRs}, and {number of performed tests}
feat_annot <- res[, list(N.sig = sum(FDR <= 0.05, na.rm=TRUE), 
                         mean_score = mean(-log10(FDR[FDR <= 0.05]), na.rm=TRUE),
                         N.tests = sum(!is.na(FDR))), by='feature']
feat_annot[is.nan(mean_score)]$mean_score <- 0
# sort rows by the metrics in order
setorder(feat_annot, +N.sig, +mean_score, +N.tests, -feature)
# factorize the "feature" column with the sorted levels
res$feature <- factor(res$feature, levels = feat_annot$feature)

R_cap <- floor(100*max(res$R, na.rm = T))/100
res$R2 <- res$R
res[R2 > R_cap]$R2 <- R_cap
res[R2 < -R_cap]$R2 <- -R_cap

nlog10FDR_cap <- 5
res$nlog10FDR <- -log10(res$FDR)
res[nlog10FDR > nlog10FDR_cap]$nlog10FDR <- nlog10FDR_cap

size_brks <- seq(2, nlog10FDR_cap)
size_labs <- size_brks
size_labs[length(size_labs)] <- paste0(nlog10FDR_cap, '+')

point.plot <- ggplot()+theme_bw()+
    geom_point(data=res[!is.na(FDR)][FDR <= 0.05], 
               aes(x=IGX, y=feature, fill = R2, size=exp(nlog10FDR), shape = prognosis), 
               stroke=0.3)+
    scale_shape_manual(values=c('Worse prognosis'=24, 'Better prognosis'=25), name = 'IGX association\nwith survival')+
    scale_size(range = c(1.5, 3.6), breaks = exp(size_brks), labels=size_labs, name='Significance,\n-log10(FDR)')+
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', limits=c(-R_cap, R_cap), name = 'Effect size\n(r value)')+
    guides(shape = guide_legend(reverse=T, override.aes = list(size=4)))+
    geom_point(data=res[is.na(FDR)], aes(x=IGX, y=feature), shape=4, stroke=0.3)+
    scale_x_discrete(drop=FALSE)+
    scale_y_discrete(drop=FALSE)+
    coord_fixed(ratio = 1)+
    theme(plot.margin = margin(t=40, l=30, r=10, b=60),
          axis.title = element_blank(),
          panel.grid = element_line(linewidth = 0.3),
          axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))



# --------------------------------------------------------------------------------------------------------
# for the BRCA-LumA IGX specifically. Figure 5H
res2 <- res[IGX == 'BRCA.LumA: 11q13.1 loss & Neutrophils low'][FDR <= 0.05]
res2$feature <- factor(as.character(res2$feature), levels = res2[order(nlog10FDR)]$feature)
sz.brks <- 1:ceiling(max(res2$nlog10FDR))
sz.brks[1] <- 1.3
R_cap <- ceiling(100*max(abs(range(res2$R2))))/100
res2[R2 > R_cap]$R2 <- R_cap
res2[R2 < -R_cap]$R2 <- -R_cap
point.plot.BRCA.LumA <- ggplot()+theme_bw()+
    geom_point(data=res2, 
               aes(x=IGX, y=feature, fill = R2, size=nlog10FDR), 
               shape=21,
               stroke=0.3)+
    scale_size_continuous(limits = range(sz.brks), breaks = sz.brks, name='Significance,\n-log10(FDR)')+
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', limits=c(-R_cap, R_cap), name = 'Effect size\n(r value)')+
    coord_fixed(ratio = 1)+
    theme(panel.grid = element_blank(),
          plot.margin = margin(l=30, r=5, t=5),
          axis.title = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
pdf(paste0(projDir, '/output/Fig-5H.pdf'), 4, 5.6)
point.plot.BRCA.LumA
invisible(dev.off())
# --------------------------------------------------------------------------------------------------------


message('Make the barplots')

# High-confidence associations were identified by comparing the IGX-positive group with the other three groups
barplt_info <- res[FDR <= 0.05][P.G <= 0.05 & P.I <= 0.05 & P.o <= 0.05]

feat2log1p <- c(unique(grep('expression \\(TPM\\)', features, value=T)), 
                       'Cytolytic Activity', 'Mutation Rate', 'SNV Neoantigens',
                       'Homologous Recombination Defects', 'Aneuploidy Score')

stopifnot(all(feat2log1p %in% features))

barplt_info$feature_label <- paste0(barplt_info$feature, ifelse(barplt_info$feature %in% feat2log1p, ' (log1p)', ''))
barplt_info$feature_label <- gsub('\\(TPM\\) \\(log1p\\)', '(log1p TPM)', barplt_info$feature_label)
barplt_info$feature_label <- gsub(' Signature Score', '', barplt_info$feature_label)

barplt_info$IGX_label <- gsub('& ', '&\n', gsub(': ', ':\n', barplt_info$IGX))
barplt_info$label <- paste0(barplt_info$feature_label, '\n', barplt_info$IGX_label)

barplt_data <- do.call(rbind, apply(barplt_info, 1, function(r){
    f <- r['feature']
    feat.dt$f <- feat.dt[[f]]    
    ct <- catalog[[as.numeric(r['ct.indx'])]]
    dt <- ct$patient_level_dt
    dt <- merge(dt[, c('patient_id', 'I', 'G', 'IGX')], feat.dt[, c('patient_id', 'f')], by='patient_id')
    # remove rows with missing data (i.e. NAs in this specific feature)
    dt <- dt[!is.na(f)]
    if(f %in% feat2log1p){ dt$f <- log1p(dt$f) }
    dt <- cbind(feature_label = r['feature_label'], 
                IGX_label = r['IGX_label'],
                label = r['label'], 
                dt)
    group_levels  <- c('Neither', 
                       'Only immune feature', 
                       'Only genomic feature', 
                       'IGX, lower risk association', 
                       'IGX, higher risk association')
    logHR_IGX <- ct$hr.dt[facet == 'IGX'][v == 'IGX']$logHR
    dt$group <- factor(group_levels[1+dt$I+2*dt$G+dt$IGX*(logHR_IGX>0)], 
                       levels = group_levels)
    return(dt)
}))

# ****** this is hardcoding - to a match to the manuscript figure ******  
displayed_features <- c('CTA Score', 'Th2 Cells', 'Lymphocyte Infiltration', 'Th1 Cells', 'TGF-beta Response',
                        paste0(c('PDCD1', 'CTLA4', 'CD274', 'LAG3'), ' expression (log1p TPM)'), 
                        'Aneuploidy Score (log1p)', 'Proliferation', 'Mutation Rate (log1p)')

# it has happened that at this point the each "feature" appears in only one significant association.
stopifnot(identical(sort(barplt_info$feature_label), sort(displayed_features)))

# reorder the rows of "barplt_info" to match to the order in "displayed_features"
barplt_info <- barplt_info[match(displayed_features, feature_label)]

grplvls <- levels(barplt_data$group)
ast <- Vectorize(function(x){
    if(x <= 0.001) return('***')
    if(x <= 0.01) return('**')
    if(x <= 0.05) return('*')
    return('ns')
})
ast_dt <- merge(barplt_data[, list(group = grplvls[1:3], y = max(f)+(1:3)*(max(f)-min(f))/10), by='label'],
                barplt_info[, list(group = grplvls[1:3], sig = ast(c(P.o, P.I, P.G))), by='label'],
                by = c('label', 'group'))
ast_dt$group <- factor(ast_dt$group, levels = grplvls[1:3])

segment_dt <- ast_dt[, list(x=1:3, xend=4, y=y[match(grplvls[1:3], group)]-abs(y[1]-y[2])/5), by='label']

bar.plot.list <- lapply(barplt_info$label, function(lb){
    info <- barplt_info[label == lb]
    dt <- barplt_data[label == lb]
    ggplot()+theme_bw()+
        geom_boxplot(data=dt, aes(x=group, y=f, fill=group), 
                     linewidth=0.3, 
                     outlier.size = 0.5, 
                     outlier.color = 'grey')+
        geom_text(data = ast_dt[label == lb], aes(x=group, y=y, label = sig))+
        geom_segment(data=segment_dt[label == lb], aes(x=x, xend=xend, y=y, yend=y), linewidth=0.3)+
        scale_fill_manual(values = c('white', '#349eeb', '#e7d153', '#74d466', '#eb5466'), 
                          name = 'Tumor samples\ngrouped by IGX status',
                          drop = FALSE)+
        guides(fill = guide_legend(reverse = TRUE))+
        ylab('Score, genomic or immunogenic characteristics')+
        ggtitle(lb)+
        theme(plot.title = element_text(size=6),
              axis.title.x = element_blank(),
              legend.title = element_text(size=10),
              panel.grid = element_blank(), 
              axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(linewidth=0.3))
})

layout <- c(area(1, 1, 1, 1),
            area(1, 2, 1, 2),
            area(1, 3, 1, 3),
            area(1, 4, 1, 4),
            area(1, 5, 1, 5),
            area(2, 1, 2, 1),
            area(2, 2, 2, 2),
            area(2, 3, 2, 3),
            area(2, 4, 2, 4),
            area(3, 1, 3, 1),
            area(3, 2, 3, 2),
            area(3, 3, 3, 3))

bar.plot <- wrap_plots(bar.plot.list, design = layout, guides = 'collect', axis_titles = 'collect')

# plot panels A and B separate pages
message('plotting')
pdf(paste0(projDir, '/output/Fig-4.pdf'), 11, 7.7)
plot(point.plot)
plot(bar.plot)
invisible(dev.off())
