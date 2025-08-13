# plot KM-HR panels for all IGXs in the catalog

# ---------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
scriptPth <- normalizePath(gsub('--file=', '', grep('--file=', args, value=T)))
projDir <- dirname(dirname(scriptPth))
# ---------------------------------------------

# import the functions for plotting
source(paste0(projDir, '/code/010-function_km_plot.R'))

# load the catalog
catalog <- readRDS(paste0(projDir, '/intermediate/catalog.rds'))

message('plotting')
pdf(paste0(projDir, '/output/Fig-S2.pdf'), 8.5, 2.9)
for(i in 1:length(catalog)){
    cat(paste0('\rProgress: ', floor(100*i/length(catalog)), '%    ')); flush.console()
    p <- get_km_plots_for_a_catalog_element(catalog[[i]])
    plot(p)
}
invisible(dev.off())
cat('\n')
