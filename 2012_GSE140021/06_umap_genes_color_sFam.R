#\begin{rbatch}
###################################################
### global vars and load libs
###################################################

source('~/scripts/R/common/commonFuns.R')
library(Seurat)
library(ggplot2)
library(sctransform)
correctBatch <- F
if (correctBatch) {
    file.rdata <- '../data/RData/seurat.genes.RData'
    file.rdata.sfam <- '../data/RData/seurat.sfam.RData'
    folder.out <- 'results'
} else {
    file.rdata <- '../data/RData/seurat.genes_noBatchCorr.RData'
    file.rdata.sfam <- '../data/RData/seurat.sfam_noBatchCorr.RData'
    folder.out <- 'results_noBatchCorr'
}
mydircreate(paste0('../', folder.out, '/figs/umap_4each_sfam'))

###################################################
### umap coloring for each subfamily
###################################################

dummyg <- load(file.rdata)
dummyt <- load(file.rdata.sfam)
sfam <- xsf@assays$SCT@scale.data
sfam <- t(sfam[, match(colnames(xg), colnames(sfam))])
sfam[is.na(sfam)] <- 0 #remove this line when all samples are back
rownames(sfam) <- colnames(xg)
colnames(sfam) <- paste0(colnames(sfam), '.norm')
xg <- AddMetaData(xg, sfam, col.name=colnames(sfam))

for (i in 1:ncol(sfam)) {
    sfam.sel <- colnames(sfam)[i]
    file <- paste0('../', folder.out, '/figs/umap_4each_sfam/', sfam.sel, '.pdf')
    pdf(file)
    cols <- c('lightgrey', 'darkblue')
    plot(FeaturePlot(xg, features=gsub('-', '.', sfam.sel), reduction='umap', min.cutoff='q1', max.cutoff='q99', 
                     cols=cols) +ggtitle(paste0(sfam.sel, ' (norm)')))
    dev.off()
}
#\end{rbatch}
