#\begin{rbatch}
###################################################
### global vars and load libs
###################################################

source('~/scripts/R/common/commonFuns.R')
library(Seurat)
library(ggplot2)
library(sctransform)
mc.cores <- 4
mydircreate('../results/figs')
mydircreate('../results/tables/diffExp')
mydircreate('../results/tables/cellTypeSpecific')
mydircreate('../results/tables_subClusters/diffExp')
mydircreate('../results/tables_subClusters/cellTypeSpecific')

###################################################
### plots
###################################################

load('../data/RData/seurat.genes.RData')

#pca and umap by day
pdf('../results/figs/umap_samplingSite.pdf', width=1.5 *7)
DimPlot(xg, reduction = "umap", group.by="samplingSite", label=T)
dev.off()
pdf('../results/figs/pca_samplingSite.pdf')
DimPlot(xg, reduction = "pca", group.by="samplingSite", label=T)
dev.off()

pdf('../results/figs/umap_authCellType.pdf', width=1.5 *7)
DimPlot(xg, reduction = "umap", group.by="authCellType", label=T)
dev.off()
pdf('../results/figs/pca_authCellType.pdf')
DimPlot(xg, reduction = "pca", group.by="authCellType", label=T)
dev.off()

pdf('../results/figs/umap_cellType.pdf', width=1.5 *7)
DimPlot(xg, reduction = "umap", group.by="cellType", label=T)
dev.off()
pdf('../results/figs/pca_cellType.pdf')
DimPlot(xg, reduction = "pca", group.by="cellType", label=T)
dev.off()

#plot a few genes
ens <- symbol2ensembl(c('NANOS3', 'GATA6', 'SOX17', 'CTCF', 'TRIM28', 'GAPDH'), 'Hs')
pdf("../results/figs/umap_ColorByGenes.pdf", width=1.5 *7)
for (i in 1:length(ens)) plot(FeaturePlot(xg, features = ens[i], reduction='umap') + ggtitle(names(ens)[i]))
dev.off() 

###################################################
### umap with new groups
###################################################

load('../data/RData/seurat.genes.RData')
newmeta <- readRDS('../data/RData/annot_umap-1.rds')

newmeta <- newmeta[match(colnames(xg), gsub('\\.', '-', newmeta$cell_name)), ]
meta <- data.frame(subCluster=newmeta$sub_cluster, row.names=colnames(xg))
xg <- AddMetaData(xg, meta)

pdf('../results/figs/umap_authCellType_subCluster.pdf', width=1.5 *7)
DimPlot(xg, reduction = "umap", group.by="subCluster", label=T)
dev.off()
