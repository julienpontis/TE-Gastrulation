#\begin{rbatch}
###################################################
### global vars and load libs
###################################################

source('~/scripts/R/common/commonFuns.R')
library(Seurat)
library(ggplot2)
library(sctransform)
mc.cores <- 4
mydircreate('../results/qc')

###################################################
### read in, qc and preprocess
###################################################

#read in the data for genes
dummy <- load('../data/RData/countTable.unique.RData')
dummy.n <- load('../data/RData/countTable.unique.names.RData')
rownames(x) <- rownames
colnames(x) <- colnames
xg <-x[grepl('^ENSG', rownames(x)), ]
anno <- annotateEnsPos(rownames(xg), org="Hs")
sel <- !is.na(anno$chr) & anno$chr == "MT"
rownames(xg)[sel] <- paste0("MT-", rownames(xg)[sel])
rm(x); gc()

save(xg, file="../data/RData/genes_rawCounts.RData", compress="gzip")

#read in the data for TEsfam
dummy <- load('../data/RData/countTable.sfam.mmaped.RData')
rownames(x.sf) <- gsub('_', '-', rownames(x.sf))

#merge them and create seurat object
x <- rbind(xg, x.sf)
x <- CreateSeuratObject(x)

#add pheno info
meta <- read.csv('../data/E-MTAB-9388.sdrf.txt', sep='\t')
stopifnot(all(gsub('-', '_', colnames(xg)) %in% meta$Source.Name))
meta <- meta[match(gsub('-', '_', colnames(xg)), meta$Source.Name), ]
samplingSite <- meta$Characteristics.sampling.site. 
authCellType <- meta$Characteristics.authors.inferred.cell.type. 
cellType <- meta$Characteristics.inferred.cell.type.
meta <- data.frame(samplingSite, authCellType, cellType, row.names=colnames(xg))
x <- AddMetaData(x, meta)

#compute mitochondrial percentage
x <- PercentageFeatureSet(x, pattern = "^MT-", col.name = "percent.mt")
percent.mt <- slot(x, 'meta.data')$percent.mt

#visualize QC metrics
pdf('../results/qc/violin_nFeat_ncCount_percMT.pdf', width=2*7)
VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#visualize feature-feature relationships
pdf('../results/qc/dotplot_compare.pdf', width=2*7)
plot1 <- FeatureScatter(x, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

#remove features
#x <- subset(x, subset = nFeature_RNA > 200 & percent.mt < 25)

#run sctransform
x <- SCTransform(x, vars.to.regress="percent.mt", verbose=F, do.correct.umi=F, return.only.var.genes=F)

#separate gene and TEsfam
xg <- x[1:nrow(xg), ]
xsf <- x[-1:-(nrow(xg)), ]

#computeing reductions
xg <- RunPCA(xg, verbose = FALSE)
xg <- RunUMAP(xg, dims=1:30, verbose=FALSE)
xsf <- RunPCA(xsf, verbose = FALSE, approx=FALSE)
xsf <- RunUMAP(xsf, dims=1:10, verbose=T)

#save
save(xg, file='../data/RData/seurat.genes.RData', compress='gzip')
save(xsf, file='../data/RData/seurat.sfam.RData', compress='gzip')
#\end{rbatch}
