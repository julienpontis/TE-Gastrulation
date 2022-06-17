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
correctBatch <- F

###################################################
### read in, qc and preprocess
###################################################

#read in the data for genes
dummy <- load('../data/RData/countTable.unique.RData')

#remove samples
col.rm <- which(grepl('T2C|_hE_|_U1_', colnames(x)))
x <- x[, -col.rm]

#add MT to mithochondrial gene
ens <- rownames(x)[grepl('^ENSG', rownames(x))]
anno <- annotateEnsPos(ens, org="Hs")
sel <- !is.na(anno$chr) & anno$chr == "MT"
rownames(x)[sel] <- paste0("MT-", rownames(x)[sel])

#remove underscore from TEsfam names
rownames(x) <- gsub('_', '-', rownames(x))

#create seurat object
x <- CreateSeuratObject(x)

#add pheno info
snames <- gsub('batch1_|batch2_|U2_', '', colnames(x))
snames <- gsub('hPGCLC', 'PGCLC', snames)
Group <- unlist(lapply(strsplit(snames, '_'), function(x) paste(x[-length(x)], collapse='_')))
Batch <- unlist(lapply(strsplit(colnames(x), '_'), function(x) x[1]))
Group <- as.factor(Group)
Batch <- as.factor(Batch)
meta <- data.frame(Group, Batch, row.names=colnames(x))
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
x <- subset(x, subset = nFeature_RNA < 8000 & percent.mt < 20)

#run sctransform
#load('~/Desktop/x.bef.sct.RData')
if (correctBatch) {
    vars.to.regress=c('percent.mt', 'Batch')
} else {
    vars.to.regress='percent.mt'
}
#sel <- rowSums(as.matrix(GetAssayData(object=x, slot="data"))>0)
x <- SCTransform(x, vars.to.regress=vars.to.regress, verbose=F, do.correct.umi=F, return.only.var.genes=F, 
                 conserve.memory=T, variable.features.n=Inf)

#separate gene and TEsfam
selg <- grepl('ENSG', rownames(x)) #important to look for ENSG anywhere and not only beginning to also get MT-ENSG
xg <- x[selg, ]
xsf <- x[!selg, ]

#computeing reductions
xg <- RunPCA(xg, verbose = FALSE)
xg <- RunUMAP(xg, dims=1:30, verbose=FALSE)
xsf <- RunPCA(xsf, verbose = FALSE)
xsf <- RunUMAP(xsf, dims=1:30, verbose=FALSE)

#save
if (correctBatch) {
    save(xg, file='../data/RData/seurat.genes.RData', compress='gzip')
    save(xsf, file='../data/RData/seurat.sfam.RData', compress='gzip')
} else {
    save(xg, file='../data/RData/seurat.genes_noBatchCorr.RData', compress='gzip')
    save(xsf, file='../data/RData/seurat.sfam_noBatchCorr.RData', compress='gzip')
}
#\end{rbatch}
