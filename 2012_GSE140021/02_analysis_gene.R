#\begin{rbatch}
###################################################
### global vars and load libs
###################################################

source('~/scripts/R/common/commonFuns.R')
library(Seurat)
library(ggplot2)
library(sctransform)
mc.cores <- 4
correctBatch <- F
if (correctBatch) {
    file.rdata <- '../data/RData/seurat.genes.RData'
    folder.out <- 'results'
} else {
    file.rdata <- '../data/RData/seurat.genes_noBatchCorr.RData'
    folder.out <- 'results_noBatchCorr'
}
mydircreate(paste0('../', folder.out, '/figs'))
mydircreate(paste0('../', folder.out, '/tables'))
mydircreate(paste0('../', folder.out, '/diffExp'))
mydircreate(paste0('../', folder.out, '/cellTypeSpecific'))

###################################################
### plots
###################################################

load(file.rdata)

#pca and umap by day
pdf(paste0('../', folder.out, '/figs/umap_Group.pdf'), width=1.5 *7)
DimPlot(xg, reduction = "umap", group.by="Group", label=T)
dev.off()
pdf(paste0('../', folder.out, '/figs/pca_Group.pdf'))
DimPlot(xg, reduction = "pca", group.by="Group", label=T)
dev.off()

#umap by batch
sel <- xg[['Batch']]=='batch1'
pdf(paste0('../', folder.out, '/figs/umap_Group_batch1.pdf'), width=1.5 *7)
DimPlot(xg[, sel], reduction = "umap", group.by="Group", label=T)
dev.off()
pdf(paste0('../', folder.out, '/figs/umap_Group_batch2.pdf'), width=1.5 *7)
DimPlot(xg[, !sel], reduction = "umap", group.by="Group", label=T)
dev.off()

#one umap for each group
Group <- as.character(xg$Group)
Idents(xg)=Group
for (Group.s in unique(Group)) {
    xg.s <- subset(x=xg, idents=Group.s)
    xg.s <- RunPCA(xg.s, verbose=FALSE)
    xg.s <- RunUMAP(xg.s, dims=1:30, verbose=FALSE)
    file <- paste0('../', folder.out, '/figs/umap_Group_', Group.s, '.pdf')
    pdf(file, width=1.5 *7)
    plot(DimPlot(xg.s, reduction="umap", group.by="Batch", label=T))
    dev.off()
}

##plot a few genes
sym <- c('NANOS3', 'SOX17', 'NANOG', 'POU5F1', 'KLF4', 'NANOG', 'GATA6', 'SOX17', 'T', 'EOMES', 
         'ZNF611', 'ZNF600', 'ZNF28', 'ZNF468', 'ZNF638',
         'GAPDH', 'TRIM28', 'RPLP0', 'PAX2', 'SOX1', 'ISL1', 'HAND1', 'SOX17', 'PRDM1', 'TFAP2C', 'GATA3', 'TFAP2A', 'CDX2')
ens <- symbol2ensembl(sym, 'Hs')
ens <- ens[ens %in% rownames(xg)]
pdf(paste0('../', folder.out, '/figs/umap_ColorByGenes.pdf'), width=1.5 *7)
for (i in 1:length(ens)) plot(FeaturePlot(xg, features = ens[i], reduction='umap') + ggtitle(names(ens)[i]))
dev.off() 

###################################################
### save norm.counts and do differential expression
###################################################

#load
load(file.rdata)

#preprocess
group <- slot(xg, 'meta.data')[, 'Group']
group <- as.factor(group)
Idents(xg)=group

#get annotation
anno <- annotateEnsBiomart(rownames(xg), 'Hs')
tss <- annotateEnsPos(rownames(xg), 'Hs', strand=F, tss=T, annoWithBiomart=T)[, c(2, 3)]
colnames(tss) <- c('chr', 'tss')
pos <- annotateEnsPos(rownames(xg), 'Hs', strand=T, tss=F, annoWithBiomart=T)[, 3:5]
anno <- data.frame(anno, tss, pos)

#export norm.counts
nc <- as.matrix(GetAssayData(object=xg, slot="data"))
xout <- data.frame(anno, nc)
#file <- paste0('../', folder.out, '/tables/norm.counts.csv')
#write.csv(xout, file=file, row.names=F)

#export average norm.counts
nc.avg <- do.call(cbind, lapply(as.list(levels(group)), function(group.s) rowMeans(nc[, group==group.s, drop=F])))
colnames(nc.avg) <- levels(group)
xout <- data.frame(anno, nc.avg)
write.csv(xout, file=paste0('../', folder.out, '/tables/norm.counts_avg.csv'), row.names=F)

##differential expression
##this will compare each level to the previous one:
#contrast <- lapply(2:nlevels(group), function(i) c(levels(group)[i], levels(group)[i -1])) 
##this will compare each group to a reference one:
#ref <- levels(group)[-1]
#contrast <- lapply(1:length(ref), function(i) c(ref[i], levels(group)[1]))
#for (i in 1:length(contrast)) {
    #cat(paste0('\nContrast ', i, ': ', paste(unlist(contrast[i]), collapse=' vs '), '\n'))
    #difexp <- FindMarkers(xg, ident.1=contrast[[i]][1], ident.2=contrast[[i]][2], min.pct=0)
    #difexp$p_val_adj <- p.adjust(difexp$p_val, 'BH')
    #xout <- data.frame(anno[match(rownames(difexp), anno$ensembl_gene_id), ], difexp)
    #file <- paste0('../', folder.out, '/tables/diffExp/DE_', contrast[[i]][1], '_VS_', contrast[[i]][2], '.csv')
    #write.csv(xout, file=file, row.names=F)
#}

##cell type specific
#for (i in 1:nlevels(group)) {
    #cat(paste0('\nCell type ', i, ': ', levels(group)[i], '\n'))
    #difexp <- FindMarkers(xg, ident.1=levels(group)[i], ident.2=NULL, only.pos=T, min.pct=0)
    #difexp$p_val_adj <- p.adjust(difexp$p_val, 'BH')
    #xout <- data.frame(anno[match(rownames(difexp), anno$ensembl_gene_id), ], difexp)
    #file <- paste0('../', folder.out, '/tables/cellTypeSpecific/CTS_', levels(group)[i], '.csv')
    #write.csv(xout, file=file, row.names=F)
#}

##summary
#for (path in c(paste0('../', folder.out, '/tables/diffExp'), paste0('../', folder.out, '/tables/cellTypeSpecific'))) {
    #files <- dir(path, '^DE_.*.csv|CTS_.*.csv') 
    #if (length(files)>0) {
        #sink(file.path(path, 'summary.txt'))
        #for (i in 1:length(files)) {
            #tmp <- read.csv(file.path(path, files[i]))
            #cnt <- sum(tmp$p_val_adj<.05 & abs(tmp$avg_logFC)>log(2))
            #txt <- paste0(gsub('^DE_|.csv$', '', files[i]), '\n', cnt, '\n\n')
            #cat(txt)
        #}
        #sink()
    #}
#}
#\end{rbatch}
