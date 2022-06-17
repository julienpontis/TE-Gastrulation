#\begin{rbatch}
###################################################
### global vars and load libs
###################################################

setwd('~/samba/projects/jpontis/2008_scrnaseq_E-MTAB-9388/scripts')
source('~/scripts/R/common/commonFuns.R')
library(Seurat)
library(ggplot2)
library(sctransform)
mc.cores <- 4
mydircreate('../results_sfam/tables')
mydircreate('../results_sfam/tables/diffExp')
mydircreate('../results_sfam/tables/cellTypeSpecific')
mydircreate('../results_sfam/tables_subClusters')
mydircreate('../results_sfam/tables_subClusters/diffExp')
mydircreate('../results_sfam/tables_subClusters/cellTypeSpecific')
mydircreate('../results_sfam/qc')
mydircreate('../results_sfam/figs/dotplots')
mydircreate('../results_sfam/figs_subClusters/dotplots/')
mydircreate('../results_sfam/figs_subClusters/boxplots/')

###################################################
### boxplots per sfamily
###################################################

dummy <- load('../data/RData/seurat.sfam.RData')
dat <- as.matrix(GetAssayData(object=xsf, slot="scale.data"))
group <- as.factor(slot(xsf, 'meta.data')[, 'authCellType'])

for (i in 1:nrow(dat)) {
    x <- dat[i, ]
    xl <- number2list(x, group)
    fam.lbl <- rownames(dat)[i]
    p <- dotplot4list(xl, ylab='Expression (scaled)', main=fam.lbl, jitter.width=.2, col.point=rgb(0, 0, .2, .2),
                                            addN=T, returnObject=T)
    p <- p + coord_flip()
    pdf(paste0('../results_sfam/figs/dotplots/', fam.lbl, '.pdf'))
    plot(p)
    dev.off()
}

###################################################
### boxplots per sfamily with subclusters
###################################################

#load
load('../data/RData/seurat.sfam.RData')

#add subclusters
newmeta <- readRDS('../data/RData/annot_umap-1.rds')
newmeta <- newmeta[match(colnames(xsf), gsub('\\.', '-', newmeta$cell_name)), ]
subCluster <- gsub(' ', '_', newmeta$sub_cluster)
meta <- data.frame(subCluster=subCluster, row.names=colnames(xsf))
xsf <- AddMetaData(xsf, meta)

#preprocess
dat <- as.matrix(GetAssayData(object=xsf, slot="scale.data"))
group <- slot(xsf, 'meta.data')[, 'subCluster']
group <- as.factor(group)

#order levels
levels <- c('Epiblast', 'Primitive_Streak', 'PGC', 'Nascent_Mesoderm', 'Emergent_Mesoderm', 'Advanced_Mesoderm', 'Axial_Mesoderm', 
            'YS_Mesoderm', 'DE(P)', 'DE(NP)', 'Hypoblast', 'YS_Endoderm', 'Ectoderm', 'Hemogenic_Endothelium', 'Myeloid_Progenitors', 
            'Erythro-Myeloid_Progenitors', 'Blood_Progenitors', 'Erythrocytes')
levels <- levels[length(levels):1]
group <- factor(group, levels=levels)

for (i in 1:nrow(dat)) {
    x <- dat[i, ]
    xl <- number2list(x, group)
    fam.lbl <- rownames(dat)[i]
    p <- dotplot4list(xl, ylab='Expression (scaled)', main=fam.lbl, jitter.width=.2, col.point=rgb(0, 0, .2, .2),
                                            addN=T, returnObject=T)
    p <- p + coord_flip()
    pdf(paste0('../results_sfam/figs_subClusters/dotplots/', fam.lbl, '.pdf'))
    plot(p)
    dev.off()
}

for (i in 1:nrow(dat)) {
    x <- dat[i, ]
    xl <- number2list(x, group)
    fam.lbl <- rownames(dat)[i]
    p <- boxplot4list(xl, ylab='Expression (scaled)', main=fam.lbl, jitter.width=.2, col.point=rgb(0, 0, .2, .2),
                                            addN=T, returnObject=T)
    p <- p + coord_flip()
    pdf(paste0('../results_sfam/figs_subClusters/boxplots/', fam.lbl, '.pdf'))
    plot(p)
    dev.off()
}
