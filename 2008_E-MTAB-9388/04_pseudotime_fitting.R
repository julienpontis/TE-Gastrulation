#docker version:
#singularity run docker://registry.c4science.ch/lvg/r-4.1.1:v3

#DESCRIPTION
#we got the code from the authors of the paper

#load libs
#install.packages(c('gam', 'foreach'), lib='R_libs')
library(foreach, lib='R_libs')
library(gam, lib='R_libs')
library(colorout)
library(Seurat)
#library(scran, quietly = TRUE, warn.conflicts = FALSE, lib='R_libs')

#get the data
{
    x <- read.csv('../data/fromAuthors/epi_nm_pseudo.csv')
    dpt <- x[, 2]; names(dpt) <- x[, 3]

    dummy.g <- load('../data/RData/seurat.genes.RData')
    cnts <- GetAssayData(xg, slot = "counts", assay = "SCT")
    cnts <- as.matrix(cnts)
    colnames(cnts) <- gsub('-', '.', colnames(cnts))

    dummy.sf <- load('../data/RData/seurat.sfam.RData')
    cnts.sf <- GetAssayData(xsf, slot = "counts", assay = "SCT")
    cnts.sf <- as.matrix(cnts.sf)
    colnames(cnts.sf) <- gsub('-', '.', colnames(cnts.sf))

    #check names of pseudotimes are in expression
    stopifnot(all(names(dpt) %in% colnames(cnts)))
    stopifnot(all(names(dpt) %in% colnames(cnts.sf)))

    #check selected cell types
    table(xg[, gsub('\\.', '-', names(dpt))][[]]$authCellType)

    #merge data
    cnts <- cnts[, names(dpt)]
    cnts.sf <- cnts.sf[, names(dpt)]
    stopifnot(all(colnames(cnts) == names(dpt))) #check they are in order
    stopifnot(all(colnames(cnts.sf) == names(dpt))) #check they are in order
}

#ensembl and sfam ids
ens <- c('ENSG00000164458', 'ENSG00000019549', 'ENSG00000204531', 'ENSG00000111704', 'ENSG00000136826',
         'ENSG00000163508', 'ENSG00000087510', 'ENSG00000181449')
names(ens) <- c('TBXT', 'SNAI2', 'OCT4', 'NANOG', 'KLF4', 'EOMES', 'TFAP2C', 'SOX2')
sfa <- c('HERVS71-int', 'HERVK-int', 'HERVK11-int', 'HERVH-int', 'L1HS', 'L1PA2', 'L1PA3', 'SVA-D')
names(sfa) <- sfa

#run gam and plot for genes 
Y <- cnts 
t <- dpt
pdf('../results/pseudotime/pseudotimeFromAuthors_genes.pdf')
for (i in 1:length(ens)) {
    vec <- c(ens[i])
    idx<-match(vec, rownames(Y))
    d <- data.frame(z=Y[idx,], t=t)
    fit <- gam(z ~ lo(t), data=d)
    #print(summary(fit))
    #print(fitted(fit))
    #
    xlab <- 'Pseudotime'
    ylab <- names(ens)[i]
    plot(fit, se=T, residuals=F, xlab=xlab, ylab=ylab)
    plot(fit, se=T, residuals=T, xlab=xlab, ylab=ylab)
}
dev.off()

#run gam and plot for TE sfam
Y <- cnts.sf 
t <- dpt
pdf('../results/pseudotime/pseudotimeFromAuthors_TESfam.pdf')
for (i in 1:length(sfa)) {
    vec <- c(sfa[i])
    idx<-match(vec, rownames(Y))
    d <- data.frame(z=Y[idx,], t=t)
    fit <- gam(z ~ lo(t), data=d)
    #print(summary(fit))
    #print(fitted(fit))
    #
    xlab <- 'Pseudotime'
    ylab <- names(sfa)[i]
    plot(fit, se=T, residuals=F, xlab=xlab, ylab=ylab)
    plot(fit, se=T, residuals=T, xlab=xlab, ylab=ylab)
}
dev.off()
