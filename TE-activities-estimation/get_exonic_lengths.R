get_exonic_lengths <- function(org='hg19') {
    ## Description: returns the cumulated exonic length for each gene
        ## Para: -org: organism
    require(GenomicFeatures)
    org <- match.arg(org)
    ens <- makeTxDbFromUCSC(genome=org,
                            tablename='ensGene',
                            url = 'http://genome-euro.ucsc.edu/cgi-bin/'
                           )
    exonic <- exonsBy(ens, by='gene')
    red.exonic <- reduce(exonic)
    exon.lengths <- sum(width(red.exonic))
    # Removing ENSEMBL gene versions
    names(exon.lengths) <- unlist(lapply(strsplit(names(exon.lengths), '\\.'), function(x) x[[1]]))
    return(exon.lengths)
}