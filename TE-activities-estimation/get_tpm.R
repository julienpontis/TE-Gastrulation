get_tpm <- function(x, org=c('hg19', 'mm9', 'rheMac8')) {
    
    ## Description: Compute TPM from a table of counts.
        ## Para: -x: table with raw counts, ensembl id as rownames
        ##       -org: organism.

    countTable = x
    require(GenomicFeatures)
    org <- match.arg(org)
    ens <- makeTxDbFromUCSC(genome=org,
                            tablename='ensGene',
                            url = 'http://genome-euro.ucsc.edu/cgi-bin/',
                            #goldenPath_url = 'http://hgdownload.cse.ucsc.edu/goldenPath',
                            )
    exonic <- exonsBy(ens, by='gene')
    red.exonic <- reduce(exonic)
    exon.lengths <- sum(width(red.exonic))
    names(exon.lengths) <- unlist(lapply(strsplit(names(exon.lengths), '\\.'), function(x) x[[1]]))
    exon.lengths.o <- exon.lengths[rownames(countTable)]
    length_norm <- countTable / exon.lengths.o
	        ans <- t(t(length_norm)/colSums(length_norm)) * 1e6 # see https://rpubs.com/bbolker/sweep_divide
	        ans
}
