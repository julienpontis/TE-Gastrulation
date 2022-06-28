preprocess_E_N_for_activities <- function(count_table, N, 
                                          count_threshold=c(10, 1), 
                                          N_threshold = 150,
                                          pseudocount_quantile=0.05,
                                          log_tpm_plot_path,
                                          import_path,
                                          preds_protected=NULL,
                                          redownload_exons=FALSE
                                          ) {
    ## Descrption: filter, normalize and process a table of raw counts with genes as rows and samples as columns to 
    ## the expression matrix E. Build a corresponding N matrix with remaining genes as rows and TEs as columns.

        ## Para: -count_table: gene count matrix (rows are genes, columns are samples).
                #-N: per gene occurences of TE integrants. Genes are rows, TE subfamilies are columns
                #-count_threhold: c(a,b) with a the minimum number of counts in at least b genes to keep a gene present in E
                #-N_threshold:  minimum sum over all genes that a TE subfamily has to reach in order to be kept as a predictor in matrix N
                #-pseudocount_quantile: quantile of non-zero values used to compute per-sample pseudocounts prior to log-transformation
                #-log_tpm_plot_path: path where histograms of expression values (in logged TPM) are drawn (useful to check quasi-normality of expression values)
                #-import_path: path pointing to the other R functions required for the pre-processing
                #-preds_protected: vector containing the names of predictors that cannot be filtered out based on N_threshold
                #-redownload_exons: TRUE: exons are redownloaded to compute TPMs. FALSE: exons are not redownloaded to compute TPMs.
        

    stopifnot(typeof(count_table)=='integer' & is.matrix(count_table) & all(count_table>=0))
    stopifnot(is.matrix(N) & all(N)>=0)
    
    # importing useful functions

    source(file.path(import_path, 'get_pseudocounts.R'))
    source(file.path(import_path, 'get_exonic_lengths.R'))
    source(file.path(import_path, 'get_tpm_local.R'))
    source(file.path(import_path, 'row_center.R'))

    # filtering genes not expressed at enough counts in at least ... samples
    genes_kept = colSums(apply(count_table,
                               MARGIN=1,
                               (function(x) x>=count_threshold[1])))>=count_threshold[2]
    
    count_table = count_table[genes_kept, ]


    # finding genes in count_table that are in N
    in_N = rownames(count_table) %in% rownames(N)
    count_table = count_table[in_N, ]

    # ordering rows of N as rows of count_table
    N = N[rownames(count_table), ]

    # Selecting cols in N with at least ... hits

    # saving lost fams to return them
    lost_cols = colSums(N)[which(colSums(N)<N_threshold & ! colnames(N) %in% preds_protected)]
    
    # removing cols of N with 
    N = N[, which(colSums(N)>=N_threshold | colnames(N) %in% preds_protected)]

    # removing genes that now have 0 hits in N:
    to_keep = rowSums(N)>0
    N = N[to_keep, ]

    # in case we lost some genes in N, removing them in the count_table
    count_table = count_table[rownames(N), ]

    # adding pseudocounts to the count_table
    count_table = get_pseudocounts(count_table, pseudocount_quantile)

    # computing log2(tpm), keeping them for the graph
    log_tpm = log2(get_tpm_local(count_table, redownload_exons=redownload_exons))
    gc()

    # plotting log tpm for quality control
    pdf(log_tpm_plot_path)
    hist(log_tpm)
    dev.off()

    # col centering E
    E = scale(log_tpm, center = T, scale = F)
    gc()

    # row centering E
    E = row_center(E)

    stopifnot(all(rownames(E)==rownames(N)))

    return(list('E_centered' = E,
                'N' = N,
                'lost_cols' = lost_cols))
}
