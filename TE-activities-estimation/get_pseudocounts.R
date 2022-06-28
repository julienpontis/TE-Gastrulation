get_pseudocounts <- function(count_table, q=0.05) {
        ## Description: compute pseudocounts on a gene count matrix (genes are rows, samples are columns)
        # to replace zero values before logging

            ## Para: -count_table: gene count table with rows as genes, columns as samples
                    #-q: quantile of non-zero values used as a pseudocount in each sample (default: quantile = 0.05)
        stopifnot(typeof(count_table)=='integer' & is.matrix(count_table))

        pseudocounts = apply(X = count_table, MARGIN = 2, FUN = function(x) quantile(x[x!=0], q))
        
        # Broadcasting only works correctly by rows in R, so when broadcasting by columns
        # one has to transpose to broadcast by rows (instead of cols), do the operation
        # and then tranpose again to get back to the original conformation. 
        return(t(t(count_table) + pseudocounts))
}
