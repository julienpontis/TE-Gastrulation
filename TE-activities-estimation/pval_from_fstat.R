pval_from_fstat <- function(res) {
    ## Description: returns the pvalue corresponding to the f-statistic of the complete versus the null model (mean only)
        ## Para: -res: result returned by getSignifActivities
    
    pval = pf(q = res$fstatistic['value'], df1 = res$fstatistic['numdf'], df2 = res$fstatistic['dendf'], lower.tail = F)
    return(pval)
}