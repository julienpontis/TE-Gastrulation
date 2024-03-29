getSignifActivities <- function(E_centered, N, treatment_group, control_group) {
    # Estimates and returns the activities of TE subfamilies by multiple linear regression
    # Each activity coefficient is tested agains the null hypothesis that it is equal to zero
    # A Benjamini-Hochberg procedure to control the FDR is employed.

    # E centered is the matrix of gene expression values with genes as rows and samples as columns
    # N is the predictor matrix with genes as rows and TE subfamilies as columns
    # treatment_group is vector with the name of treated samples (correspond to columns in E_centered)
    # control_group: same for control samples

    print('Starting new sample')

    # defensive programming
    stopifnot(dim(E_centered)[1]==dim(N)[1])

    stopifnot(all(sapply(treatment_group, function(x) x %in% colnames(E_centered))))
    stopifnot(all(sapply(control_group, function(x) x %in% colnames(E_centered))))

    stopifnot(all(sapply(control_group, function(x) ! x %in% treatment_group)))

    # generating the difference in expression vector

    control = c(sapply(control_group, (function(x) E_centered[, x])))
    treatment = c(sapply(treatment_group, (function(x) E_centered[, x])))
    delta = treatment-control
               
    # generating the augmented N matrix to match the dimensions of delta
    n_replicates = length(control_group)
    N_augmented = N

    i = 1
    while (i < n_replicates) {
        N_augmented = rbind(N_augmented, N)  
        i = i+1
    }

                       
    # assembling the design matrix, regression delta as a function of all columns in N_augmented
    design_matrix = as.data.frame(cbind(delta, N_augmented))
    fit = lm(delta ~ ., data = design_matrix)
    coefs = as.data.frame(summary(fit)$coefficients)

    # FDR adjustment

    # ignoring the intercept (-1)
    p_adj = p.adjust(coefs[-1, 4], method='BH')

    # adding an NA to match the dimensions of `coefs` and appending
    p_adj = append(p_adj, NA, 0)
    coefs$p_adj = p_adj

    res = list()
    res$coefs = coefs
    res$control_group = control_group
    res$treatment_group = treatment_group
    res$r.squared = summary(fit)$r.squared
    res$ajd.r.squared = summary(fit)$adj.r.squared
    res$aliased = summary(fit)$aliased
    res$fstatistic = summary(fit)$fstatistic                                
    res$rss = deviance(fit)
    res$p = ncol(N)
    res$N = nrow(E_centered)
    return(res)
}
