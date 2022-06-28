corr_heatmap <- function(X) {
    ## Description: compute the pairwise correlation of matrix X columns 
    ## and plot it as a heatmap.
    ## PARA: matrix X
    library(gplots)
    corr.pair <- cor(X)
    palette <- colorRampPalette(colors = c("blue", "white", "red"))
    ans <- heatmap.2(corr.pair, Rowv = NULL, Colv = NULL, dendrogram = "none", 
              scale = "none", trace = "none",
              col = palette, key.title = "Pearson cor. coef.")
    ans
}


get_tpm <- function(x, org=c('hg19', 'mm9', 'rheMac8')) {
    ## Description: Compute TPM from a table of counts.
    ## Para: -x: table with raw counts, ensembl id as rownames
    ##       -org: organism.
    countTable = x
    require(GenomicFeatures)
    org <- match.arg(org)
    ens <- makeTxDbFromUCSC(genome=org, tablename='ensGene',
                            url = 'http://genome-euro.ucsc.edu/cgi-bin/',
                            goldenPath_url = 'http://hgdownload.cse.ucsc.edu/goldenPath')
    exonic <- exonsBy(ens, by='gene')
    red.exonic <- reduce(exonic)
    exon.lengths <- sum(width(red.exonic))
    names(exon.lengths) <- unlist(lapply(strsplit(names(exon.lengths), '\\.'), function(x) x[[1]]))
    exon.lengths.o <- exon.lengths[rownames(countTable)]
    length_norm <- countTable / exon.lengths.o
    ans <- t(t(length_norm)/colSums(length_norm)) * 1e6 # see https://rpubs.com/bbolker/sweep_divide
    ans
}

# doesnt work yet
get_tpm_mappable <- function(x, org="hg19", readLength=50, 
                             pathToSamba="~/samba") {
    ## Description: normalize by mappable length (TPM) from a table of counts.
    ## Para: -x: table with raw counts, ensembl id as rownames
    ##       -org: organism.


    # Future improvements
    stopifnot(org=="hg19", "get_tpm_mappable not implemented for genomes other than hg19")
    stopifnot(readLength == 50, "get_tpm_mappable not implemented for other read lengths")


    print("loading and computing mappable lengths")

    TE_idx_table <- read.csv(pathToSamba + "/projects/didier/" +
                     "1510_mapAbility_UCSCtrack/data/" + 
                     "hg19_REPEATS_repmask_LTRm_s_20140131_UCSC_unique50.bed",
                 as.is=T, sep="\t", header=F)

    colnames(TE_idx_table) <- c('chr', 'start', 'end_cropped', 'idx')
    TE_idx_table$TE_idx <- paste(TE_idx_table$chr, TE_idx_table$start,
                                 TE_idx_table$end_cropped + 50, sep="-")

    print("subsetting on TEs from the count table")
    TE_idx_table = TE_idx_table[TE_idx_table$TE_idx == rownames(x), ]
    

    print("loading bigwig")
    mappable_lengths <- read.csv(pathToSamba + "/projects/didier/" +
                     "1510_mapAbility_UCSCtrack/data/" + 
                     "avgOverBed50.bed",
                 as.is=T, sep="\t", header=F)
    mappable_lengths <- mappable_lengths[, -c(3, 5, 7)]
    colnames(mappable_lengths) <- 


    stopifnot(all(table(tmp[, 2] == (tmp[, 9] - tmp[, 8]))))


    length_norm <- countTable / exon.lengths.o
    ans <- t(t(length_norm)/colSums(length_norm)) * 1e6 # see https://rpubs.com/bbolker/sweep_divide
    ans
}


get_high_leverage_lm <- function(X, thresh){
    ## Description: identify indices of high leverage observations in linear regression
    ## Para: -X: predictor matrix (predictors = columns)
    ##      -thresh: high leverage threshold (mult. of expected lev. stat.)

    # expected leverage (p+1)/n, ref ISLR)
    p = ncol(X)
    n = nrow(X)

    exp_lev = (p + 1) / n
    thresh_lev = thresh * exp_lev

    # identifying high leverage observations
    levs = hat(X)
    high_levs = which(levs>thresh_lev)
    high_levs
}

get_aliased_preds_lm <- function(model) {
    ## Description: identify aliased predictors (causing singularities) in a linear
    ## regression model trained by lm()
    ## Para: -model: linear regression model trained with lm()

    aliases = alias(model)
    aliased_preds = row.names(aliases$Complete)
    aliased_preds
}

# this function is pure bullshit and is only there to keep a souvenir of my
# previous misunderstanding of 1) high leverage points and 2) collinear predictors
# in least square regression with or without regularization. In particular, 
# high leverage points do NOT require calculating the model (of course, they only
# depend on the predictor matrix). Thus I do NOT need to care whether there are
# singularities in the model to get my leverage points. I fooled myself because of 
# blind trust into the hatvalues() function (presented in ISLR), when I can instead
# use hat().
purge_aliases_and_high_leverage_lm <- function(f, reg_mat, thresh){
    ## Description: identify aliased predictors (causing singularities), remove them,
    ## identify high leverage observations, remove them, identify new aliased predictors,
    ## remove them, identify high leverage observations, remove them until no 
    ## aliases or high leverage observations are left. Returns a list for
    ## high leverage points and another one for aliases.
    ## Para: -f: formula to train lm()
    ##       -reg_mat: dataframe containing observations and predictors as columns
    ##       -thresh: high leverage threshold (mult. of expected lev. stat.)
    
    #initialization
    i = 1
    aliased_preds = c()
    high_levs = c()
    repeat{
        # identification and removal of aliased predictors
        cat("Training model nr", i, "\n")
        model = lm(f, data = reg_mat, singular.ok = T)
        new_aliased_preds = get_aliased_preds_lm(model) 
        aliased_preds = append(aliased_preds, new_aliased_preds)
        cat("New aliased preds: ", new_aliased_preds, "\n")
        #updating reg_mat
        reg_mat = reg_mat[, which(!names(reg_mat) %in% aliased_preds)]

        #identification and removal of high leverage observations
        model = lm(f, data = reg_mat, singular.ok = F)
        new_high_levs = get_high_leverage_lm(model, thresh)
        cat("New high leverage observations: ", names(new_high_levs), "\n")

        reg_mat = reg_mat[-new_high_levs, ]
        high_levs = append(high_levs, names(new_high_levs))

        if(length(new_high_levs) == 0) {
            break
        }
        i = i+1
    }
    ans = list(aliased_preds, high_levs)
    ans 
}

getActivities <- function(expr, X, alpha = 0.1, k = 5,
                           ncores = 1, log_id = '') {
    ## Description: estimate transcription factor activities by regressing each
    ## column of the expression matrix (1 col per sample) onto per promoter
    ## transcription factor binding sites. Multiple linear regression with the
    ## elastic net regularization is performed by the glmnet package.


    ## Para: -expr: row-centered TPM values in a gene (rows) x sample (columns) matrix.
    ##       -X: per promoter (rows) TFBS (columns) counts. Column-centered.
    ##       -k: number of folds for cross-validation
    ##       -alpha: trades off between ridge (alpha = 0) and lasso (alpha = 1) penalties
    ##       -ncores: number of cores for multiprocessing (1 -> no multiprocessing)


    ## Author: Cyril Pulver (cyril.pulver@epfl.ch)



    require(glmnet)
    require(Biobase)
    require(doParallel)
    print("Going parallel !")
    cl <- makeCluster(ncores, type="FORK", outfile=paste(c("./getActivities_log", log_id, ".txt"), collapse=''))

    n = ncol(expr)
    cat("Regressing activities of TFs for ", n, " samples\n")
    # regularized least squares on single expression vector
    f <- function(y) {
        print("started new sample")
        cv.out = cv.glmnet(x = X, y = y, alpha = alpha, 
                           nfolds = k)
        coefs = predict(cv.out$glmnet.fit, type = "coefficients", s = cv.out$lambda.min)
        r_squared = cv.out$glmnet.fit$dev.ratio[which(cv.out$glmnet.fit$lambda == cv.out$lambda.min)]
        lambda_min = cv.out$lambda.min
        cvm_lambda_min = cv.out$cvm[which(cv.out$lambda == cv.out$lambda.min)]
        cvm_sd = cv.out$cvsd[which(cv.out$lambda == cv.out$lambda.min)]
        list("coefs" = coefs, 
             "r_squared" = r_squared, 
             "lambda_min" = lambda_min,
             "cvm_lambda_min" = cvm_lambda_min,
             "cvm_sd" = cvm_sd)
    }
    res = list()
    res <- parApply(cl, expr, MARGIN = 2, FUN = f)

    # shutting down parallel clusters
    stopCluster(cl)

    # renaming activities as they were entered at first (removing the "X"
    # added by default by glmnet
    #names(res) = substr(names(res), 2, max(sapply(names(res), nchar)))

    # reshaping results into a five-list object: ans$bestlambda, ans$r2,
    # ans$test_err_est, ans$test_err_se and ans$activities
    ans <- list()
    ans$r2 <- as.vector(subListExtract(res,
                                  "r_squared",
                                  simplify = T))
    ans$lambda_min <- as.vector(subListExtract(res,
                                  "lambda_min",
                                  simplify = T))
    ans$cvm_lambda_min <- as.vector(subListExtract(res,
                                  "cvm_lambda_min",
                                  simplify = T))
    ans$cvm_se_lambda_min <- as.vector(subListExtract(res,
                                  "cvm_sd",
                                  simplify = T))


    # Building the TF activities (rows) by samples (columns) table
    ans$activities <- t(do.call(rbind, lapply(subListExtract(res, "coefs"), as.vector)))
    rownames(ans$activities) <- rownames(subListExtract(res, "coefs")[[1]])
    
    return(ans)
}

plotActivities <- function(act, TF, metadata, sort_by, decreasing=F, 
                           palette_type='sequential', 
                           palette_add_black = T, palette_black_idx=1,
                           title=NULL, act_lim=NULL) {


    ## Description: scatter plot of TF activities per sample 
    
    ## Para: -act: TF x samples matrix, with row and col names set
    ##       -TF: string for TF of interest, e.g "MYC"
    ##       -metadata: samples x metadata dataframe. Metadata should be set
    ## as factors
    ##       -sort_by: string corresponding to a column of metadata
    ##       -decreasing: bool to sort samples by metadata[sort_by] in decr. order.
    ##       -palette_type: string, one of c("sequential", "diverging", "categorical")
    ##       -palette_add_black: bool to enable setting a color to black (useful
    ## for baselines)
    ##       -palette_black_idx: category to plot in black
    ##       -title: title if it is to be different than the TF name
    ##       -act_lim: limits of activities to plot, useful to get comparable y axes

    ## Author: Cyril Pulver (cyril.pulver@epfl.ch)



    require(RColorBrewer)

    # to reset after fun
    old_pal = palette()
    old_par = par(no.readonly=T)

    # matching TF name to idx
    tf_idx = which(rownames(act)==TF)

    # sorting RNA seq samples by selected metadata
    sample_order_idx = seq(1, dim(act)[2])
    if(!is.null(sort_by)){
        sample_order_idx = order(metadata[sort_by], decreasing=decreasing)
    }

    # generating a palette
    new_pal = palette("default")
    # sorting by categorical factor -> n colors in the palette = n factors
    if(is.factor(metadata[[sort_by]])) {
        if (palette_type=='sequential'){
                new_pal = brewer.pal(n = length(levels(metadata[[sort_by]])), name = "OrRd")
            } else if (palette_type=='diverging') {
            new_pal = brewer.pal(n = length(levels(metadata[[sort_by]])), name = "Spectral")
            } else if (palette_type=='qualitative') {
            new_pal = brewer.pal(n = length(levels(metadata[[sort_by]])), name = "Set1")
        }
    }

    # In case of ref level chosen as black
    if (palette_add_black){
        new_pal[palette_black_idx] = "black"
    }
    palette(new_pal)

    
    # enlarging margins for legend out of plot, forcing square plot area
    par(mai=c(0.2, 0.8, 0.4, 1.2), pty='s')
        

    # title
    main_title=title
    if (is.null(title)) {
        main_title=TF
    }

    plot(x = seq(from = 1, to = ncol(act)),
         y = act[tf_idx, sample_order_idx],
         col = metadata[[sort_by]][sample_order_idx],
         main = main_title,
         xlab = paste("Samples (ordered by ", sort_by, ")", sep=""),
         ylab = "Activity (deviation from average sample activity)",
         ylim=act_lim,
         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
    legend("right", 
           legend = levels(metadata[[sort_by]]), 
          col = new_pal, pch = 'o', inset=c(-0.26, 0), 
            box.lwd = 0, bty = "n", xpd=T, cex=1.2)
    abline(a = 0, b = 0, col = "blue", lty=2, lwd = 2)

    # resetting palette and par
    #palette(old_pal)
    #par(old_par)
}
