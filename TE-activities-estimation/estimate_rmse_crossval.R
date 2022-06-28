estimate_rmse_crossval <- function(preprocessed, control_group, treatment_group,
                                  N_list, N_metadata,
                                  folds, f) {
        ## Description: computes the RMSE in the context of a parallelizable cross-validation scheme
    
            #-Paras: -preprocessed: data produced by preprocess_E_N_for_activities()
                    #-control_group: vector containing the names of the control samples (column names in preprocessed$E_centered)
                    #-treatment_group: vector containing the sames of the treatment samples (column names in preprocessed$E_centered)
                    #-N_list: list of the different N matrices. Rows must already correspond to those of preprocessed$E_centered
                    #-N_metadata: dataframe with a number of rows equal to the number of elements in N_list
                    #-folds: vector of the same length as the number of observations in E_centered, with integers from 1:nFolds, each defining a fold
                    #-f: fold index for which to estimate the RMSE

    idx_training = folds==f
    idx_test = !idx_training
    # initializing for temporary results
    rmse = N_metadata
    rmse$fold = as.integer(0)
    rmse$rmse = 0       
    rmse['fold'] = f
    for(i in 1:nrow(N_metadata)) {
        rmse[[i, 'rmse']] = estimateRMSE(preprocessed$E_centered, N=N_list[[i]], treatment_group = treatment_group, control_group=control_group, idx_training = idx_training, idx_val = idx_test)
    }
    return(rmse)
}