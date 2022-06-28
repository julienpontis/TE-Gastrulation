row_center <- function(X) {
    
    ## Description: centers the matrix X by rows
        # Paras: X: matrix with rows to be centered.
    stopifnot(is.matrix(X))

    return (t(scale(x = t(X), center = T, scale = F)))
}
