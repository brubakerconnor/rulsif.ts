#' Relative unconstrained least squares importance fitting
#'
#' @param Xnu Samples from numerator probability density
#' @param Xde Samples from denomenator probability density
#' @param Xce Matrix of centers
#' @param sigma Scalar or vector; Gaussian kernel bandwidth(s)
#' @param lambda Scalar or vector; regularization parameter(s)
#' @param alpha Scalar; relative parameter in [0, 1)
#' @param k Integer; number of basis functions
#' @param n_folds Integer; number of folds to use in cross fold validation
#'
#' @return List of 3
#' - opt_sigma: chosen sigma from CV
#' - opt_lambda: chosen lambda from CV
#' - rPE: relative Pearson divergence
#' @export
#'
#' @examples
RelULSIF <- function(Xnu, Xde, Xce = NULL, sigma = NULL, lambda = NULL,
                   alpha = 0.01, k = 100, n_folds = 5) {
    # compatibility checks on data types
    if (is.vector(Xnu)) {
        Xnu <- as.matrix(Xnu)
    }
    if (is.vector(Xde)) {
        Xde <- as.matrix(Xde)
    }

    # dimension compatibility checks
    dim_nu <- dim(Xnu)[1]
    dim_de <- dim(Xde)[1]
    n_nu <- dim(Xnu)[2]
    n_de <- dim(Xde)[2]
    if (dim_nu != dim_de) {
        stop("Number of rows must match between Xnu and Xde.")
    }

    # ensure that k does not exceed dimension of the data
    k <- min(k, dim_nu)

    # initialize centers for Gaussian kernel if NULL
    if (is.null(Xce)) {
        Xce <- Xnu[ , sample(n_nu, size = k, replace = TRUE), drop = FALSE]
    } else{
        # check compatibility
        if (dim(Xce)[1] != dim(Xnu)[1] || dim(Xce)[1] != dim(Xde)[1]) {
            stop("Dimensions of Xce and Xnu or Xde incompatible.")
        }
    }

    # distance matrices
    dist_nu <- comp_dist(Xnu, Xce)
    dist_de <- comp_dist(Xde, Xce)

    # check if search on sigma and lambda is needed
    if (length(sigma) == 1 && length(lambda) == 1) {
        # proceed with user specified sigma and lambda values
        opt_sigma <- sigma
        opt_lambda <- lambda
    } else {
        # first check if user specified sigma/lambda search ranges
        sigma_search_vec <- sigma
        lambda_search_vec <- lambda
        if (is.null(sigma)) {
            med <- comp_median(cbind(Xde, Xnu))
            sigma_search_vec <- med * seq(0.6, 1.4)
        }
        if (is.null(lambda)) {
            lambda_search_vec <- 10 ^ seq(-3, 1, 1)
        }
        cv_out <- sigma_lambda_grid_search(n_nu, n_de, k, dist_nu, dist_de,
            sigma_search_vec, lambda_search_vec, alpha, n_folds)
        opt_sigma <- cv_out$opt_sigma
        opt_lambda <- cv_out$opt_lambda
    }

    # compute results
    ker_Xnu <- gaussian_kernel(t(dist_nu), opt_sigma)
    ker_Xde <- gaussian_kernel(t(dist_de), opt_sigma)

    H <- ((1 - alpha) / n_de) * (ker_Xde %*% t(ker_Xde)) +
        (alpha / n_nu) * (ker_Xnu %*% t(ker_Xnu))
    h <- rowMeans(ker_Xnu)
    theta <- solve(H + diag(k) * opt_lambda) %*% h

    g_nu <- t(theta) %*% ker_Xnu
    g_de <- t(theta) %*% ker_Xde
    rPE <- mean(g_nu) - 0.5 * ((alpha * mean(g_nu ^ 2)) +
        (1 - alpha) * mean(g_de ^ 2)) - 0.5

    # return
    return(list(
        opt_sigma = opt_sigma,
        opt_lambda = opt_lambda,
        rPE = rPE
    ))
}
