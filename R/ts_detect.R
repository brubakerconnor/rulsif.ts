#' Change point detection
#'
#' Detect change points in a time series using RelULSIF.
#'
#' @param ts Time series to detect change points in. Assumes this is a D by N matrix
#'           where D is the dimension of the time series and N is the number of
#'           time points.
#' @param window_size The length of the sub-sequences generated from the series. Default `5`.
#' @param step How many sub-sequences forward and backward to from a time point
#'             to compute a score from. Default is 10% of the length of the series
#'             if not specified.
#' @param alpha Relative parameter. Default `0.05`. Setting to `0` recovers
#'              ordinary unconstrained least squares importance fitting.
#' @param k Number of basis functions. Default is minimum of `100` and dimension of
#'          the time series.
#' @param n_folds Number of folds to use in determining optimal kernel bandwidth
#'                and lambda parameter in RULSIF.
#'
#' @return A vector of change point scores for each time point starting at t = step
#'         and ending at t = (number of time points) - window_size.
#' @export
#'
#' @examples
#' t <- c(rnorm(150, mean = 0), rnorm(150, mean = 5), rnorm(150, mean = 1))
#' ts_detect(t)
ts_detect <- function(ts, window_size = 5, step = NULL,
                      alpha = 0.05, k = 100, n_folds = 5) {
    # parameters
    dim_ts <- dim(X)[1]
    N_time_points <- dim(X)[2]

    # compatibility checks on input variables
    # window_size
    if (window_size <= 0 || w %% 1 != 0) {
        stop("Parameter window_size must be a positive integer.")
    } else if (window_size > (N_time_points - 1)) {
        stop("Parameter window_size too large. Try reducing it.")
    }

    # step
    if (is.null(step)) {
        step <- 0.1 * N_time_points
    } else if (step <= 0 || step %% 1 != 0) {
        stop("Parameter step must be a positive integer.")
    } else if (step > N_time_points - window_size) {
        stop("Parameter step is too large. Try reducing it.")
    }

    # alpha
    if (alpha <= 0 || alpha > 1 || length(alpha) != 1) {
        stop("Parameter alpha must be in (0, 1].")
    }

    # k
    if (k <= 0 || k %% 1 != 0) {
        stop("Parameter k must be a positive integer.")
    }

    # n_folds
    if (n_folds <= 0 || n_folds %% 1 != 0) {
        stop("Parameter n_folds must be a positive integer.")
    } else if (n_folds > N_time_points) {
        stop("Parameter n_folds exceeds number of points in the input time series.")
    }

    # constructing sliding window
    sw <- sliding_window(X = ts, window_size = window_size)
    n_samples <- dim(sw)[2]

    # compute change-point scores
    scores <- vector()
    t <- step + 1
    while(t + step - 1 <= n_samples) {
        y <- sw[ , (t - step):(step + t - 1), drop = FALSE]
        y_nu <- y[ , 1:step]
        y_de <- y[ , (step + 1):ncol(y)]

        out <- RelULSIF(y_nu, y_de)
        scores <- c(scores, out$rPE)
        t <- t + 1
    }

    # return
    return(scores)
}
