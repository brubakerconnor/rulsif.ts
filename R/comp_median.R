comp_med <- function(x) {
    # dimension parameters
    d <- dim(x)[1]
    n <- dim(x)[2]

    # calculate median
    G <- matrix(rep(colSums(x ^ 2), times = n), nrow = n, byrow = T)
    D <- G - (2 * t(x) %*% x) + t(G)
    D[upper.tri(D)] <- 0

    # return
    return(sqrt(0.5 * median(D[D > 0])))
}
