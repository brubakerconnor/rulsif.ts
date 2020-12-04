gaussian_kernel <- function(distance_matrix, sigma) {
    exp(-distance_matrix / (2 * sigma ^ 2))
}
