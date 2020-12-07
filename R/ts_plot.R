ts_plot <- function(ts, step, scores, change_points) {
    # dimension parameters
    ts_dim <- dim(ts)[1]
    n_time_points <- dim(ts)[2]

    # graphical parameters
    graphics::par(mfrow = c(ts_dim + 1, 1),
        oma = rep(1, 4),
        mar = rep(1.2, 4))

    # plot each time series
    for (s in 1:ts_dim) {
        # plot the time series
        plot(1:n_time_points, ts[s, ], type = 'l',
             xlab = NULL, ylab = paste("y", s, sep = ''),
             xlim = c(1, n_time_points),
             ylim = c(1.5 * min(ts[s, ]), 1.5 * max(ts[s, ])))

        # plot rectangles to highlight change points
        graphics::rect(xleft = change_points,
             ybottom = rep(1.5 * min(ts), length(change_points)),
             xright = change_points + 1,
             ytop = rep(1.5 * max(ts), length(change_points)),
             col = '#DD666666', border = NA)
    }

    # plot scores
    plot(1:length(scores) + step, scores,
         type = 'l',
         xlab = 'time', ylab = 'rPE',
         xlim = c(1, n_time_points),
         col = 'red')
}
