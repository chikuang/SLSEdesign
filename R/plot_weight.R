#' Plot the weight distribution of the optimal design
#'
#' @param u The discretized design space
#' @param w_hat The estimated weight of each design point
#' @param S The design space
#'
#' @details TODO
#'
#' @importFrom graphics abline legend lines points
#'
#' @return The plot that shows the given optimal design
#'
#' @export


plot_weight <- function(u, w_hat, S){
  plot(u, w_hat, type = "p", xlim = S, ylim = c(0, 1),
       xlab  = "Design space", ylab = "weight")
  lines(u, w_hat, type = "l")
  points(u[which(w_hat > 1E-4)], w_hat[which(w_hat > 1E-4)],
         col = "red", cex = 2)
}
