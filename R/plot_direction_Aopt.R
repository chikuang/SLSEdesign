#' Calculate the loss function of the A-optimal design
#'
#' @param design The resulted design that contains the design points and the associated weights
#' @param FUN The function to calculate the derivative of the given model.
#' @param tt The level of skewness
#' @param theta The parameter value of the model
#' @param A The calculated covariance matrix
#' @param phi The loss loss function on each design point
#'
#' @details This function produces the figure for the directional derivative of the given A-optimal design of the compact supports. According to the general equivalence theorem, for an optimal design, all the directional derivative should be below zero line.
#'
#' @import CVXR
#' @importFrom pracma blkdiag
#'
#' @return The plot of the directional derivative of a A-optimal design
#'
#' @example
#' \dontrun{
#'   TODO
#' }
#'
#' @export

plot_direction_Aopt <- function(design, tt, FUN, A, phi, theta){

  u <- design$location
  w_hat <- design$weight
  N <- length(u)
  S <- u[c(1, length(u))]

  n <- length(theta)

  y <- rep(0, N)
  fx <- function(x){
    FUN(x, theta)
  }

  AI <- solve(A)
  C <- blkdiag(matrix(0), diag(1, n))
  q <- sum(diag(C * AI * t(C)))
  mini <- min(phi - q)
  ff <- function(x){
    I <- rbind(cbind(1, sqrt(tt) * t(fx(x))),
               cbind(sqrt(tt) * fx(x), fx(x) %*% t(fx(x))))
    sum(diag(I %*% AI %*% t(C) %*% C %*% AI)) - q
  }

  for(i in 1:N){
    y[i] <- ff(u[i])
  }

  plot(u, y, type = "l",  col = "black", xlim  = S, ylim  = c(mini*1.1, 1),
       xlab = "Design space", ylab = expression(d(x, theta)))
  points(u, phi - q,  col = "blue")
  points(u[which(w_hat > 1E-4)], (phi - q)[which(w_hat > 1E-4)], col = "red", cex = 2)
  abline(u, 0)
  legend("bottomright",
         legend=c("Reference Line", expression(d(x, theta)), "discretized point", "Support point"),
         col = c("black", "blue", "blue", "red"),
         lty = c(1, 1, NA, NA), pch = c(NA, NA, 1, 1), cex = 0.8)
}
