#' Calculate the loss function of the D-optimal design
#'
#' @param N The number of sample points in the design space
#' @param w_hat The weight of each design points
#' @param FUN The function to calculate the derivative of the given model.
#' @param u The discretized design space
#' @param tt The level of skewness
#' @param theta The parameter value of the model
#' @param A The calculated covariance matrix
#' @param q The penality term
#' @param S The design space
#' @param phi The loss loss function on each design point
#'
#' @details TODO
#'
#' @import CVXR
#' @importFrom pracma blkdiag
#' 
#' @return The loss of the model at each design points
#'
#' @export


plot_direction_Dopt <- function(u, w_hat, tt, FUN, A, phi, theta, q, S, N){
  n <- length(theta)
  mini <- min(phi - q)
  y <- rep(0, N)
  fx <- function(x){
    FUN(x, theta)
  }

  ff <- function(x){
    I <- rbind(cbind(1, sqrt(tt) * t(fx(x))),
               cbind(sqrt(tt) * fx(x), fx(x) %*% t(fx(x))))
    sum(diag(solve(A, I))) - q
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
