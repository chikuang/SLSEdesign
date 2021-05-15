#' Calculate the loss function of the D-optimal design
#'
#' @param design The resulted design that contains the design points and the associated weights
#' @param FUN The function to calculate the derivative of the given model.
#' @param tt The level of skewness
#' @param theta The parameter value of the model
#' @param A The calculated covariance matrix
#' @param phi The loss loss function on each design point
#'
#' @details TODO
#'
#' @import CVXR
#' @importFrom pracma blkdiag
#' 
#' @return The plot of the directional derivative of a D-optimal design
#'
#' @export


plot_direction_Dopt <- function(design, tt, FUN, A, phi, theta){
  
  u <- design$location
  w_hat <- design$weight
  N <- length(u)
  S <- u[c(1, length(u))]
  
  n <- length(theta)
  q <- n + 1
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