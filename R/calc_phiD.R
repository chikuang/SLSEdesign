#' Calculate the loss function of the D-optimal design
#'
#' @param N The number of sample points in the design space
#' @param w_hat The weight of each design points
#' @param FUN The function to calculate the derivative of the given model.
#' @param u The discretized design space
#' @param tt The level of skewness
#' @param theta The parameter value of the model
#' @param A The calculated covariance matrix
#'
#' @details TODO
#'
#' @import CVXR
#' 
#' @return The loss of the model at each design points
#'
#' @export

calc_phiD <- function(N, u, w_hat, FUN, tt, A, theta){
  w_hat <- c(w_hat)
  n <- length(theta)
  g1 <- matrix(0, n, 1)
  G2 <- matrix(0, n, n)
  phi_D <- rep(0, N)
  for( i in 1:N){
    f <- FUN(u[i], theta)
    I <- rbind(cbind(1, sqrt(tt) * t(f)),
               cbind(sqrt(tt) * f, f %*% t(f)))
    phi_D[i] = sum(diag(solve(A, I)));
  }
  phi_D
}
