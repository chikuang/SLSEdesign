#' Calculate the loss function of the A-optimal design
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
#' @importFrom pracma blkdiag
#'
#' @return The loss of the model at each design points
#'
#' @export

calc_phiA <- function(N, u, w_hat, theta, FUN, tt, A){
  w_hat <- c(w_hat)
  n <- length(theta)
  g1 <- matrix(0, n, 1)
  G2 <- matrix(0, n, n)
  phi_A <- rep(0, N)
  BI <- solve(A)
  C <- blkdiag(matrix(0), diag(1, n))

  for( i in 1:N){
    f <- FUN(u[i], theta)
    I <- rbind(cbind(1, sqrt(tt) * t(f)),
               cbind(sqrt(tt) * f, f %*% t(f)))
    phi_A[i] <- sum(diag(I %*% BI %*% t(C) %*% C %*% BI))
  }
  phi_A
}
