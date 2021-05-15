#' Calculate the loss function of the A-optimal design
#'
#' @param design The resulted design that contains the design points and the associated weights
#' @param FUN The function to calculate the derivative of the given model.
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

calc_phiA <- function(design, theta, FUN, tt, A){
  u <- design$location
  w_hat <- design$weight
  N <- length(u)
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
