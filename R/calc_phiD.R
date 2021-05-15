#' Calculate the loss function of the D-optimal design
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
#' 
#' @return The loss of the model at each design points
#'
#' @export

calc_phiD <- function(design, FUN, tt, A, theta){
  u <- design$location
  w_hat <- design$weight
  N <- length(u)
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
