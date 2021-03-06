#' Calculate the loss function of the A-optimal design
#'
#' @param design The resulted design that contains the design points and the associated weights
#' @param FUN The function to calculate the derivative of the given model.
#' @param tt The level of skewness
#' @param theta The parameter value of the model
#' @param A The calculated covariance matrix
#'
#' @details This function calculates the loss function of the design problem under the A-optimality. The loss function under A-optimality is defined as the trace of the inverse of the Fisher information matrix
#'
#' @import CVXR
#' @importFrom pracma blkdiag
#'
#' @return The loss of the model at each design points
#'
#' @examples
#' \dontrun{
#' my_design <- tibble(location = c(0, 180), weight = c(1/2, 1/2))
#' theta <- c(0.05, 0.5)
#' peleg <- function(xi, theta){
#'    deno <- (theta[1] + xi * theta[2])^2
#'    rbind(-xi/deno, -xi^2/deno)
#' }
#' A <- matrix(c(1, 0, 0, 0, 0.2116, 1.3116, 0, 1.3116, 15.462521), byrow = T, ncol = 3)
#' res <- calc_phiA(my_design, theta, peleg, 0, A)
#' res
#' }
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
