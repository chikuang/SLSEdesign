#' Calculate the D-optimal design under the SLSE
#'
#' @param FUN The function to calculate the derivative of the given model.
#' @param N The number of sample points in the design space
#' @param u The discretized design space
#' @param tt The level of skewness
#' @param theta The parameter value of the model
#'
#' @details TODO
#'
#' @import CVXR
#'
#' @return CVX object that contains the solved optimal design
#'
#' @export

Dopt <- function(N, u, tt, FUN, theta){
  w <- Variable(N)
  del <- Variable(1)
  n <- length(theta)
  # Set up constraints --------------------------------------------------------------------------
  constraint1 <- lapply(1:N,
                        function(x){
                          -w[x] <= 0
                        })
  constraint2 <- list({
    for (i in 1:N) {
      f <- FUN(u[i], theta)
      g1 <- g1 + w[i] * f
      G2 <- G2 + w[i] * f %*% t(f)
    }

    # ff <- lapply(1:N, function(x){
    #  FUN(u[x], theta)
    # }) %>% list.cbind()
    # g1 <- colSums(w * t(ff))
    # G2 <- (w * ff) %*% t(ff)
    B <- rbind(cbind(1, sqrt(tt) * t(g1)),
               cbind(sqrt(tt) * g1, G2))
    -log_det(B) <= del
  })

  constraint3 <- list( t(w) %*% rep(1, N) == 1)

  # Solve the optimization problem --------------------------------------------------------------
  objective <- Minimize(del)
  problem <- Problem(objective,
                     c(constraint1, constraint2, constraint3))
  solve(problem, num_iter = 10000)

}

