#' Calculate the A-optimal design under the SLSE
#'
#' @param FUN The function to calculate the derivative of the given model.
#' @param N The number of sample points in the design space
#' @param u The discretized design space
#' @param tt The level of skewness
#' @param theta The parameter value of the model
#' @param num_iter Maximum number of iteration
#'
#' @details TODO
#'
#' @import CVXR
#' @importFrom pracma blkdiag
#' @importFrom tibble tibble
#' 
#' @return A list that contains 1. Value of the objective function at solution. 2. Status. 3. Optimal design
#' 
#' @export

Aopt <- function(N, u, tt, FUN, theta, num_iter = 2500){
  
  n <- length(theta)
  g1 <- matrix(0, n, 1)
  G2 <- matrix(0, n, n)
  obj_val <- 0
  C <- rbind(0, diag(1, n, n))
  
  w <- Variable(N)
  del <- Variable(1)
  
  # Set up constraints 
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

    B <- rbind(cbind(1, sqrt(tt) * t(g1)),
               cbind(sqrt(tt) * g1, G2))

    # C <- blkdiag(matrix(0), diag(1, n))

    for(k in 1:n){
      obj_val <- obj_val + matrix_frac(C[, k], B)
    }
    obj_val <= del
  })

  constraint3 <- list( t(w) %*% rep(1, N) == 1)

  # Solve the optimization problem 
  objective <- Minimize(del)
  problem <- Problem(objective,
                     c(constraint1, constraint2, constraint3))
  res <- CVXR::solve(problem, num_iter = num_iter)

  # figure out the location of the design points
  tb <- tibble(location = u,
                weight = c(res$getValue(w)))
  list(val = res$value, status = res$status, design = tb)
}
