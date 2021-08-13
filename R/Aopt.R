#' Calculate the A-optimal design under the second-order Least squares estimator
#'
#' @param FUN The function to calculate the derivative of the given model.
#' @param N The number of sample points in the design space
#' @param u The discretized design space
#' @param tt The level of skewness
#' @param theta The parameter value of the model
#' @param num_iter Maximum number of iteration
#'
#' @details This function calculates the loss function of the design problem under the A-optimality. The loss function under A-optimality is defined as the trace of the inverse of the Fisher information matrix
#'
#' @import CVXR
#' @importFrom pracma blkdiag
#' @importFrom tibble tibble
#'
#' @return A list that contains 1. Value of the objective function at solution. 2. Status. 3. Optimal design
#'
#' @examples
#' #' @example
#' \dontrun{
#'
#' peleg <- function(xi, theta){
#'    deno <- (theta[1] + xi * theta[2])^2
#'    rbind(-xi/deno, -xi^2/deno)
#' }
#' my_design <- Aopt(N = 31, u = seq(0, 180, length.out = 31), tt = 0, FUN = peleg,
#'     theta = c(0.05, 0.5), num_iter = 500)
#' my_design$design
#' my_design$val
#' }
#'
#' @export

Aopt <- function(N, u, tt, FUN, theta, num_iter = 1000){

  n <- length(theta)
  g1 <- matrix(0, n, 1)
  G2 <- matrix(0, n, n)
  obj_val <- 0
  C <- rbind(0, diag(1, n, n))

  w <- CVXR::Variable(N)
  del <- CVXR::Variable(1)

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

    C <- pracma::blkdiag(matrix(0), diag(1, n))

    for(k in 1:n){
      obj_val <- obj_val + CVXR::matrix_frac(C[, k], B)
    }
    obj_val <= del
  })

  constraint3 <- list( t(w) %*% rep(1, N) == 1)

  # Solve the optimization problem
  objective <- CVXR::Minimize(del)
  problem <- CVXR::Problem(objective,
                     c(constraint1, constraint2, constraint3))
  res <- CVXR::solve(problem, num_iter = num_iter)

  # figure out the location of the design points
  tb <- tibble::tibble(location = u,
                weight = c(res$getValue(w)))
  list(val = res$value, status = res$status, design = tb)
}
