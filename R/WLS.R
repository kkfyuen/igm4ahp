#' Numerical Weighted Least Squares (WLS) Optimization Solver for AHP
#'
#' Solves the non-linear Weighted Least Squares (WLS) optimization problem for a
#' Pairwise Reciprocal Matrix (PRM) in the Analytic Hierarchy Process (AHP), originally
#' formulated by Chu et al. (1979). Uses stochastic or derivative-free algorithms
#' via the \code{nloptr} package (default: \code{"NLOPT_GN_ISRES"}) to find the global
#' minimizer of the Weighted Least Squares Error (WLSE), serving as a numerical benchmark
#' against exact closed-form Inverse Gram Matrix (IGM) methods.
#'
#' @param A A numeric square matrix of dimension \code{n x n} representing a Pairwise Reciprocal Matrix (PRM).
#' @param maxeval Integer. Maximum number of function evaluations allowed for the optimization solver (default: \code{500000}).
#' @param xtol_rel Numeric. Relative convergence stopping tolerance for the objective function (default: \code{1.0e-6}).
#' @param algorithm Character string. The algorithm code recognized by \code{nloptr} (default: \code{"NLOPT_GN_ISRES"}).
#'
#' @return An object of class \code{nloptr} containing the optimization results, including the
#'   optimal priority weight vector (\code{$solution}) and the objective function value (\code{$objective}).
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Chu, A. T. W., Kalaba, R. E., & Spingarn, K. (1979). A comparison of two methods for determining the weights of belonging to fuzzy sets. \emph{Journal of Optimization Theory and Applications}, 27(4), 531-538.
#'
#' Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. \emph{arXiv preprint arXiv:2401.01190}. \url{https://doi.org/10.48550/arXiv.2401.01190}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Define a sample 3x3 Pairwise Reciprocal Matrix
#'   A <- matrix(c(1, 2, 5,
#'                 1/2, 1, 3,
#'                 1/5, 1/3, 1), nrow = 3, byrow = TRUE)
#'
#'   # Run WLS numerical optimization
#'   res <- WLSopt(A)
#'   print(res$solution)
#' }
WLSopt <- function(A, maxeval = 500000, xtol_rel = 1.0e-6, algorithm = "NLOPT_GN_ISRES") {

  # Extract matrix order (dimension n)
  n <- nrow(A)

  # Internal wrapper function evaluating the Weighted Least Squares Error (WLSE)
  WLSE.Obj <- function(W) {
    WLSE(W, A)
  }

  # Equality constraint function enforcing unit simplex normalization: sum(w_i) - 1 = 0
  eqc <- function(W) {
    sum(W) - 1
  }

  # Assemble control parameter list for nloptr solver
  opts <- list(
    "algorithm" = algorithm,
    "xtol_rel"  = xtol_rel,
    "maxeval"   = maxeval
  )

  # Execute derivative-free non-linear numerical optimization via nloptr.
  # Initial starting vector (x0) is initialized using rounded PNIGM (Pseudo-Normalized Inverse Gram Matrix)
  # closed-form weights for stable global convergence on higher-dimensional simplexes (Yuen, 2024).
  result <- nloptr::nloptr(
    x0         = round(PNIGM(A, detail = FALSE), 8), # Initial search point from PNIGM (8-decimal precision)
    eval_f     = WLSE.Obj,                                           # Objective function to minimize: f(W) = sum_{i,j} (w_i - a_ij * w_j)^2
    lb         = rep(0, n),                                           # Lower bounds: w_i >= 0
    ub         = rep(1, n),                                           # Upper bounds: w_i <= 1
    eval_g_eq  = eqc,                                                 # Equality constraint function: sum(w_i) = 1
    opts       = opts                                                 # Solver options list
  )

  return(result)
}


#' Weighted Least Squares Error (WLSE) Objective Function
#'
#' Computes the total Weighted Least Squares Error (WLSE) objective function
#' value \eqn{f(W) = \sum_{i=1}^n \sum_{j=1}^n (w_i - a_{ij} w_j)^2} for a given
#' candidate priority weight vector \code{W} and Pairwise Reciprocal Matrix \code{A}
#' (Chu et al., 1979; Yuen, 2024).
#'
#' @param W Numeric vector of length \code{n} representing candidate priority weights.
#' @param A Numeric square matrix of dimension \code{n x n} representing a Pairwise Reciprocal Matrix (PRM).
#'
#' @return Numeric scalar. The total Weighted Least Squares Error value.
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Chu, A. T. W., Kalaba, R. E., & Spingarn, K. (1979). A comparison of two methods for determining the weights of belonging to fuzzy sets. \emph{Journal of Optimization Theory and Applications}, 27(4), 531-538.
#'
#' Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. \emph{arXiv preprint arXiv:2401.01190}. \url{https://doi.org/10.48550/arXiv.2401.01190}
#'
#' @export
#'
#' @examples
#' # Define a 3x3 Pairwise Reciprocal Matrix
#' A <- matrix(c(1, 2, 5,
#'               1/2, 1, 3,
#'               1/5, 1/3, 1), nrow = 3, byrow = TRUE)
#'
#' # Candidate priority weight vector
#' W <- c(0.58, 0.28, 0.14)
#'
#' # Calculate WLSE objective error
#' err <- WLSE(W, A)
#' print(err)
WLSE <- function(W, A) {
  # Extract matrix dimension n
  n <- nrow(A)

  # Construct n x n matrix of ones for vector broadcasting
  J <- matrix(1, n, n)

  # Compute component-wise WLS residual matrix and return total sum of squared residuals:
  # W * J creates matrix where entry (i,j) is w_i
  # t(W * J) creates matrix where entry (i,j) is w_j
  # A * t(W * J) creates matrix where entry (i,j) is a_ij * w_j
  # Matrix difference computes (a_ij * w_j - w_i), whose square equals (w_i - a_ij * w_j)^2
  # Summing across all entries yields total WLSE = sum_{i=1}^n sum_{j=1}^n (w_i - a_ij * w_j)^2
  sum((A * t(W * J) - W * J)^2)
}
