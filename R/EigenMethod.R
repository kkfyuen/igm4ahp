#' Check if Input Object is an Upper-Triangular Matrix (UTM) Ratings Vector
#'
#' Evaluates whether an input object \code{X} is a 1D atomic vector whose length \code{m}
#' corresponds to the exact upper-triangular element count \eqn{m = n(n-1)/2} of a square
#' Pairwise Reciprocal Matrix (PRM) of integer dimension \code{n}.
#'
#' @param X A vector, matrix, or array candidate object to be evaluated.
#'
#' @return A logical scalar: \code{TRUE} if \code{X} is a vector matching a valid
#'   upper-triangular pairwise comparison ratings length for an integer matrix order \code{n},
#'   or \code{FALSE} otherwise.
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. \emph{arXiv preprint arXiv:2401.01190}. \url{https://doi.org/10.48550/arXiv.2401.01190}
#'
#' @export
#'
#' @examples
#' # Length 3 vector corresponds to a 3x3 PRM upper triangle: 3*(3-1)/2 = 3
#' is.SizeOfUTM(c(2, 5, 3))     # TRUE
#'
#' # Length 4 vector does not correspond to an integer matrix order n
#' is.SizeOfUTM(c(1, 2, 3, 4))   # FALSE
is.SizeOfUTM <-function(X)
{
  # Default evaluation result initialized to FALSE
  Result=FALSE

  # Solve quadratic formula n^2 - n - 2*length(X) = 0 to derive matrix order: n = (1 + sqrt(1 + 8*m)) / 2
  n=(1+sqrt(1+8*length(X)))/2;

  # Verify if X is an atomic vector and the inferred matrix dimension n is a whole integer
  if(is.vector(X)&isWholeNo(n)) Result=TRUE;

  # Return logical evaluation status
  Result
}


#' Principal Eigenvector Priority Derivation (Saaty's EVM)
#'
#' Computes the priority weight vector \eqn{W} from a Pairwise Reciprocal Matrix (PRM)
#' using Saaty's Eigenvalue Method (EVM) by extracting and normalizing the principal
#' eigenvector corresponding to the maximum real eigenvalue \eqn{\lambda_{\max}} (Saaty, 1977, 1980).
#' Automatically converts an upper-triangular ratings vector into a full \code{n x n} PRM if supplied.
#'
#' @param PRM A numeric square matrix of dimension \code{n x n} representing a Pairwise Reciprocal Matrix,
#'   or an upper-triangular pairwise comparison ratings vector.
#'
#' @return A numeric vector of normalized priority weights \code{W} summing to unity.
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Saaty, T. L. (1977). A scaling method for priorities in hierarchical structures. \emph{Journal of Mathematical Psychology}, 15(3), 234-281.
#'
#' Saaty, T. L. (1980). \emph{The Analytic Hierarchy Process: Planning, Priority Setting, Resource Allocation}. McGraw-Hill, New York.
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
#' # Compute EVM priority weight vector
#' w_evm <- SaatyEigen(A)
#' print(w_evm)
SaatyEigen  <- function(PRM)
{
  # Check if input PRM is an upper-triangular ratings vector; convert to matrix if true
  if(is.SizeOfUTM(PRM))
  {
    PRM=getPRMFromRatings(PRM);
  }

  # Compute spectral decomposition (eigenvalues and eigenvectors) of PRM
  PRM.eigen=eigen(PRM);

  # Identify column index of the maximum real eigenvalue (lambda_max)
  PRM.eigen.max.index=which.max(Re(PRM.eigen$values));

  # Extract real component of the principal eigenvector corresponding to lambda_max
  PRM.top.vector=Re(PRM.eigen$vectors[,PRM.eigen.max.index]);

  # Normalize principal eigenvector to sum to unity (unit simplex constraint)
  W=PRM.top.vector/sum(PRM.top.vector);

  # Return normalized priority weight vector
  return(W)
}


#' Calculate Consistency Index (CI) for Pairwise Reciprocal Matrix
#'
#' Computes Saaty's Consistency Index \eqn{CI = (\lambda_{\max} - n) / (n - 1)} for an
#' \code{n x n} Pairwise Reciprocal Matrix (PRM), where \eqn{\lambda_{\max}} is the maximum
#' real eigenvalue and \eqn{n} is the matrix order (Saaty, 1977, 1980).
#'
#' @param PRM A numeric square matrix of dimension \code{n x n} representing a Pairwise Reciprocal Matrix.
#'
#' @return A numeric scalar representing the Consistency Index (CI), rounded to 10 decimal places.
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Saaty, T. L. (1977). A scaling method for priorities in hierarchical structures. \emph{Journal of Mathematical Psychology}, 15(3), 234-281.
#'
#' Saaty, T. L. (1980). \emph{The Analytic Hierarchy Process: Planning, Priority Setting, Resource Allocation}. McGraw-Hill, New York.
#'
#' Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. \emph{arXiv preprint arXiv:2401.01190}. \url{https://doi.org/10.48550/arXiv.2401.01190}
#'
#' @export
#'
#' @examples
#' A <- matrix(c(1, 2, 5,
#'               1/2, 1, 3,
#'               1/5, 1/3, 1), nrow = 3, byrow = TRUE)
#' ci <- ConsistencyIndex(A)
#' print(ci)
ConsistencyIndex <- function(PRM)
{
  # Extract dimension (number of evaluation criteria) n of the square PRM
  n=nrow(PRM);

  # Compute spectral decomposition of the PRM
  PRM.eigen=eigen(PRM);

  # Extract maximum real eigenvalue (lambda_max)
  PRM.eigen.max=max(Re(PRM.eigen$values));

  # Calculate Consistency Index: CI = (lambda_max - n) / (n - 1)
  CI=((PRM.eigen.max-n) /(n-1));

  # Return CI rounded to 10 decimal places
  return(round(CI,10))
}


#' Calculate Consistency Ratio (CR) for Pairwise Reciprocal Matrix
#'
#' Computes Saaty's Consistency Ratio \eqn{CR = CI / RI} for a Pairwise Reciprocal Matrix (PRM),
#' where \eqn{CI} is the Consistency Index and \eqn{RI} is the average Random Index for matrix
#' size \eqn{n \le 15} (Saaty, 1980, p. 21).
#'
#' @param PRM A numeric square matrix of dimension \code{n x n} representing a Pairwise Reciprocal Matrix.
#' @param randomIndex Optional numeric vector providing custom Random Index (RI) values for
#'   matrix dimensions \code{1} to \code{15}. Defaults to Saaty's standard RI table.
#'
#' @return A numeric scalar representing the calculated Consistency Ratio (CR), or \code{NA_real_}
#'   with a warning if \code{n > 15} without custom RI values.
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Saaty, T. L. (1980). \emph{The Analytic Hierarchy Process: Planning, Priority Setting, Resource Allocation}. McGraw-Hill, New York (p. 21).
#'
#' Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. \emph{arXiv preprint arXiv:2401.01190}. \url{https://doi.org/10.48550/arXiv.2401.01190}
#'
#' @export
#'
#' @examples
#' A <- matrix(c(1, 2, 5,
#'               1/2, 1, 3,
#'               1/5, 1/3, 1), nrow = 3, byrow = TRUE)
#' cr <- ConsistencyRatio(A)
#' print(cr)
ConsistencyRatio <- function(PRM,randomIndex=NULL)
{
  # Extract dimension n of the PRM
  n=nrow(PRM);

  # If custom Random Index vector is not supplied, default to Saaty's standard RI lookup table
  if(is.null(randomIndex))
  {
    # Return warning and NA if dimension n exceeds standard table threshold (n > 15)
    if(n>15) {
      warning("n should be less than and equal to 15");
      return(NA_real_)
    }
    # Saaty (1980, p. 21) Random Index (RI) table for matrix orders n = 1 to 15
    randomIndex = c(0, 0, 0.58, 0.90, 1.12, 1.24, 1.32, 1.41, 1.45, 1.49, 1.51, 1.48,
                    1.56, 1.57, 1.59);
  }

  # Compute Consistency Ratio: CR = CI / RI[n]
  CR=ConsistencyIndex(PRM)/ randomIndex[n];

  # Return calculated Consistency Ratio
  return(CR)
}
