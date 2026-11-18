#' Construct Full Design (Coefficient) Matrix for Pairwise Comparisons
#'
#' Converts an \code{n x n} Pairwise Reciprocal Matrix (PRM) \code{A} into the full
#' design coefficient matrix \code{D} of dimension \code{(n(n-1) + 1) x n} via Algorithm 1
#' (Yuen, 2024). Stacks all \code{n(n-1)} pairwise weighted differential equations
#' (\eqn{w_i - a_{ij} w_j = 0}) and appends the unit simplex normalization row vector
#' (\eqn{e = [1, \dots, 1]}, representing \eqn{\sum_{i=1}^n w_i = 1}).
#'
#' @param A Numeric square matrix of dimension \code{n x n} representing a Pairwise Reciprocal Matrix (PRM).
#'
#' @return A numeric design matrix of dimension \code{(n(n-1) + 1) x n} with column names \code{w1, ..., wn}.
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. \emph{arXiv preprint arXiv:2401.01190}. \url{https://doi.org/10.48550/arXiv.2401.01190}
#'
#' @export
#'
#' @examples
#' A <- matrix(c(1, 2, 1/2, 1), nrow = 2)
#' D <- fullCoefficientMatrix(A)
#' print(D)
fullCoefficientMatrix <- function(A)
{
  # maybe use ratings score in upper instead of A in another version
  n = nrow(A)
  # ME = rep(1,n)  # equation matrix
  ME = matrix(0, n * (n - 1), n)
  k = 1
  for (i in 1:(n - 1)){
    for (j in (i + 1):n){
      # upper triangular ones: w_i - a_ij * w_j = 0
      ME[k, i] = 1
      ME[k, j] = -A[i, j]
      k = k + 1
      # lower triangular ones: -a_ji * w_i + w_j = 0 (using a_ji = 1/a_ij)
      ME[k, i] = -(1 / A[i, j])
      ME[k, j] = 1
      k = k + 1

    }  }

  # Append unit simplex normalization row vector of ones (e = [1, ..., 1])
  ME = rbind(ME, rep(1, n))
  colnames(ME) <- paste0("w", 1:n)
  ME
}



#' Pseudo-Normalized Inverse Gram Matrix (PNIGM) Prioritization
#'
#' Computes the priority weight vector \eqn{W} from a Pairwise Reciprocal Matrix (PRM) \code{A}
#' using the Pseudo Inverse Gram Matrix (PIGM) closed-form solution for Weighted Least Squares
#' (WLS) prioritization (Theorem 1 & Lemma 1; Yuen, 2024).
#'
#' @param A Numeric square matrix of dimension \code{n x n} representing a Pairwise Reciprocal Matrix (PRM).
#' @param detail Logical. If \code{TRUE}, returns a list containing intermediate matrices
#'   (\code{D}, \code{G}, \code{Ginv}), unnormalized row sum vector \code{V}, and final priority vector \code{W}
#'   (default: \code{FALSE}).
#'
#' @return If \code{detail = FALSE}, a numeric vector of normalized priority weights \code{W}.
#'   If \code{detail = TRUE}, a list of intermediate calculation steps and matrices.
#'
#' @author Kevin Kam Fung Yuen
#' @references
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
#' # Compute priority weights via PNIGM
#' W <- PNIGM(A)
#' print(W)
PNIGM <- function(A, detail = FALSE)
{

  # Construct full design coefficient matrix D via Algorithm 1
  D = fullCoefficientMatrix(A)

  # Compute Gram Matrix G = D^T * D (Symmetric Positive-Definite; Lemma 1)
  G = t(D) %*% D

  # Invert Gram Matrix G^(-1)
  Ginv = solve(G)

  # Calculate unnormalized weight vector V via row sums of G^(-1): V = e * G^(-1)
  V = rowSums(Ginv)

  # Normalize priority weights to sum to unity: W = V / sum(V)
  W = V / sum(V)
  result = W

  if (detail)
  {
    result = list(D = D, G = G, Ginv = Ginv, V = V, W = W)

  }
  result
}


#' Pseudo Inverse Gram Matrix (PIGM) Prioritization
#'
#' Computes the priority weight vector \eqn{W} from a Pairwise Reciprocal Matrix (PRM) \code{A}
#' using the closed-form Pseudo Inverse Gram Matrix (PIGM) formula \eqn{W = \frac{e G^{-1}}{e G^{-1} e^T}}
#' where \eqn{G = D^T D} (Theorem 1; Yuen, 2024).
#'
#' @param A Numeric square matrix of dimension \code{n x n} representing a Pairwise Reciprocal Matrix (PRM).
#'
#' @return A numeric vector of normalized priority weights \code{W}.
#'
#' @author Kevin Kam Fung Yuen
#' @references
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
#' # Compute priority weights via pseudoIGM
#' W <- pseudoIGM(A)
#' print(W)
pseudoIGM <- function(A)
{
  # Extract matrix dimension n
  n = ncol(A)

  # Construct full design coefficient matrix D via Algorithm 1
  D = fullCoefficientMatrix(A)

  # Compute Gram Matrix G = D^T * D
  G = t(D) %*% (D)

  # Invert Gram Matrix G^(-1)
  Ginv = solve(G)

  #q = c(rep(0,n*(n-1)),1)
  #V = Ginv %*% t(D) %*% q  # same as Moore-Penrose pseudoinverse D+ = G^(-1) * D^T

  # Construct vector of ones (e = [1, ..., 1])
  ones = rep(1,n)

  # Calculate unnormalized weight vector V via vector-matrix multiplication: ones %*% Ginv
  V =   ones %*% Ginv

  # Normalize priority weights to sum to unity: W = V / sum(V) (Theorem 1, Eq. 23)
  W = V / sum(V)

  # Return final priority weight vector as a plain numeric vector
  as.vector(W)
}


#' Pseudo Inverse Gram Matrix (PIGM) Prioritization
#'
#' Computes the priority weight vector \eqn{W} from a Pairwise Reciprocal Matrix (PRM) \code{A}
#' by calling \code{\link{pseudoIGM}}.
#'
#' @param A Numeric square matrix of dimension \code{n x n} representing a Pairwise Reciprocal Matrix (PRM).
#'
#' @return A numeric vector of normalized priority weights \code{W}.
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. \emph{arXiv preprint arXiv:2401.01190}. \url{https://doi.org/10.48550/arXiv.2401.01190}
#'
#' @export
#'
#' @examples
#' A <- matrix(c(1, 2, 5,
#'               1/2, 1, 3,
#'               1/5, 1/3, 1), nrow = 3, byrow = TRUE)
#' W <- PIGM(A)
#' print(W)
PIGM <- function(A)
{
  pseudoIGM(A)
}


#' Element-wise Gram Matrix Construction
#'
#' Constructs the \code{n x n} Gram matrix \eqn{G = (g_{ij})} directly from element-wise formulas
#' derived from the Weighted Least Squares (WLS) objective function (Theorem 2; Yuen, 2024),
#' bypassing the construction of design matrix \code{D}, with an optional scalar shift parameter \code{r}.
#'
#' @param A Numeric square matrix of dimension \code{n x n} representing a Pairwise Reciprocal Matrix (PRM).
#' @param r Numeric scalar shift parameter added to all matrix entries (default: \code{0}).
#'
#' @return A numeric square matrix of dimension \code{n x n} representing the constructed Gram matrix.
#'
#' @author Kevin Kam Fung Yuen
#' @references
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
#' # Compute element-form Gram matrix (r = 0)
#' G <- eleGram(A, r = 0)
#' print(G)
#element form of Gram matrix based on objective function if r = 0
eleGram <- function(A,r=0)
{
  # Extract dimension n of the PRM
  n = nrow(A)

  # Initialize n x n matrix g with NA values
  g = matrix(NA,n,n)

  # Populate diagonal and off-diagonal elements of the Gram matrix (Theorem 2, Eq. 33)
  for (i in 1:n)
    for (j in 1:n)
    {
      if ( i == j)
      {
        # Diagonal elements: g_ii = (n - 1) + sum_{k != i} a_{ki}^2
        g[i,j] = (n-1) + sum(A[-j,i]^2)

      } else
      {
        # Off-diagonal elements: g_ij = -(a_ij + a_ji)
        g[i,j] = -(A[i,j]+A[j,i])
      }
    }

  # Return Gram matrix shifted by scalar r (Eq. 49)
  g + r
}


#' Element-form Normalized Inverse Gram Matrix (eNIGM) Prioritization
#'
#' Computes the priority weight vector \eqn{W} from a Pairwise Reciprocal Matrix (PRM) \code{A}
#' using the element-wise Normalized Inverse Gram Matrix (NIGM) closed-form solution (Theorem 3 & Remark 1; Yuen, 2024).
#' Parameter \code{r != 0} (e.g., \code{r = 1}) guarantees invertibility even for perfectly consistent PRMs.
#'
#' @param A Numeric square matrix of dimension \code{n x n} representing a Pairwise Reciprocal Matrix (PRM).
#' @param r Numeric scalar parameter for the Gram matrix shift (default: \code{1}).
#' @param detail Logical. If \code{TRUE}, returns a list containing intermediate results
#'   (\code{W}, \code{r}, \code{G}, \code{G_inv}, \code{rs_G_inv}) (default: \code{FALSE}).
#'
#' @return If \code{detail = FALSE}, a numeric vector of normalized priority weights \code{W}.
#'   If \code{detail = TRUE}, a list of intermediate calculation steps and matrices.
#'
#' @author Kevin Kam Fung Yuen
#' @references
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
#' # Compute priority weights via eNIGM (r = 1)
#' W <- eNIGM(A, r = 1)
#' print(W)
# element form of NIGM
eNIGM <-function(A, r=1, detail = FALSE)
{
  # Construct element-form Gram matrix G with regularization parameter r (Theorem 2 & 3)
  G = eleGram(A,r)

  # Compute matrix inverse G^(-1) = (G + r*E)^(-1)
  G_inv = solve(G)

  # Compute row sums (unnormalized priority vector) of G^(-1)
  rs_G_inv = rowSums(G_inv) # or colSums

  # Normalize priority weight vector to sum to unity (Eq. 42; parameter r cancels out during division)
  result = rs_G_inv / sum(rs_G_inv)

  if(detail)
  {
    result = list (W = result,r=r, G = G,
                   G_inv = G_inv,
                   rs_G_inv = rs_G_inv)
  }
  result
}

#' Construct Parametrized Design Matrix for NIGM
#'
#' Converts an \code{n x n} Pairwise Reciprocal Matrix (PRM) \code{A} into a design
#' coefficient matrix \code{D} of dimension \code{(n(n-1) + 1) x n}, incorporating
#' pairwise comparison differential equations and appending a parametrized shift row
#' vector based on scalar parameter \code{r} (Section 3.2, Theorem 3; Yuen, 2024).
#'
#' @param A Numeric square matrix of dimension \code{n x n} representing a Pairwise Reciprocal Matrix (PRM).
#' @param r Numeric scalar parameter determining the weight shift factor appended to the design matrix.
#'
#' @return A numeric design matrix of dimension \code{(n(n-1) + 1) x n} with column names \code{w1, ..., wn}.
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. \emph{arXiv preprint arXiv:2401.01190}. \url{https://doi.org/10.48550/arXiv.2401.01190}
#'
#' @export
#'
#' @examples
#' A <- matrix(c(1, 2, 1/2, 1), nrow = 2)
#' D <- rDesignMatrix(A, r = 1)
#' print(D)
rDesignMatrix <- function(A,r)
{
  # maybe use ratings score in upper instead of A in another version
  n = nrow(A)
  #ME = rep(1,n)  #equation matrix
  ME = matrix(0,n*(n-1),n)
  k=1
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      # upper triangular ones
      ME[k,i] =1
      ME[k,j] = -A[i,j]
      k=k+1
      # lower triangular ones
      ME[k,i]= -(1/A[i,j])
      ME[k,j]= 1
      k=k+1
    }  }

  # Append parameterized row vector based on scalar r (Eqs. 47-48)
  if (r>= 0)
  {
    ME = rbind(ME,rep(sqrt(r),n))
  } else
  {
    ME = rbind(ME,rep(-sqrt(abs(r)),n))
  }
  colnames(ME) <- paste0("w",1:n)
  ME
}


#' Matrix-form Normalized Inverse Gram Matrix (NIGM) Prioritization
#'
#' Computes the priority weight vector \eqn{W} from a Pairwise Reciprocal Matrix (PRM) \code{A}
#' using the matrix formulation of the Normalized Inverse Gram Matrix (NIGM) closed-form solution
#' (Theorem 3; Yuen, 2024).
#'
#' @param A Numeric square matrix of dimension \code{n x n} representing a Pairwise Reciprocal Matrix (PRM).
#' @param r Numeric scalar parameter for the Gram matrix shift (default: \code{1}).
#' @param detail Logical. If \code{TRUE}, returns a list containing intermediate matrices
#'   (\code{W}, \code{r}, \code{G}, \code{G_inv}, \code{rs_G_inv}) (default: \code{FALSE}).
#'
#' @return If \code{detail = FALSE}, a numeric vector of normalized priority weights \code{W}.
#'   If \code{detail = TRUE}, a list of intermediate calculation steps and matrices.
#'
#' @author Kevin Kam Fung Yuen
#' @references
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
#' # Compute priority weights via matrix NIGM (r = 1)
#' W <- NIGM(A, r = 1)
#' print(W)
# matrix form
NIGM <-function(A, r=1, detail = FALSE)
{
  # Construct parametrized design matrix D
  D = rDesignMatrix(A,r)

  # Compute Gram Matrix G = D^T * D depending on sign of parameter r (Eqs. 47-48)
  if(r>=0)
  {
    G = t(D) %*% D
  } else
  {
    nD = nrow(D)
    eD = D
    eD[nD,] = -D[nD,]
    G = t(D) %*% eD
  }

  # Compute matrix inverse G^(-1)
  G_inv = solve(G)

  # Calculate unnormalized weight vector via row sums of G^(-1)
  rs_G_inv = rowSums(G_inv) # or colSums

  # Normalize priority weights to sum to unity: W = rs_G_inv / sum(rs_G_inv)
  result = rs_G_inv / sum(rs_G_inv)

  if(detail)
  {
    result = list (W = result,r=r, G = G,
                   G_inv = G_inv,
                   rs_G_inv = rs_G_inv)
  }
  result
}


#' Lagrangian Inverse Gram Matrix (LIGM) Prioritization
#'
#' Computes the priority weight vector \eqn{W} and exact Lagrange multiplier \eqn{\lambda} from a
#' Pairwise Reciprocal Matrix (PRM) \code{A} using the augmented Lagrangian Inverse Gram Matrix
#' (LIGM) closed-form KKT system solution \eqn{[W, \lambda] = q \bar{G}^{-1}} (Theorem 4; Yuen, 2024),
#' eliminating the requirement for post-hoc vector normalization.
#'
#' @param A Numeric square matrix of dimension \code{n x n} representing a Pairwise Reciprocal Matrix (PRM).
#' @param r Numeric scalar parameter for the Gram matrix shift (default: \code{0}).
#' @param detail Logical. If \code{TRUE}, returns a list containing intermediate matrices
#'   (\code{W}, \code{G}, \code{G_inv}) (default: \code{FALSE}).
#'
#' @return If \code{detail = FALSE}, a numeric vector of length \code{n + 1} containing the
#'   priority weights \eqn{W} and the Lagrange multiplier \eqn{\lambda}. If \code{detail = TRUE}, a list
#'   containing intermediate calculation steps and matrices.
#'
#' @author Kevin Kam Fung Yuen
#' @references
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
#' # Compute priority weights via lagrangeIGM (r = 0)
#' W <- lagrangeIGM(A, r = 0)
#' print(W)
lagrangeIGM <- function(A,r=0,detail = FALSE)
{
  # Extract dimension n of the PRM
  n = nrow(A)

  # Initialize augmented (n+1) x (n+1) block matrix with NA values
  g = matrix(NA,n+1,n+1)

  # Construct the core n x n Gram matrix elements (Theorem 2, Eq. 33)
  for (i in 1:n)
    for (j in 1:n)
    {
      if ( i == j)
      {
        # Diagonal elements: g_ii = (n - 1) + sum_{k != i} a_{ki}^2
        g[i,j] = (n-1) + sum(A[-j,i]^2)

      } else
      {
        # Off-diagonal elements: g_ij = -(a_ij + a_ji)
        g[i,j] = -(A[i,j]+A[j,i])
      }
    }

  # Apply scalar shift r to the core n x n Gram block (Eq. 72)
  g = g + r

  # Set unit simplex equality constraint borders: row n+1 and column n+1 (Theorem 4, Eq. 63)
  g[n+1,] = g[,n+1] = 1
  g[n+1,n+1] = 0

  # Compute matrix inverse of the augmented KKT Gram matrix G^(-1)
  G_inv = solve(g)

  # Compute unnormalized row sums representing priority vector and Lagrange multiplier: [W, lambda] = q * G^(-1)
  result = rowSums(solve(g))

  if(detail)
  {
    result = list (W = result, G = g,
                   G_inv = G_inv)
  }
  result
}

#' @rdname lagrangeIGM
#' @export
LIGM <- function(A, r = 0, detail = FALSE) {
  lagrangeIGM(A = A, r = r, detail = detail)
}
