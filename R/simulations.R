#' Check if a Value is a Whole Number
#'
#' Tests whether numerical values in \code{x} are whole numbers (integers) within
#' a specified floating-point machine precision tolerance threshold. Used to verify
#' that the matrix order \code{n} derived from upper-triangular pairwise ratings count
#' \eqn{m = n(n-1)/2} is a valid integer.
#'
#' @param x Numeric vector or scalar to be evaluated.
#' @param tolerance Numeric scalar specifying the numerical tolerance threshold
#'   (defaults to machine precision \code{sqrt(.Machine$double.eps)}).
#'
#' @return A logical vector or scalar of the same dimensions as \code{x}, indicating
#'   \code{TRUE} for whole numbers and \code{FALSE} otherwise.
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. \emph{arXiv preprint arXiv:2401.01190}. \url{https://doi.org/10.48550/arXiv.2401.01190}
#'
#' @export
#'
#' @examples
#' isWholeNo(4.0)  # TRUE
#' isWholeNo(4.1)  # FALSE
isWholeNo <-function(x, tolerance  = sqrt(.Machine$double.eps))
{
  # Evaluate if absolute difference between x and round(x) is within precision tolerance
  result=abs(x - round(x)) < tolerance;

  # Return logical evaluation result
  return(result)
}

#' Generate Random Pairwise Comparison Ratings Vector
#'
#' Generates a random upper-triangular pairwise comparison ratings vector sampled
#' from Saaty's 1/9 to 9 ratio scale (or a custom evaluation scale) for a Pairwise
#' Reciprocal Matrix (PRM) of dimension \code{n} (Section 4; Yuen, 2024).
#'
#' @param n Integer. Dimension of the Pairwise Reciprocal Matrix (PRM), representing
#'   the number of evaluation criteria. Defaults to a random integer sampled between 3 and 15.
#' @param scale Numeric vector. Set of discrete rating scale values from which pairwise
#'   comparisons are sampled. Defaults to Saaty's standard AHP scale \code{c(1/9, 1/8, ..., 1/2, 1, 2, ..., 9)}.
#'
#' @return A numeric vector of length \code{n(n-1)/2} containing randomly sampled pairwise comparison ratings.
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. \emph{arXiv preprint arXiv:2401.01190}. \url{https://doi.org/10.48550/arXiv.2401.01190}
#'
#' @export
#'
#' @examples
#' # Generate a random ratings vector for a 4x4 matrix
#' ratings <- randomPRR(n = 4)
#' print(ratings)
randomPRR <- function(n=sample(3:15,1),scale=c(1/(9:1),2:9))
{
  # Determine total number of scale elements available
  z = length(scale)

  # Calculate number of upper triangular pairwise comparisons: m = n*(n-1)/2
  m=n*(n-1)/2

  # Randomly sample m scale indices with replacement
  ramInd = sample(1:z,m, replace = TRUE)

  # Extract pairwise ratings matching sampled indices
  ratingVector = scale[ramInd]

  # Return upper-triangular ratings vector
  ratingVector
}


#' Simulation Verification for Inverse Gram Matrix (IGM) Methods
#'
#' Executes a Monte Carlo simulation over \code{size} randomly generated Pairwise Reciprocal
#' Matrices (PRMs) to verify exact numerical weight vector equivalence across closed-form Inverse Gram
#' Matrix variants—Pseudo-IGM (\code{\link{pseudoIGM}}), Element-form NIGM (\code{\link{eNIGM}}),
#' and Lagrangian-IGM (\code{\link{lagrangeIGM}})—logging verification metrics and Weighted Least
#' Squares Error (WLSE) values to a CSV file (Section 4.1; Yuen, 2024).
#'
#' @param size Integer. Total number of random PRM simulation iterations to execute.
#' @param seedNo Integer. Seed number for the pseudo-random number generator to ensure reproducibility (default: \code{100}).
#' @param errorFile Character string. Output CSV file path for logging simulation results and discrepancies (default: \code{"igmData.csv"}).
#'
#' @return Logical. Returns \code{FALSE} if all simulation iterations pass without discrepancy, or \code{TRUE} if a discrepancy is detected.
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. \emph{arXiv preprint arXiv:2401.01190}. \url{https://doi.org/10.48550/arXiv.2401.01190}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Run a small simulation of 10 random PRM instances
#'   runSimIGM(size = 10, seedNo = 123, errorFile = tempfile(fileext = ".csv"))
#' }
#' @importFrom MASS fractions
runSimIGM <-function(size,seedNo = 100,
                     errorFile = "igmData.csv")
{
  # Initialize CSV file header for simulation logging
  cat(
    "ID;",
    "Error;",
    "r;",
    "A;",
    "n;",
    "w;","WLSE","\n",
    file = errorFile)

  # Set pseudo-random seed for simulation reproducibility
  set.seed(seedNo)
  i = 1
  go = TRUE

  # Main simulation loop over requested sample size
  while (go) {
    print(paste0("i:",i))

    # Generate random ratings and construct Pairwise Reciprocal Matrix (PRM)
    rA = randomPRR()
    A =getPRMFromRatings(rA)
    #print(A)

    # Sample random non-zero scalar shift r in range [-1000, 1000]
    r = sample(c((-1000):(-1), 1:1000),1)

    # Compute priority weight vectors using Pseudo-IGM, eNIGM, and Lagrangian-IGM
    w_pseudoIGM  = round(pseudoIGM(A),8)
    w_eNIGM = round(eNIGM(A, r = r ),8)
    w_lagrangeIGM = lagrangeIGM(A)
    n1 = length(w_lagrangeIGM)
    w_lagrangeIGM = w_lagrangeIGM[-n1]
    w_rlagrangeIGM = lagrangeIGM(A,r=r)[-n1]

    #  wlsOpt = WLSopt(A)
    #   w_WLSopt = round(wlsOpt$solution,2)

    # Calculate absolute discrepancies between Pseudo-IGM and other IGM variants
    er = sum(abs(w_pseudoIGM - w_eNIGM)) +
      sum(abs(w_pseudoIGM - w_lagrangeIGM)) +
      sum(abs(w_pseudoIGM - w_rlagrangeIGM))

    # +
    #   sum(abs(round(w_pseudoIGM,2) - w_WLSopt))
    er = round(er, 6)
    # print(er)

    # Flag error if total discrepancy deviates from zero
    error = (er != 0) # there is an error

    if(!error)
    {
      # print(A)
      # print(w_pseudoIGM)

      # Compute Weighted Least Squares Error (WLSE)
      wlse = round(WLSE(w_pseudoIGM,A),4)

      # Log successful iteration results to CSV output file
      cat(
        i,";",
        er,";",
        r,";",
        paste(MASS::fractions(rA),collapse = ","),";",
        length(w_pseudoIGM),";",
        paste(round(w_pseudoIGM,4),collapse =","),";",
        wlse,"\n",
        file = errorFile, append = TRUE)
    } else
    {
      # Halt execution on discrepancy detection and output diagnostic report
      if (error) go = FALSE
      print("Error is found.")
      cat("-----------------------------------\n",
          file = errorFile, append = TRUE)

      cat("An error is found. Stop.\n",
          file = errorFile, append = TRUE)

      cat("rA",paste(MASS::fractions(rA),collapse = ","),".\n",
          "r:",r,".\n",
          "w_pseudoIGM:", w_pseudoIGM,".\n",
          "w_eNIGM:", w_eNIGM,".\n",
          "w_lagrangeIGM:", w_lagrangeIGM,".\n",
          # "w_WLSopt:", w_WLSopt,".\n",
          # "wlse_WLSopt:",WLSE(w_WLSopt,A),".\n",
          # "wlse_WLSopt_obj:",wlsOpt$objective,".\n",
          "wlse_pseudoIGM:",WLSE(w_pseudoIGM,A),".\n",
          file = errorFile,append = TRUE)

      cat("-----------------------------------\n",
          file = errorFile, append = TRUE)

    }

    i = i + 1
    if (i>size) go = FALSE


  }

  #if true/ same, 1; if False, 0.
  # easier to check how many cases are 1 and for debug
  print("end")
  error
}


#' Simulation Verification Comparing Pseudo-IGM and Non-Linear WLS Optimization
#'
#' Executes a Monte Carlo simulation benchmarking the closed-form priority weight vectors
#' derived by Pseudo-IGM (\code{\link{pseudoIGM}}) against the numerical global optimization
#' solver (\code{\link{WLSopt}}) via \code{nloptr}, verifying that closed-form IGM formulas
#' consistently achieve the exact global minimum Weighted Least Squares Error (WLSE) without
#' relying on iterative numerical solvers (Section 4.2; Yuen, 2024).
#'
#' @param size Integer. Total number of random PRM simulation instances to evaluate (default: \code{10}).
#' @param seedNo Integer. Seed number for pseudo-random number generator initialization to ensure reproducibility (default: \code{100}).
#' @param sigDigits Integer. Number of significant decimal places used for rounding priority weight vectors during equality verification (default: \code{4}).
#' @param errorFile Character string. Output text file path for logging simulation results and discrepancy metrics (default: \code{"wlseOptData.txt"}).
#'
#' @return Logical. Returns \code{FALSE} if all simulation iterations pass without discrepancy, or \code{TRUE} if a discrepancy is detected.
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. \emph{arXiv preprint arXiv:2401.01190}. \url{https://doi.org/10.48550/arXiv.2401.01190}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Run simulation comparing Pseudo-IGM against numerical solver for 5 random PRMs
#'   runWLSEOpt(size = 5, seedNo = 123, errorFile = tempfile(fileext = ".txt"))
#' }
#' @importFrom MASS fractions
runWLSEOpt <-function(size=10,seedNo = 100, sigDigits = 4,
                      errorFile = "wlseOptData.txt")
{
  # Initialize text file header for simulation data logging
  cat(
    "ID;",
    "Error;",
    "A;",
    "n;",
    "w_PIGM;",
    "w_opt;",
    "wlse_PIGM;",
    "wlse_opt;",
    "wlse_sig;",
    "Sum(w_opt)",
    "\n",
    file = errorFile)

  # Set pseudo-random seed for simulation reproducibility
  set.seed(seedNo)
  i = 1
  go = TRUE

  # Main simulation loop over requested sample size
  while (go) {
    print(paste0("i:",i))

    # Generate random ratings and construct Pairwise Reciprocal Matrix (PRM)
    rA = randomPRR()
    A =getPRMFromRatings(rA)

    # Calculate exact closed-form priority weight vector using Pseudo-IGM
    w_pseudoIGM  = pseudoIGM(A)
    # w_eNIGM = round(eNIGM(A, r = r ),8)
    # w_lagrangeIGM = lagrangeIGM(A)
    # n1 = length(w_lagrangeIGM)
    # w_lagrangeIGM = w_lagrangeIGM[-n1]
    # w_rlagrangeIGM = lagrangeIGM(A,r=r)[-n1]

    # Calculate numerical optimization solution vector via nloptr WLS solver
    wlsOpt = WLSopt(A)
    w_WLSopt = wlsOpt$solution

    # Calculate rounded absolute discrepancy between Pseudo-IGM and numerical optimization solutions
    er =  sum(abs(round(w_pseudoIGM,sigDigits) - round(w_WLSopt,sigDigits)))
    er = round(er, sigDigits-1)
    # print(er)

    # Flag error if discrepancy deviates from zero tolerance
    error = (er != 0) # there is an error

    if(!error)
    {
      # print(A)
      # print(w_pseudoIGM)

      # Evaluate Weighted Least Squares Error (WLSE) for both weight vectors
      wlse_PIGM = WLSE(w_pseudoIGM,A)
      wlse_opt = WLSE(w_WLSopt,A)

      # Log successful benchmarking iteration to output file
      cat(
        i,";",
        er,";",
        paste(MASS::fractions(rA),collapse = ","),";",
        length(w_pseudoIGM),";",
        paste(round(w_pseudoIGM,sigDigits),collapse =","),";",
        paste(round(w_WLSopt,sigDigits),collapse =","),";",
        wlse_PIGM,";",
        wlse_opt,";",
        round(wlse_PIGM,sigDigits),";",
        sum(w_WLSopt),"\n",
        file = errorFile, append = TRUE)
    } else
    {
      # Halt simulation loop on discrepancy detection and append detailed error log
      if (error) go = FALSE
      print("Error is found")
      cat("-----------------------------------\n",
          file = errorFile, append = TRUE)

      cat("An error is found. Stop.\n",
          file = errorFile, append = TRUE)

      cat(
        "er",er,".\n",
        "rA",paste(MASS::fractions(rA),collapse = ","),".\n",
        "w_pseudoIGM:", w_pseudoIGM,".\n",
        "w_WLSopt:", w_WLSopt,".\n",
        "wlse_WLSopt:",WLSE(w_WLSopt,A),".\n",
        "wlse_WLSopt_obj:",wlsOpt$objective,".\n",
        "wlse_pseudoIGM:",WLSE(w_pseudoIGM,A),".\n",
        "sum(w_WLSopt):",sum(w_WLSopt),"\n",
        file = errorFile,append = TRUE)

      cat("-----------------------------------\n",
          file = errorFile, append = TRUE)

    }

    i = i + 1
    if (i>size) go = FALSE


  }

  #if true/ same, 1; if False, 0.
  # easier to check how many cases are 1 and for debug
  print("end")
  error
}






#' Construct Pairwise Reciprocal Matrix from Ratings Vector
#'
#' Reconstructs an \code{n x n} Pairwise Reciprocal Matrix (PRM) from an upper-triangular
#' pairwise ratings vector \code{V} of length \code{n(n-1)/2}.
#'
#' @param V Numeric vector containing upper-triangular pairwise comparison ratings of length \code{n(n-1)/2}.
#' @param n Numeric scalar representing the dimension of the square PRM.
#'   Calculated automatically from the length of \code{V} by solving \code{n(n-1)/2 = length(V)}.
#'
#' @return A numeric square matrix of dimension \code{n x n} representing the Pairwise Reciprocal Matrix (PRM),
#'   or \code{NA_real_} with a warning if the input vector length is invalid.
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. \emph{arXiv preprint arXiv:2401.01190}. \url{https://doi.org/10.48550/arXiv.2401.01190}
#'
#' @export
#'
#' @examples
#' # Vector of 3 upper-triangular ratings for a 3x3 matrix
#' V <- c(2, 5, 3)
#' A <- getPRMFromRatings(V)
#' print(A)
getPRMFromRatings <- function(V,n=(1+sqrt(1+8*length(V)))/2)
{
  # Check if inferred dimension n is a whole number
  if(!isWholeNo(n))
  {
    warning("Not Valid: length(V) == (n*(n-1)/2): return NA")
    return(NA_real_)
  }

  # Validate that length of vector V equals n*(n-1)/2
  if (!(length(V) == (n*(n-1)/2)))
  {
    warning("Not Valid: length(V) == (n*(n-1)/2): return NA")
    return(NA_real_)
  }

  # Compute reciprocal ratings vector RV = 1/V
  RV=1/V;

  # Initialize n x n identity matrix for the PRM
  RPM = diag(x = 1, n, n);
  k=0;

  # Populate upper-triangular with V and lower-triangular with reciprocal values RV
  for(i in 1:(n-1))
    for(j in (i+1):n)
    {
      k=k+1;
      RPM[i,j] = V[k];
      RPM[j,i] = RV[k];
    }

  # Return constructed Pairwise Reciprocal Matrix
  return(RPM)
}






# startTime =Sys.time()
# runSimIGM(1000000,seedNo = 2000)
# fileIGMData = "igmData.Rdata"
# endTime = Sys.time()
# timeSlot = endTime - startTime
# save(startTime,endTime,timeSlot,file = fileIGMData)



#' Run and Benchmark Inverse Gram Matrix Simulation
#'
#' Wraps the \code{\link{runSimIGM}} simulation function with automated execution timing,
#' recording the start time, end time, and total elapsed duration, and returning/saving
#' a consolidated list containing simulation outcome metrics and timing metadata.
#'
#' @param size Integer. Total number of random PRM simulation iterations to execute (default: \code{1000000}).
#' @param seedNo Integer. Seed number for pseudo-random number generator initialization (default: \code{2000}).
#' @param fileIGMData Character string. Output RData file path to save the timing and result metadata (default: \code{"igmData.Rdata"}).
#' @param errorFile Character string. CSV file path for logging simulation results and discrepancies (default: \code{"igmData.csv"}).
#'
#' @return A named \code{list} containing:
#'   \item{error}{Logical status flag returned by \code{\link{runSimIGM}} (\code{FALSE} if no errors found, \code{TRUE} if discrepancy detected).}
#'   \item{startTime}{\code{POSIXct} timestamp recorded when execution started.}
#'   \item{endTime}{\code{POSIXct} timestamp recorded when execution completed.}
#'   \item{timeSlot}{\code{difftime} object representing the total elapsed execution time.}
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. \emph{arXiv preprint arXiv:2401.01190}. \url{https://doi.org/10.48550/arXiv.2401.01190}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Run timed simulation and access timing metrics from result list
#'   res <- runSimWithTiming(size = 10, seedNo = 2000, fileIGMData = tempfile(fileext = ".Rdata"))
#'   print(res$timeSlot)
#' }
runSimWithTiming <- function(size = 1000000, seedNo = 2000,
                             fileIGMData = "igmData.Rdata",
                             errorFile = "igmData.csv")
{
  # Record simulation start time
  startTime = Sys.time()

  # Execute Inverse Gram Matrix Monte Carlo simulation
  error = runSimIGM(size = size, seedNo = seedNo, errorFile = errorFile)

  # Record simulation end time and compute total elapsed duration
  endTime = Sys.time()
  timeSlot = endTime - startTime

  # Construct structured result list containing simulation status and timing metadata
  result = list(error = error,
                startTime = startTime,
                endTime = endTime,
                timeSlot = timeSlot)

  # Save result list and individual timing variables to RData file
  save(result, startTime, endTime, timeSlot, file = fileIGMData)

  # Return structured result list
  return(result)
}


# startTime =Sys.time()
# runWLSEOpt(10000,seedNo = 2000)
# filewlseOptData = "wlseOptData.Rdata"
# endTime = Sys.time()
# timeSlot = endTime - startTime
# save(startTime,endTime,timeSlot,file = filewlseOptData)


#' Run and Benchmark Non-Linear WLS Optimization Simulation
#'
#' Wraps the \code{\link{runWLSEOpt}} simulation function with automated execution timing,
#' recording the start time, end time, and total elapsed duration, and returning/saving
#' a consolidated list containing simulation outcome metrics and timing metadata.
#'
#' @param size Integer. Total number of random PRM simulation instances to evaluate (default: \code{10000}).
#' @param seedNo Integer. Seed number for pseudo-random number generator initialization (default: \code{2000}).
#' @param filewlseOptData Character string. Output RData file path to save timing and result metadata (default: \code{"wlseOptData.Rdata"}).
#' @param sigDigits Integer. Number of significant decimal places used for rounding priority weight vectors during equality verification (default: \code{4}).
#' @param errorFile Character string. Output text file path for logging simulation results and discrepancy metrics (default: \code{"wlseOptData.txt"}).
#'
#' @return A named \code{list} containing:
#'   \item{error}{Logical status flag returned by \code{\link{runWLSEOpt}} (\code{FALSE} if no errors found, \code{TRUE} if discrepancy detected).}
#'   \item{startTime}{\code{POSIXct} timestamp recorded when execution started.}
#'   \item{endTime}{\code{POSIXct} timestamp recorded when execution completed.}
#'   \item{timeSlot}{\code{difftime} object representing total elapsed execution duration.}
#'
#' @author Kevin Kam Fung Yuen
#' @references
#' Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy Process: Explainability of Weighted Least Squares Optimization Method. \emph{arXiv preprint arXiv:2401.01190}. \url{https://doi.org/10.48550/arXiv.2401.01190}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Run timed WLS optimization benchmark simulation
#'   tmp <- tempfile(fileext = ".Rdata")
#'   res <- runWLSEOptWithTiming(
#'     size = 100,
#'     seedNo = 2000,
#'     filewlseOptData = tmp
#'   )
#'   print(res$timeSlot)
#' }
runWLSEOptWithTiming <- function(size = 10000, seedNo = 2000,
                                 filewlseOptData = "wlseOptData.Rdata",
                                 sigDigits = 4,
                                 errorFile = "wlseOptData.txt")
{
  # Record simulation start time
  startTime = Sys.time()

  # Execute WLS optimization Monte Carlo simulation
  error = runWLSEOpt(size = size, seedNo = seedNo, sigDigits = sigDigits, errorFile = errorFile)

  # Record simulation end time and compute total elapsed duration
  endTime = Sys.time()
  timeSlot = endTime - startTime

  # Construct structured result list containing simulation status and timing metadata
  result = list(error = error,
                startTime = startTime,
                endTime = endTime,
                timeSlot = timeSlot)

  # Save result list and individual timing variables to RData file
  save(result, startTime, endTime, timeSlot, file = filewlseOptData)

  # Return structured result list
  return(result)
}
