## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

## ----load-package-------------------------------------------------------------
library(igm4ahp)

## ----single-instance-demo-----------------------------------------------------
# Set seed for reproducible single-instance demonstration
set.seed(2024)

# Generate a random 5x5 Pairwise Reciprocal Matrix
n_demo <- 5
rA_demo <- randomPRR(n = n_demo)
A_demo <- getPRMFromRatings(rA_demo)

print("--- Generated 5x5 Pairwise Reciprocal Matrix A ---")
print(round(A_demo, 3))

# 1. Calculate exact closed-form priorities via PIGM
w_PIGM <- PIGM(A_demo)

# 2. Calculate numerical optimization priorities via WLSopt (nloptr ISRES solver)
# Note: Lower maxeval used here for quick vignette compilation
opt_res <- WLSopt(A_demo, maxeval = 50000)
w_WLSopt <- opt_res$solution

# 3. Calculate Weighted Least Squares Error (WLSE) for both weight vectors
wlse_PIGM <- WLSE(w_PIGM, A_demo)
wlse_WLSopt <- WLSE(w_WLSopt, A_demo)

# Assemble comparative summary table (rounded to 4 decimal places, eps = 4)
df_comp <- data.frame(
  PIGM_ClosedForm = round(w_PIGM, 4),
  WLSopt_Numerical = round(w_WLSopt, 4),
  Abs_Difference = round(abs(w_PIGM - w_WLSopt), 4)
)
rownames(df_comp) <- paste0("w", 1:n_demo)

print("--- Priority Weight Comparison (4 Decimal Places) ---")
print(df_comp)

cat("\n--- Objective Function Error Comparison f(W) ---\n")
cat("PIGM Closed-Form WLSE :", round(wlse_PIGM, 6), "\n")
cat("WLSopt Numerical WLSE :", round(wlse_WLSopt, 6), "\n")
cat("Residual Difference   :", round(abs(wlse_PIGM - wlse_WLSopt), 6), "\n")

## ----run-sim-wlseopt----------------------------------------------------------
# Ensure local ./data directory exists for saving output logs
data_dir <- "./data"
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

txt_file <- file.path(data_dir, "wlseOptData.txt")

# Execute Monte Carlo Simulation 2 for N = 10 sample iterations
# Matching precision threshold eps = 4 decimal places
sim_error <- runWLSEOpt(
  size = 10,
  seedNo = 2000,
  sigDigits = 4,
  errorFile = txt_file
)

cat("Simulation Discrepancy Status (FALSE = No errors detected):", sim_error, "\n")

## ----inspect-wlseopt-log------------------------------------------------------
# Read simulation text log
sim_data <- read.csv(txt_file, sep = ";", stringsAsFactors = FALSE)
colnames(sim_data) <- trimws(colnames(sim_data))

# Display simulation metrics table
knitr::kable(
  sim_data,
  caption = "Simulation 2 Benchmarking Log Output (./data/wlseOptData.txt)"
)

## ----run-sim-wlseopt-timing---------------------------------------------------
rdata_file <- file.path(data_dir, "wlseOptData.Rdata")

# Execute timed Simulation 2 for N = 5 sample iterations
timing_results <- runWLSEOptWithTiming(
  size = 5,
  seedNo = 2024,
  filewlseOptData = rdata_file,
  sigDigits = 4,
  errorFile = txt_file
)

cat("--- Simulation 2 Execution Timing Summary ---\n")
cat("Start Time       :", format(timing_results$startTime), "\n")
cat("End Time         :", format(timing_results$endTime), "\n")
cat("Elapsed Duration :", round(as.numeric(timing_results$timeSlot, units = "secs"), 3), "seconds\n")
cat("Discrepancy Status:", timing_results$error, "\n")

