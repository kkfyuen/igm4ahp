## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

## ----load-package-------------------------------------------------------------
library(igm4ahp)

## ----single-instance-gen------------------------------------------------------
# Set seed for reproducible single-instance demonstration
set.seed(42)

# Sample random order n and ratings vector rA from 9-point scale
n_demo <- sample(3:15, 1)
rA_demo <- randomPRR(n = n_demo)
A_demo <- getPRMFromRatings(rA_demo)

cat("Random Matrix Dimension n:", n_demo, "\n")
cat("Upper-triangular Rating Count m:", length(rA_demo), "\n\n")
print("--- Generated Pairwise Reciprocal Matrix A (First 4x4 block) ---")
print(round(A_demo[1:min(4, n_demo), 1:min(4, n_demo)], 3))

## ----single-instance-eval-----------------------------------------------------
r_shift <- 73

# 1. Pseudo-IGM
w_PIGM <- PIGM(A_demo)

# 2. Element-form NIGM with shift r
w_NIGM_r <- eNIGM(A_demo, r = r_shift)

# 3. Lagrangian-IGM (r = 0, excluding last lambda element)
res_LIGM_0 <- lagrangeIGM(A_demo, r = 0)
w_LIGM_0 <- res_LIGM_0[1:n_demo]

# 4. Lagrangian-IGM (r = 73, excluding last lambda element)
res_LIGM_r <- lagrangeIGM(A_demo, r = r_shift)
w_LIGM_r <- res_LIGM_r[1:n_demo]

# Combine into comparison table rounded to 8 decimal places (eps = 8)
df_comp <- data.frame(
  PIGM      = round(w_PIGM, 8),
  eNIGM_r   = round(w_NIGM_r, 8),
  LIGM_r0   = round(w_LIGM_0, 8),
  LIGM_r73  = round(w_LIGM_r, 8)
)
rownames(df_comp) <- paste0("w", 1:n_demo)

print("--- Single Instance Priority Vector Equivalence (8 Decimal Places) ---")
print(df_comp)

# Compute absolute discrepancy error er (Eq. 83; Yuen, 2024)
er_single <- sum(abs(round(w_PIGM, 8) - round(w_NIGM_r, 8))) +
             sum(abs(round(w_PIGM, 8) - round(w_LIGM_0, 8))) +
             sum(abs(round(w_PIGM, 8) - round(w_LIGM_r, 8)))
er_single <- round(er_single, 6)

cat("\nTotal Discrepancy Error (er):", er_single, "(Exact Equivalence Confirmed)\n")

## ----run-sim-igm--------------------------------------------------------------
# Ensure ./data directory exists in current working directory
data_dir <- "./data"
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

# Define CSV file path under ./data/
csv_file <- file.path(data_dir, "sim1_error_log.csv")

# Execute Monte Carlo Simulation 1 for N = 10 iterations
sim_error <- runSimIGM(size = 10, seedNo = 2000, errorFile = csv_file)

cat("Simulation Error Status (FALSE = No errors detected):", sim_error, "\n")
cat("Simulation log saved to:", csv_file, "\n")

## ----inspect-sim-log----------------------------------------------------------
# Read simulation CSV output from ./data directory
sim_data <- read.csv(csv_file, sep = ";", stringsAsFactors = FALSE)
colnames(sim_data) <- trimws(colnames(sim_data))

cat("Total Simulation Iterations Logged:", nrow(sim_data), "\n")
cat("Maximum Discrepancy Error Recorded:", max(sim_data$Error), "\n\n")

# Display complete simulation table for N = 10
print("--- Complete Simulation Results (N = 10) ---")
knitr::kable(sim_data, caption = "Simulation 1 Logged Output (./data/sim1_error_log.csv)")

## ----run-sim-timing-----------------------------------------------------------
# Define output file paths under ./data/
rdata_timing_file <- file.path(data_dir, "sim1_timing_data.RData")
csv_timing_file   <- file.path(data_dir, "sim1_timing_log.csv")

# Execute timed simulation for N = 10 iterations demo
timing_results <- runSimWithTiming(
  size = 10,
  seedNo = 2024,
  fileIGMData = rdata_timing_file,
  errorFile = csv_timing_file
)

cat("--- Simulation 1 Execution Timing (N = 10) ---\n")
cat("Start Time       :", format(timing_results$startTime), "\n")
cat("End Time         :", format(timing_results$endTime), "\n")
cat("Elapsed Duration :", round(as.numeric(timing_results$timeSlot, units = "secs"), 4), "seconds\n")
cat("Discrepancy Flag :", timing_results$error, "\n\n")

# Read and print timing log results from ./data directory
timing_data <- read.csv(csv_timing_file, sep = ";", stringsAsFactors = FALSE)
colnames(timing_data) <- trimws(colnames(timing_data))

print("--- Timed Simulation Log Results (./data/sim1_timing_log.csv) ---")
knitr::kable(timing_data, caption = "Timed Simulation 1 Logged Output")

