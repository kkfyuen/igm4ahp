## ----setup, include=FALSE---------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
options(width = 120)

## ----load-package-----------------------------------------------------------------------------------------------------
library(igm4ahp)

## ----example1-matrix--------------------------------------------------------------------------------------------------
# Construct 4x4 Pairwise Reciprocal Matrix A1
A1 <- matrix(c(
  1,   2,   4,   8,
  1/2, 1,   2,   4,
  1/4, 1/2, 1,   2,
  1/8, 1/4, 1/2, 1
), nrow = 4, byrow = TRUE)

rownames(A1) <- colnames(A1) <- paste0("C", 1:4)
print("--- Pairwise Reciprocal Matrix A1 ---")
print(A1)

# Verify Consistency Ratio (CR)
cr_A1 <- ConsistencyRatio(A1)
cat("Consistency Ratio (CR):", cr_A1, "(Perfectly Consistent)\n")

## ----example1-design--------------------------------------------------------------------------------------------------
D1 <- fullCoefficientMatrix(A1)
cat("Dimension of Design Matrix D1:", dim(D1)[1], "x", dim(D1)[2], "\n\n")
print("--- Full Design Matrix D1 (13x4) ---")
print(D1)

## ----example1-gram----------------------------------------------------------------------------------------------------
# Compute Gram Matrix G = D^T * D
G1 <- t(D1) %*% D1
print("--- Gram Matrix G1 = D1^T * D1 ---")
print(round(G1, 3))

# Compute Inverse Gram Matrix G^(-1)
G1_inv <- solve(G1)
print("--- Inverse Gram Matrix G1^(-1) ---")
print(round(G1_inv, 4))

## ----example1-pigm----------------------------------------------------------------------------------------------------
# Compute unnormalized row sum vector v = e * G^(-1)
v1 <- rowSums(G1_inv)
cat("Unnormalized Row Sum Vector v:", round(v1, 4), "\n")
cat("Sum of Row Sums (e * G^(-1) * e^T):", sum(v1), "\n\n")

# Closed-form priority vector via PIGM
w_PIGM_A1 <- PIGM(A1)
names(w_PIGM_A1) <- paste0("C", 1:4)

print("--- PIGM Priority Weight Vector (Section 4.1) ---")
print(round(w_PIGM_A1, 4))

# Analytical fractions
cat("\nExact Analytical Fractions:", paste(MASS::fractions(w_PIGM_A1), collapse = ", "), "\n")

## ----example1-blankmeyer----------------------------------------------------------------------------------------------
# Construct unshifted Gram matrix C1
C1 <- eleGram(A1, r = 0)
print("--- Blankmeyer Unshifted Gram Matrix C1 (r = 0) ---")
print(round(C1, 3))

# Check determinant of C1
det_C1 <- det(C1)
cat("\nDeterminant of Unshifted Matrix C1:", det_C1, "\n")
cat("Invertibility Status: Singular (C1^(-1) does NOT exist for consistent PRM)\n")

## ----example1-nigm----------------------------------------------------------------------------------------------------
# eNIGM with r = 5
res_nigm_r5 <- eNIGM(A1, r = 5, detail = TRUE)

print("--- Shifted Gram Matrix G (r = 5) ---")
print(round(res_nigm_r5$G, 3))

print("--- Shifted Inverse Gram Matrix G^(-1) (r = 5) ---")
print(round(res_nigm_r5$G_inv, 4))

cat("\nUnnormalized Row Sums v (r = 5):", round(res_nigm_r5$rs_G_inv, 4), "\n")
cat("Sum of Row Sums (r = 5):", sum(res_nigm_r5$rs_G_inv), "\n")

w_NIGM_r5 <- res_nigm_r5$W
cat("Normalized Priority Vector W (r = 5):", round(w_NIGM_r5, 4), "\n")

## ----example1-ligm----------------------------------------------------------------------------------------------------
# Lagrangian-IGM with r = 0 (using C1 in top-left block)
res_LIGM_r0 <- lagrangeIGM(A1, r = 0, detail = TRUE)

print("--- Augmented KKT Gram Matrix G_bar (5x5) ---")
print(round(res_LIGM_r0$G, 3))

print("--- KKT Inverse Gram Matrix G_bar^(-1) (5x5) ---")
print(round(res_LIGM_r0$G_inv, 4))

w_LIGM_r0 <- res_LIGM_r0$W[1:4]
lambda_r0 <- res_LIGM_r0$W[5]

cat("\nExtracted Priority Weights W:", round(w_LIGM_r0, 4), "\n")
cat("Extracted Lagrange Multiplier lambda (r = 0):", round(lambda_r0, 4), "\n")

## ----example1-comparison----------------------------------------------------------------------------------------------
w_EVM <- SaatyEigen(A1)

df_ex1_summary <- data.frame(
  PIGM = round(w_PIGM_A1, 4),
  NIGM_r1 = round(eNIGM(A1, r = 1), 4),
  NIGM_r5 = round(w_NIGM_r5, 4),
  LIGM_r0 = round(w_LIGM_r0, 4),
  Saaty_EVM = round(w_EVM, 4)
)
rownames(df_ex1_summary) <- paste0("Criterion ", 1:4)

print("--- Example 1 Method Comparison ---")
print(df_ex1_summary)

## ----example2-matrix--------------------------------------------------------------------------------------------------
# Construct 6x6 inconsistent PRM A2
A2 <- matrix(c(
  1,   4,   3, 1,   3,   4,
  1/4, 1,   7, 3,   1/5, 1,
  1/3, 1/7, 1, 1/5, 1/5, 1/6,
  1,   1/3, 5, 1,   1,   1/3,
  1/3, 5,   5, 1,   1,   3,
  1/4, 1,   6, 3,   1/3, 1
), nrow = 6, byrow = TRUE)

criteria_names <- c("Learning", "Friends", "School Life", "Vocational", "College Prep", "Music Classes")
rownames(A2) <- colnames(A2) <- criteria_names

print("--- Pairwise Reciprocal Matrix A2 ---")
print(round(A2, 3))

# Consistency Evaluation
ci_A2 <- ConsistencyIndex(A2)
cr_A2 <- ConsistencyRatio(A2)

cat("\nConsistency Index (CI):", round(ci_A2, 3), "\n")
cat("Consistency Ratio (CR):", round(cr_A2, 3), "(Inconsistent since CR > 0.10)\n")

## ----example2-design--------------------------------------------------------------------------------------------------
D2 <- fullCoefficientMatrix(A2)
cat("Dimension of Design Matrix D2:", dim(D2)[1], "x", dim(D2)[2], "\n")
cat("First 4 differential rows and bottom row e of D2:\n")
print(round(rbind(head(D2, 4), tail(D2, 1)), 3))

## ----example2-gram----------------------------------------------------------------------------------------------------
# Gram matrix G2 = D2^T * D2
G2 <- t(D2) %*% D2
print("--- Gram Matrix G2 (6x6) ---")
print(round(G2, 3))

# Inverse Gram Matrix G2^(-1)
G2_inv <- solve(G2)
print("--- Inverse Gram Matrix G2^(-1) (6x6) ---")
print(round(G2_inv, 3))

## ----example2-pigm----------------------------------------------------------------------------------------------------
# Unnormalized row sums v = e * G2^(-1)
v2 <- rowSums(G2_inv)
cat("Unnormalized Row Sum Vector v2:", round(v2, 3), "\n")
cat("Sum of Row Sums (e * G2^(-1) * e^T):", round(sum(v2), 3), "\n\n")

# Normalized PIGM priority vector
w_PIGM_A2 <- PIGM(A2)
names(w_PIGM_A2) <- criteria_names

print("--- PIGM Priority Weights & Rankings ---")
df_pigm_rank <- data.frame(
  Weight = round(w_PIGM_A2, 3),
  Rank   = rank(-w_PIGM_A2)
)
print(df_pigm_rank)

## ----example2-nigm----------------------------------------------------------------------------------------------------
# Unshifted NIGM (r = 0)
res_nigm_r0 <- eNIGM(A2, r = 0, detail = TRUE)

print("--- Unshifted Inverse Matrix C2^(-1) (r = 0) ---")
print(round(res_nigm_r0$G_inv, 3))

cat("\nUnnormalized Row Sums v_C (r = 0):", round(res_nigm_r0$rs_G_inv, 3), "\n")
cat("Sum of Row Sums (r = 0):", round(sum(res_nigm_r0$rs_G_inv), 3), "\n")

w_NIGM_r0 <- res_nigm_r0$W
cat("Normalized Priority Vector W (r = 0):", round(w_NIGM_r0, 3), "\n")

# Shifted NIGM (r = 1)
w_NIGM_r1 <- eNIGM(A2, r = 1)
cat("Normalized Priority Vector W (r = 1):", round(w_NIGM_r1, 3), "\n")

## ----example2-ligm----------------------------------------------------------------------------------------------------
# LIGM with r = 0
res_LIGM_r0_A2 <- lagrangeIGM(A2, r = 0, detail = TRUE)
w_LIGM_r0_A2 <- res_LIGM_r0_A2$W[1:6]
lambda_r0_A2 <- res_LIGM_r0_A2$G_inv[7, 7]

# LIGM with r = 1 (matching manuscript Passage 98)
res_LIGM_r1_A2 <- lagrangeIGM(A2, r = 1, detail = TRUE)

print("--- Augmented KKT Inverse Gram Matrix G_bar^(-1) (7x7, r = 1) ---")
print(round(res_LIGM_r1_A2$G_inv, 3))

w_LIGM_r1_A2 <- res_LIGM_r1_A2$W[1:6]
lambda_r1_A2 <- res_LIGM_r1_A2$G_inv[7, 7]

cat("\nExtracted Priority Weights W (r = 1):", round(w_LIGM_r1_A2, 3), "\n")
cat("Extracted Lagrange Multiplier lambda (r = 0):", round(lambda_r0_A2, 3), "\n")
cat("Extracted Lagrange Multiplier lambda (r = 1):", round(lambda_r1_A2, 3), "(Matches manuscript Passage 98)\n")

## ----example2-wlse----------------------------------------------------------------------------------------------------
# WLSE for IGM closed-form weights
wlse_PIGM_A2 <- WLSE(w_PIGM_A2, A2)

# WLSE for Saaty EVM weights
w_EVM_A2 <- SaatyEigen(A2)
wlse_EVM_A2 <- WLSE(w_EVM_A2, A2)

cat("WLSE (PIGM / IGM Closed-Form) :", round(wlse_PIGM_A2, 4), "\n")
cat("WLSE (Saaty Eigenvalue Method):", round(wlse_EVM_A2, 4), "\n")
cat("WLS Error Advantage of IGM   :", round(wlse_EVM_A2 - wlse_PIGM_A2, 4), "\n")

## ----summary-table----------------------------------------------------------------------------------------------------
summary_table <- data.frame(
  Example = c("Example 1 (4x4 Consistent)", "Example 2 (6x6 Inconsistent)"),
  Dimension = c("4 x 4", "6 x 6"),
  Consistency_CR = c(round(cr_A1, 3), round(cr_A2, 3)),
  PIGM_Top_Rank = c("C1 (0.533)", "Learning (0.415)"),
  WLSE_Error = c(0.0000, round(wlse_PIGM_A2, 4))
)

knitr::kable(summary_table, caption = "Summary of Section 4 Numerical Demonstrations")

