library(testthat)
library(igm4ahp)

test_that("PIGM computes exact priority weights for consistent 4x4 PRM (Example 1)", {
  # Construct consistent 4x4 matrix A1
  A1 <- matrix(c(
    1,   2,   4,   8,
    1/2, 1,   2,   4,
    1/4, 1/2, 1,   2,
    1/8, 1/4, 1/2, 1
  ), nrow = 4, byrow = TRUE)
  
  expected_w <- c(8/15, 4/15, 2/15, 1/15)
  
  # Check PIGM calculation
  w_pigm <- PIGM(A1)
  expect_equal(w_pigm, expected_w, tolerance = 1e-6)
  
  # Check equivalence with SaatyEigen and eNIGM
  expect_equal(SaatyEigen(A1), expected_w, tolerance = 1e-6)
  expect_equal(eNIGM(A1, r = 1), expected_w, tolerance = 1e-6)
  
  # Check consistency ratio
  expect_equal(ConsistencyRatio(A1), 0, tolerance = 1e-6)
})

test_that("PIGM, eNIGM, and lagrangeIGM compute valid weights, lambda, and WLSE for inconsistent 6x6 PRM (Example 2)", {
  # Construct 6x6 inconsistent PRM A2
  A2 <- matrix(c(
    1,   4,   3, 1,   3,   4,
    1/4, 1,   7, 3,   1/5, 1,
    1/3, 1/7, 1, 1/5, 1/5, 1/6,
    1,   1/3, 5, 1,   1,   1/3,
    1/3, 5,   5, 1,   1,   3,
    1/4, 1,   6, 3,   1/3, 1
  ), nrow = 6, byrow = TRUE)
  
  w_pigm <- PIGM(A2)
  
  # Check unit sum normalization
  expect_equal(sum(w_pigm), 1.0, tolerance = 1e-6)
  
  # Check equivalence with eNIGM
  w_nigm <- eNIGM(A2, r = 1)
  expect_equal(w_pigm, w_nigm, tolerance = 1e-5)
  
  # Check consistency ratio > 0.10
  cr <- ConsistencyRatio(A2)
  expect_gt(cr, 0.10)
  
  # Check Lagrangian-IGM (lagrangeIGM) weights and output
  res_ligm <- lagrangeIGM(A2, r = 0)
  w_ligm <- res_ligm[1:6]
  lambda_val <- res_ligm[7]
  
  expect_equal(w_pigm, w_ligm, tolerance = 1e-5)
  expect_equal(lambda_val, 0.3665, tolerance = 1e-3)
  
  # Check Weighted Least Squares Error (WLSE)
  wlse_pigm <- WLSE(w_pigm, A2)
  expect_equal(wlse_pigm, 0.6335, tolerance = 1e-3)
})

test_that("is.SizeOfUTM correctly validates upper-triangular rating vector lengths", {
  expect_true(is.SizeOfUTM(c(2, 5, 3)))       # n = 3 -> 3*2/2 = 3
  expect_false(is.SizeOfUTM(c(1, 2, 3, 4)))   # 4 is not a valid n(n-1)/2
})
