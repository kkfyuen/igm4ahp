# Package index

## Core Inverse Gram Matrix (IGM) Methods

Closed-form Weighted Least Squares optimization solvers for AHP PRMs.

- [`PIGM()`](https://kkfyuen.github.io/igm4ahp/reference/PIGM.md) :
  Pseudo Inverse Gram Matrix (PIGM) Prioritization
- [`pseudoIGM()`](https://kkfyuen.github.io/igm4ahp/reference/pseudoIGM.md)
  : Pseudo Inverse Gram Matrix (PIGM) Prioritization
- [`eNIGM()`](https://kkfyuen.github.io/igm4ahp/reference/eNIGM.md) :
  Element-form Normalized Inverse Gram Matrix (eNIGM) Prioritization
- [`lagrangeIGM()`](https://kkfyuen.github.io/igm4ahp/reference/lagrangeIGM.md)
  [`LIGM()`](https://kkfyuen.github.io/igm4ahp/reference/lagrangeIGM.md)
  : Lagrangian Inverse Gram Matrix (LIGM) Prioritization
- [`NIGM()`](https://kkfyuen.github.io/igm4ahp/reference/NIGM.md) :
  Matrix-form Normalized Inverse Gram Matrix (NIGM) Prioritization
- [`PNIGM()`](https://kkfyuen.github.io/igm4ahp/reference/PNIGM.md) :
  Pseudo-Normalized Inverse Gram Matrix (PNIGM) Prioritization
- [`eleGram()`](https://kkfyuen.github.io/igm4ahp/reference/eleGram.md)
  : Element-wise Gram Matrix Construction

## Saaty EVM & Evaluation Utilities

Saaty’s Principal Eigenvector method and consistency ratios.

- [`SaatyEigen()`](https://kkfyuen.github.io/igm4ahp/reference/SaatyEigen.md)
  : Principal Eigenvector Priority Derivation (Saaty's EVM)
- [`ConsistencyIndex()`](https://kkfyuen.github.io/igm4ahp/reference/ConsistencyIndex.md)
  : Calculate Consistency Index (CI) for Pairwise Reciprocal Matrix
- [`ConsistencyRatio()`](https://kkfyuen.github.io/igm4ahp/reference/ConsistencyRatio.md)
  : Calculate Consistency Ratio (CR) for Pairwise Reciprocal Matrix
- [`WLSE()`](https://kkfyuen.github.io/igm4ahp/reference/WLSE.md) :
  Weighted Least Squares Error (WLSE) Objective Function

## Numerical Solvers & Optimization Benchmarks

Non-linear WLS optimization solvers for benchmarking.

- [`WLSopt()`](https://kkfyuen.github.io/igm4ahp/reference/WLSopt.md) :
  Numerical Weighted Least Squares (WLS) Optimization Solver for AHP

## Matrix & Simulation Helpers

Helper functions for design coefficient matrices, matrix reconstruction,
and simulation benchmarking.

- [`fullCoefficientMatrix()`](https://kkfyuen.github.io/igm4ahp/reference/fullCoefficientMatrix.md)
  : Construct Full Design (Coefficient) Matrix for Pairwise Comparisons
- [`getPRMFromRatings()`](https://kkfyuen.github.io/igm4ahp/reference/getPRMFromRatings.md)
  : Construct Pairwise Reciprocal Matrix from Ratings Vector
- [`rDesignMatrix()`](https://kkfyuen.github.io/igm4ahp/reference/rDesignMatrix.md)
  : Construct Parametrized Design Matrix for NIGM
- [`randomPRR()`](https://kkfyuen.github.io/igm4ahp/reference/randomPRR.md)
  : Generate Random Pairwise Comparison Ratings Vector
- [`runSimIGM()`](https://kkfyuen.github.io/igm4ahp/reference/runSimIGM.md)
  : Simulation Verification for Inverse Gram Matrix (IGM) Methods
- [`runSimWithTiming()`](https://kkfyuen.github.io/igm4ahp/reference/runSimWithTiming.md)
  : Run and Benchmark Inverse Gram Matrix Simulation
- [`runWLSEOpt()`](https://kkfyuen.github.io/igm4ahp/reference/runWLSEOpt.md)
  : Simulation Verification Comparing Pseudo-IGM and Non-Linear WLS
  Optimization
- [`runWLSEOptWithTiming()`](https://kkfyuen.github.io/igm4ahp/reference/runWLSEOptWithTiming.md)
  : Run and Benchmark Non-Linear WLS Optimization Simulation

## Validation Utilities

Validation functions for matrix dimensions and input checks.

- [`isWholeNo()`](https://kkfyuen.github.io/igm4ahp/reference/isWholeNo.md)
  : Check if a Value is a Whole Number
- [`is.SizeOfUTM()`](https://kkfyuen.github.io/igm4ahp/reference/is.SizeOfUTM.md)
  : Check if Input Object is an Upper-Triangular Matrix (UTM) Ratings
  Vector
