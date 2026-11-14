# Section 4: Numerical Demonstration of Inverse Gram Matrix (IGM) Methods

## 1. Overview & Theoretical Framework

This vignette provides a comprehensive, step-by-step numerical
reproduction of **Section 4 (“Numerical Examples”)** from *Yuen (2024)*
(*“Inverse Gram Matrix Methods for Prioritization in Analytic Hierarchy
Process: Explainability of Weighted Least Squares Optimization
Method”*).

The primary objective is to demonstrate the practical application of the
closed-form **Inverse Gram Matrix (IGM)** prioritization
framework—specifically **Pseudo-IGM (`PIGM`)**, **Normalized-IGM
(`eNIGM` / `NIGM`)**, and **Lagrangian-IGM (`lagrangeIGM`)**—on two
benchmark decision scenarios:

1.  **Section 4.1 (Example 1)**: A perfectly consistent \\4 \times 4\\
    Pairwise Reciprocal Matrix (PRM) \\A_1\\, demonstrating that `PIGM`
    resolves the matrix singularity issue inherent in earlier unshifted
    formulation attempts (Blankmeyer, 1987) while matching Saaty’s
    Eigenvalue Method (`SaatyEigen`).
2.  **Section 4.2 (Example 2)**: An inconsistent \\6 \times 6\\ PRM
    \\A_2\\ evaluating high school selection criteria (Saaty, 1980),
    demonstrating exact closed-form Weighted Least Squares (WLS)
    prioritization without relying on iterative numerical optimization
    algorithms (`WLSopt`).

------------------------------------------------------------------------

## 2. Package Setup

Load the **`igm4ahp`** package containing all closed-form IGM
algorithms, Saaty Eigenvalue utilities, and WLS objective error
functions.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`igm4ahp`](https://kkfyuen.github.io/igm4ahp)`)`

------------------------------------------------------------------------

## 3. Section 4.1: Example 1 — Perfectly Consistent 4x4 PRM

### 3.1. Problem Setup & Matrix Formulation

Consider a perfectly consistent \\4 \times 4\\ Pairwise Reciprocal
Matrix \\A_1\\ (Yuen, 2024, Section 4.1):

\\A_1 = \begin{pmatrix} 1 & 2 & 4 & 8 \\ 1/2 & 1 & 2 & 4 \\ 1/4 & 1/2 &
1 & 2 \\ 1/8 & 1/4 & 1/2 & 1 \end{pmatrix}\\

\
`# Construct 4x4 Pairwise Reciprocal Matrix A1`\
`A1`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(`\
`  ``1``,   ``2``,   ``4``,   ``8``,`\
`  ``1``/``2``, ``1``,   ``2``,   ``4``,`\
`  ``1``/``4``, ``1``/``2``, ``1``,   ``2``,`\
`  ``1``/``8``, ``1``/``4``, ``1``/``2``, ``1`\
`)``, nrow ``=`` ``4``, byrow ``=`` ``TRUE``)`\
\
[`rownames`](https://rdrr.io/r/base/colnames.html)`(``A1``)`` ``<-`` `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``A1``)`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"C"``, ``1``:``4``)`\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Pairwise Reciprocal Matrix A1 ---"``)`

    ## [1] "--- Pairwise Reciprocal Matrix A1 ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(``A1``)`

    ##       C1   C2  C3 C4
    ## C1 1.000 2.00 4.0  8
    ## C2 0.500 1.00 2.0  4
    ## C3 0.250 0.50 1.0  2
    ## C4 0.125 0.25 0.5  1

\
`# Verify Consistency Ratio (CR)`\
`cr_A1`` ``<-`` `[`ConsistencyRatio`](https://kkfyuen.github.io/igm4ahp/reference/ConsistencyRatio.md)`(``A1``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Consistency Ratio (CR):"``, ``cr_A1``, ``"(Perfectly Consistent)\n"``)`

    ## Consistency Ratio (CR): 0 (Perfectly Consistent)

### 3.2. Step 1: Design Coefficient Matrix D Construction

Using `fullCoefficientMatrix(A1)` (Algorithm 1; Yuen, 2024), we
construct the \\(4(4-1)+1) \times 4 = 13 \times 4\\ design coefficient
matrix \\D = \begin{bmatrix} D_e \\ e \end{bmatrix}\\. The first 12 rows
(\\D_e\\) represent the pairwise weighted differential equations
(\\w_i - a\_{ij} w_j = 0\\), and the 13th row \\e = \[1, 1, 1, 1\]\\
represents the unit simplex normalization constraint (\\\sum\_{i=1}^4
w_i = 1\\).

\
`D1`` ``<-`` `[`fullCoefficientMatrix`](https://kkfyuen.github.io/igm4ahp/reference/fullCoefficientMatrix.md)`(``A1``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Dimension of Design Matrix D1:"``, `[`dim`](https://rdrr.io/r/base/dim.html)`(``D1``)``[``1``]``, ``"x"``, `[`dim`](https://rdrr.io/r/base/dim.html)`(``D1``)``[``2``]``, ``"\n\n"``)`

    ## Dimension of Design Matrix D1: 13 x 4

\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Full Design Matrix D1 (13x4) ---"``)`

    ## [1] "--- Full Design Matrix D1 (13x4) ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(``D1``)`

    ##           w1    w2   w3 w4
    ##  [1,]  1.000 -2.00  0.0  0
    ##  [2,] -0.500  1.00  0.0  0
    ##  [3,]  1.000  0.00 -4.0  0
    ##  [4,] -0.250  0.00  1.0  0
    ##  [5,]  1.000  0.00  0.0 -8
    ##  [6,] -0.125  0.00  0.0  1
    ##  [7,]  0.000  1.00 -2.0  0
    ##  [8,]  0.000 -0.50  1.0  0
    ##  [9,]  0.000  1.00  0.0 -4
    ## [10,]  0.000 -0.25  0.0  1
    ## [11,]  0.000  0.00  1.0 -2
    ## [12,]  0.000  0.00 -0.5  1
    ## [13,]  1.000  1.00  1.0  1

### 3.3. Step 2 & 3: Gram Matrix G and Inverse Gram Matrix G^(-1)

The \\4 \times 4\\ Gram matrix \\G_1 = D_1^T D_1\\ is constructed via
matrix multiplication or directly using the element-wise formula
(Theorem 2; Yuen, 2024): - **Diagonal entries**: \\g\_{ii} = (n-1) +
\sum\_{k \neq i} a\_{ki}^2 + 1\\ - **Off-diagonal entries**: \\g\_{ij} =
1 - a\_{ij} - a\_{ji}\\

\
`# Compute Gram Matrix G = D^T * D`\
`G1`` ``<-`` `[`t`](https://rdrr.io/r/base/t.html)`(``D1``)`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``D1`\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Gram Matrix G1 = D1^T * D1 ---"``)`

    ## [1] "--- Gram Matrix G1 = D1^T * D1 ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``G1``, ``3``)``)`

    ##        w1     w2    w3     w4
    ## w1  4.328 -1.500 -3.25 -7.125
    ## w2 -1.500  8.312 -1.50 -3.250
    ## w3 -3.250 -1.500 24.25 -1.500
    ## w4 -7.125 -3.250 -1.50 88.000

\
`# Compute Inverse Gram Matrix G^(-1)`\
`G1_inv`` ``<-`` `[`solve`](https://rdrr.io/r/base/solve.html)`(``G1``)`\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Inverse Gram Matrix G1^(-1) ---"``)`

    ## [1] "--- Inverse Gram Matrix G1^(-1) ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``G1_inv``, ``4``)``)`

    ##        w1     w2     w3     w4
    ## w1 0.3574 0.0874 0.0554 0.0331
    ## w2 0.0874 0.1450 0.0215 0.0128
    ## w3 0.0554 0.0215 0.0504 0.0061
    ## w4 0.0331 0.0128 0.0061 0.0146

### 3.4. Step 4: Unnormalized Row Sums & Closed-Form PIGM Priority Vector

By Theorem 1 (Yuen, 2024), the unnormalized weight vector \\v = e
G^{-1}\\ is obtained by summing the rows of \\G_1^{-1}\\:

\\W = \frac{e G^{-1}}{e G^{-1} e^T} = \frac{e (D^T D)^{-1}}{e (D^T
D)^{-1} e^T}\\

\
`# Compute unnormalized row sum vector v = e * G^(-1)`\
`v1`` ``<-`` `[`rowSums`](https://rdrr.io/r/base/colSums.html)`(``G1_inv``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Unnormalized Row Sum Vector v:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``v1``, ``4``)``, ``"\n"``)`

    ## Unnormalized Row Sum Vector v: 0.5333 0.2667 0.1333 0.0667

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Sum of Row Sums (e * G^(-1) * e^T):"``, `[`sum`](https://rdrr.io/r/base/sum.html)`(``v1``)``, ``"\n\n"``)`

    ## Sum of Row Sums (e * G^(-1) * e^T): 1

\
`# Closed-form priority vector via PIGM`\
`w_PIGM_A1`` ``<-`` `[`PIGM`](https://kkfyuen.github.io/igm4ahp/reference/PIGM.md)`(``A1``)`\
[`names`](https://rdrr.io/r/base/names.html)`(``w_PIGM_A1``)`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"C"``, ``1``:``4``)`\
\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- PIGM Priority Weight Vector (Section 4.1) ---"``)`

    ## [1] "--- PIGM Priority Weight Vector (Section 4.1) ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``w_PIGM_A1``, ``4``)``)`

    ##     C1     C2     C3     C4 
    ## 0.5333 0.2667 0.1333 0.0667

\
`# Analytical fractions`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"\nExact Analytical Fractions:"``, `[`paste`](https://rdrr.io/r/base/paste.html)`(``MASS``::`[`fractions`](https://rdrr.io/pkg/MASS/man/fractions.html)`(``w_PIGM_A1``)``, collapse ``=`` ``", "``)``, ``"\n"``)`

    ## 
    ## Exact Analytical Fractions: 8/15, 4/15, 2/15, 1/15

### 3.5. Blankmeyer (1987) Unshifted Gram Matrix C Singularity Analysis

For a perfectly consistent matrix, earlier formulations (Blankmeyer,
1987) failed because the unshifted Gram matrix \\C_1 = G_1 - E = D_e^T
D_e\\ was singular (\\\det C_1 = 0\\):

\
`# Construct unshifted Gram matrix C1`\
`C1`` ``<-`` `[`eleGram`](https://kkfyuen.github.io/igm4ahp/reference/eleGram.md)`(``A1``, r ``=`` ``0``)`\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Blankmeyer Unshifted Gram Matrix C1 (r = 0) ---"``)`

    ## [1] "--- Blankmeyer Unshifted Gram Matrix C1 (r = 0) ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``C1``, ``3``)``)`

    ##        [,1]   [,2]  [,3]   [,4]
    ## [1,]  3.328 -2.500 -4.25 -8.125
    ## [2,] -2.500  7.312 -2.50 -4.250
    ## [3,] -4.250 -2.500 23.25 -2.500
    ## [4,] -8.125 -4.250 -2.50 87.000

\
`# Check determinant of C1`\
`det_C1`` ``<-`` `[`det`](https://rdrr.io/r/base/det.html)`(``C1``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"\nDeterminant of Unshifted Matrix C1:"``, ``det_C1``, ``"\n"``)`

    ## 
    ## Determinant of Unshifted Matrix C1: 6.096824e-12

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Invertibility Status: Singular (C1^(-1) does NOT exist for consistent PRM)\n"``)`

    ## Invertibility Status: Singular (C1^(-1) does NOT exist for consistent PRM)

### 3.6. Normalized Inverse Gram Matrix (eNIGM) with Parameter Shifts

`eNIGM` resolves the singularity of consistent PRMs by adding a
regularization scalar shift \\r \neq 0\\ (Theorem 3; Yuen, 2024). We
demonstrate calculations for \\r = 1\\ and \\r = 5\\:

\\\mathbf{G}\_r = C_1 + r \mathbf{E}, \quad W = \frac{e
\mathbf{G}\_r^{-1}}{e \mathbf{G}\_r^{-1} e^T}\\

\
`# eNIGM with r = 5`\
`res_nigm_r5`` ``<-`` `[`eNIGM`](https://kkfyuen.github.io/igm4ahp/reference/eNIGM.md)`(``A1``, r ``=`` ``5``, detail ``=`` ``TRUE``)`\
\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Shifted Gram Matrix G (r = 5) ---"``)`

    ## [1] "--- Shifted Gram Matrix G (r = 5) ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``res_nigm_r5``$``G``, ``3``)``)`

    ##        [,1]   [,2]  [,3]   [,4]
    ## [1,]  8.328  2.500  0.75 -3.125
    ## [2,]  2.500 12.312  2.50  0.750
    ## [3,]  0.750  2.500 28.25  2.500
    ## [4,] -3.125  0.750  2.50 92.000

\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Shifted Inverse Gram Matrix G^(-1) (r = 5) ---"``)`

    ## [1] "--- Shifted Inverse Gram Matrix G^(-1) (r = 5) ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``res_nigm_r5``$``G_inv``, ``4``)``)`

    ##         [,1]    [,2]    [,3]    [,4]
    ## [1,]  0.1299 -0.0263 -0.0015  0.0047
    ## [2,] -0.0263  0.0881 -0.0070 -0.0014
    ## [3,] -0.0015 -0.0070  0.0361 -0.0010
    ## [4,]  0.0047 -0.0014 -0.0010  0.0111

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"\nUnnormalized Row Sums v (r = 5):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``res_nigm_r5``$``rs_G_inv``, ``4``)``, ``"\n"``)`

    ## 
    ## Unnormalized Row Sums v (r = 5): 0.1067 0.0533 0.0267 0.0133

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Sum of Row Sums (r = 5):"``, `[`sum`](https://rdrr.io/r/base/sum.html)`(``res_nigm_r5``$``rs_G_inv``)``, ``"\n"``)`

    ## Sum of Row Sums (r = 5): 0.2

\
`w_NIGM_r5`` ``<-`` ``res_nigm_r5``$``W`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Normalized Priority Vector W (r = 5):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``w_NIGM_r5``, ``4``)``, ``"\n"``)`

    ## Normalized Priority Vector W (r = 5): 0.5333 0.2667 0.1333 0.0667

The parameter \\r\\ cancels out completely during post-hoc vector
normalization (\\0.1067 / 0.2000 = 0.5333\\).

### 3.7. Augmented KKT System via Lagrangian-IGM (lagrangeIGM)

By Theorem 4 (Yuen, 2024), `lagrangeIGM` solves the augmented KKT block
system \\\bar{G} \begin{bmatrix} W \\ \lambda \end{bmatrix} = q^T\\:

\\\bar{G} = \begin{pmatrix} C & e^T \\ e & 0 \end{pmatrix}, \quad
\begin{bmatrix} W \\ \lambda \end{bmatrix} = q \bar{G}^{-1}\\

\
`# Lagrangian-IGM with r = 0 (using C1 in top-left block)`\
`res_LIGM_r0`` ``<-`` `[`lagrangeIGM`](https://kkfyuen.github.io/igm4ahp/reference/lagrangeIGM.md)`(``A1``, r ``=`` ``0``, detail ``=`` ``TRUE``)`\
\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Augmented KKT Gram Matrix G_bar (5x5) ---"``)`

    ## [1] "--- Augmented KKT Gram Matrix G_bar (5x5) ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``res_LIGM_r0``$``G``, ``3``)``)`

    ##        [,1]   [,2]  [,3]   [,4] [,5]
    ## [1,]  3.328 -2.500 -4.25 -8.125    1
    ## [2,] -2.500  7.312 -2.50 -4.250    1
    ## [3,] -4.250 -2.500 23.25 -2.500    1
    ## [4,] -8.125 -4.250 -2.50 87.000    1
    ## [5,]  1.000  1.000  1.00  1.000    0

\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- KKT Inverse Gram Matrix G_bar^(-1) (5x5) ---"``)`

    ## [1] "--- KKT Inverse Gram Matrix G_bar^(-1) (5x5) ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``res_LIGM_r0``$``G_inv``, ``4``)``)`

    ##         [,1]    [,2]    [,3]    [,4]   [,5]
    ## [1,]  0.0730 -0.0548 -0.0158 -0.0024 0.5333
    ## [2,] -0.0548  0.0738 -0.0141 -0.0050 0.2667
    ## [3,] -0.0158 -0.0141  0.0326 -0.0028 0.1333
    ## [4,] -0.0024 -0.0050 -0.0028  0.0102 0.0667
    ## [5,]  0.5333  0.2667  0.1333  0.0667 0.0000

\
`w_LIGM_r0`` ``<-`` ``res_LIGM_r0``$``W``[``1``:``4``]`\
`lambda_r0`` ``<-`` ``res_LIGM_r0``$``W``[``5``]`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"\nExtracted Priority Weights W:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``w_LIGM_r0``, ``4``)``, ``"\n"``)`

    ## 
    ## Extracted Priority Weights W: 0.5333 0.2667 0.1333 0.0667

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Extracted Lagrange Multiplier lambda (r = 0):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``lambda_r0``, ``4``)``, ``"\n"``)`

    ## Extracted Lagrange Multiplier lambda (r = 0): 1

### 3.8. Example 1 Method Comparison

\
`w_EVM`` ``<-`` `[`SaatyEigen`](https://kkfyuen.github.io/igm4ahp/reference/SaatyEigen.md)`(``A1``)`\
\
`df_ex1_summary`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`  PIGM ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(``w_PIGM_A1``, ``4``)``,`\
`  NIGM_r1 ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(`[`eNIGM`](https://kkfyuen.github.io/igm4ahp/reference/eNIGM.md)`(``A1``, r ``=`` ``1``)``, ``4``)``,`\
`  NIGM_r5 ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(``w_NIGM_r5``, ``4``)``,`\
`  LIGM_r0 ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(``w_LIGM_r0``, ``4``)``,`\
`  Saaty_EVM ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(``w_EVM``, ``4``)`\
`)`\
[`rownames`](https://rdrr.io/r/base/colnames.html)`(``df_ex1_summary``)`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"Criterion "``, ``1``:``4``)`\
\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Example 1 Method Comparison ---"``)`

    ## [1] "--- Example 1 Method Comparison ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(``df_ex1_summary``)`

    ##               PIGM NIGM_r1 NIGM_r5 LIGM_r0 Saaty_EVM
    ## Criterion 1 0.5333  0.5333  0.5333  0.5333    0.5333
    ## Criterion 2 0.2667  0.2667  0.2667  0.2667    0.2667
    ## Criterion 3 0.1333  0.1333  0.1333  0.1333    0.1333
    ## Criterion 4 0.0667  0.0667  0.0667  0.0667    0.0667

All methods yield **exact numerical equivalence** (\\W = \[0.5333,
0.2667, 0.1333, 0.0667\]\\).

------------------------------------------------------------------------

## 4. Section 4.2: Example 2 — Inconsistent 6x6 PRM

### 4.1. Problem Setup: High School Selection Criteria

Section 4.2 evaluates an inconsistent \\6 \times 6\\ PRM \\A_2\\
representing six high school selection criteria (Saaty, 1980; Yuen,
2024): 1. **Learning** (\\C_1\\) 2. **Friends** (\\C_2\\) 3. **School
Life** (\\C_3\\) 4. **Vocational Training** (\\C_4\\) 5. **College
Preparation** (\\C_5\\) 6. **Music Classes** (\\C_6\\)

\\A_2 = \begin{pmatrix} 1 & 4 & 3 & 1 & 3 & 4 \\ 1/4 & 1 & 7 & 3 & 1/5 &
1 \\ 1/3 & 1/7 & 1 & 1/5 & 1/5 & 1/6 \\ 1 & 1/3 & 5 & 1 & 1 & 1/3 \\ 1/3
& 5 & 5 & 1 & 1 & 3 \\ 1/4 & 1 & 6 & 3 & 1/3 & 1 \end{pmatrix}\\

\
`# Construct 6x6 inconsistent PRM A2`\
`A2`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(`\
`  ``1``,   ``4``,   ``3``, ``1``,   ``3``,   ``4``,`\
`  ``1``/``4``, ``1``,   ``7``, ``3``,   ``1``/``5``, ``1``,`\
`  ``1``/``3``, ``1``/``7``, ``1``, ``1``/``5``, ``1``/``5``, ``1``/``6``,`\
`  ``1``,   ``1``/``3``, ``5``, ``1``,   ``1``,   ``1``/``3``,`\
`  ``1``/``3``, ``5``,   ``5``, ``1``,   ``1``,   ``3``,`\
`  ``1``/``4``, ``1``,   ``6``, ``3``,   ``1``/``3``, ``1`\
`)``, nrow ``=`` ``6``, byrow ``=`` ``TRUE``)`\
\
`criteria_names`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Learning"``, ``"Friends"``, ``"School Life"``, ``"Vocational"``, ``"College Prep"``, ``"Music Classes"``)`\
[`rownames`](https://rdrr.io/r/base/colnames.html)`(``A2``)`` ``<-`` `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``A2``)`` ``<-`` ``criteria_names`\
\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Pairwise Reciprocal Matrix A2 ---"``)`

    ## [1] "--- Pairwise Reciprocal Matrix A2 ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``A2``, ``3``)``)`

    ##               Learning Friends School Life Vocational College Prep Music Classes
    ## Learning         1.000   4.000           3        1.0        3.000         4.000
    ## Friends          0.250   1.000           7        3.0        0.200         1.000
    ## School Life      0.333   0.143           1        0.2        0.200         0.167
    ## Vocational       1.000   0.333           5        1.0        1.000         0.333
    ## College Prep     0.333   5.000           5        1.0        1.000         3.000
    ## Music Classes    0.250   1.000           6        3.0        0.333         1.000

\
`# Consistency Evaluation`\
`ci_A2`` ``<-`` `[`ConsistencyIndex`](https://kkfyuen.github.io/igm4ahp/reference/ConsistencyIndex.md)`(``A2``)`\
`cr_A2`` ``<-`` `[`ConsistencyRatio`](https://kkfyuen.github.io/igm4ahp/reference/ConsistencyRatio.md)`(``A2``)`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"\nConsistency Index (CI):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``ci_A2``, ``3``)``, ``"\n"``)`

    ## 
    ## Consistency Index (CI): 0.284

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Consistency Ratio (CR):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``cr_A2``, ``3``)``, ``"(Inconsistent since CR > 0.10)\n"``)`

    ## Consistency Ratio (CR): 0.229 (Inconsistent since CR > 0.10)

### 4.2. Step 1: Design Matrix D Construction (31x6)

Using `fullCoefficientMatrix(A2)`, we construct the \\(6(6-1)+1) \times
6 = 31 \times 6\\ design coefficient matrix \\D_2\\:

\
`D2`` ``<-`` `[`fullCoefficientMatrix`](https://kkfyuen.github.io/igm4ahp/reference/fullCoefficientMatrix.md)`(``A2``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Dimension of Design Matrix D2:"``, `[`dim`](https://rdrr.io/r/base/dim.html)`(``D2``)``[``1``]``, ``"x"``, `[`dim`](https://rdrr.io/r/base/dim.html)`(``D2``)``[``2``]``, ``"\n"``)`

    ## Dimension of Design Matrix D2: 31 x 6

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"First 4 differential rows and bottom row e of D2:\n"``)`

    ## First 4 differential rows and bottom row e of D2:

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(`[`rbind`](https://rdrr.io/r/base/cbind.html)`(`[`head`](https://rdrr.io/r/utils/head.html)`(``D2``, ``4``)``, `[`tail`](https://rdrr.io/r/utils/head.html)`(``D2``, ``1``)``)``, ``3``)``)`

    ##           w1 w2 w3 w4 w5 w6
    ##        1.000 -4  0  0  0  0
    ##       -0.250  1  0  0  0  0
    ##        1.000  0 -3  0  0  0
    ##       -0.333  0  1  0  0  0
    ## [31,]  1.000  1  1  1  1  1

### 4.3. Step 2 & 3: Gram Matrix G and Inverse Gram Matrix G^(-1)

The \\6 \times 6\\ Gram matrix \\G_2 = D_2^T D_2\\ and its inverse
\\G_2^{-1}\\ are computed below:

\
`# Gram matrix G2 = D2^T * D2`\
`G2`` ``<-`` `[`t`](https://rdrr.io/r/base/t.html)`(``D2``)`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``D2`\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Gram Matrix G2 (6x6) ---"``)`

    ## [1] "--- Gram Matrix G2 (6x6) ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``G2``, ``3``)``)`

    ##        w1     w2      w3     w4     w5     w6
    ## w1  7.347 -3.250  -2.333 -1.000 -2.333 -3.250
    ## w2 -3.250 48.132  -6.143 -2.333 -4.200 -1.000
    ## w3 -2.333 -6.143 150.000 -4.200 -4.200 -5.167
    ## w4 -1.000 -2.333  -4.200 26.040 -1.000 -2.333
    ## w5 -2.333 -4.200  -4.200 -1.000 16.191 -2.333
    ## w6 -3.250 -1.000  -5.167 -2.333 -2.333 32.139

\
`# Inverse Gram Matrix G2^(-1)`\
`G2_inv`` ``<-`` `[`solve`](https://rdrr.io/r/base/solve.html)`(``G2``)`\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Inverse Gram Matrix G2^(-1) (6x6) ---"``)`

    ## [1] "--- Inverse Gram Matrix G2^(-1) (6x6) ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``G2_inv``, ``3``)``)`

    ##       w1    w2    w3    w4    w5    w6
    ## w1 0.166 0.016 0.005 0.012 0.033 0.021
    ## w2 0.016 0.023 0.002 0.004 0.009 0.004
    ## w3 0.005 0.002 0.007 0.002 0.003 0.002
    ## w4 0.012 0.004 0.002 0.040 0.006 0.005
    ## w5 0.033 0.009 0.003 0.006 0.072 0.010
    ## w6 0.021 0.004 0.002 0.005 0.010 0.035

### 4.4. Step 4: PIGM Priority Weights & Criteria Rankings

\
`# Unnormalized row sums v = e * G2^(-1)`\
`v2`` ``<-`` `[`rowSums`](https://rdrr.io/r/base/colSums.html)`(``G2_inv``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Unnormalized Row Sum Vector v2:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``v2``, ``3``)``, ``"\n"``)`

    ## Unnormalized Row Sum Vector v2: 0.254 0.057 0.021 0.069 0.134 0.077

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Sum of Row Sums (e * G2^(-1) * e^T):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(`[`sum`](https://rdrr.io/r/base/sum.html)`(``v2``)``, ``3``)``, ``"\n\n"``)`

    ## Sum of Row Sums (e * G2^(-1) * e^T): 0.612

\
`# Normalized PIGM priority vector`\
`w_PIGM_A2`` ``<-`` `[`PIGM`](https://kkfyuen.github.io/igm4ahp/reference/PIGM.md)`(``A2``)`\
[`names`](https://rdrr.io/r/base/names.html)`(``w_PIGM_A2``)`` ``<-`` ``criteria_names`\
\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- PIGM Priority Weights & Rankings ---"``)`

    ## [1] "--- PIGM Priority Weights & Rankings ---"

\
`df_pigm_rank`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`  Weight ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(``w_PIGM_A2``, ``3``)``,`\
`  Rank   ``=`` `[`rank`](https://rdrr.io/r/base/rank.html)`(``-``w_PIGM_A2``)`\
`)`\
[`print`](https://rdrr.io/r/base/print.html)`(``df_pigm_rank``)`

    ##               Weight Rank
    ## Learning       0.415    1
    ## Friends        0.094    5
    ## School Life    0.035    6
    ## Vocational     0.112    4
    ## College Prep   0.219    2
    ## Music Classes  0.125    3

The resulting PIGM priority vector is: - **Learning**: \\0.415\\
(Rank 1) - **College Preparation**: \\0.219\\ (Rank 2) - **Music
Classes**: \\0.125\\ (Rank 3) - **Vocational Training**: \\0.112\\ (Rank
4) - **Friends**: \\0.094\\ (Rank 5) - **School Life**: \\0.035\\ (Rank
6)

### 4.5. Step 5: Unshifted NIGM (r = 0) vs. Shifted NIGM (r = 1)

Because matrix \\A_2\\ is inconsistent, the unshifted Gram matrix \\C_2
= G_2 - E = D_0^T D_0\\ is non-singular, allowing unshifted NIGM (\\r =
0\\) evaluation:

\
`# Unshifted NIGM (r = 0)`\
`res_nigm_r0`` ``<-`` `[`eNIGM`](https://kkfyuen.github.io/igm4ahp/reference/eNIGM.md)`(``A2``, r ``=`` ``0``, detail ``=`` ``TRUE``)`\
\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Unshifted Inverse Matrix C2^(-1) (r = 0) ---"``)`

    ## [1] "--- Unshifted Inverse Matrix C2^(-1) (r = 0) ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``res_nigm_r0``$``G_inv``, ``3``)``)`

    ##       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
    ## [1,] 0.333 0.053 0.019 0.057 0.121 0.072
    ## [2,] 0.053 0.032 0.005 0.014 0.029 0.015
    ## [3,] 0.019 0.005 0.008 0.006 0.011 0.006
    ## [4,] 0.057 0.014 0.006 0.052 0.030 0.019
    ## [5,] 0.121 0.029 0.011 0.030 0.118 0.036
    ## [6,] 0.072 0.015 0.006 0.019 0.036 0.050

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"\nUnnormalized Row Sums v_C (r = 0):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``res_nigm_r0``$``rs_G_inv``, ``3``)``, ``"\n"``)`

    ## 
    ## Unnormalized Row Sums v_C (r = 0): 0.655 0.148 0.055 0.177 0.346 0.198

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Sum of Row Sums (r = 0):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(`[`sum`](https://rdrr.io/r/base/sum.html)`(``res_nigm_r0``$``rs_G_inv``)``, ``3``)``, ``"\n"``)`

    ## Sum of Row Sums (r = 0): 1.579

\
`w_NIGM_r0`` ``<-`` ``res_nigm_r0``$``W`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Normalized Priority Vector W (r = 0):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``w_NIGM_r0``, ``3``)``, ``"\n"``)`

    ## Normalized Priority Vector W (r = 0): 0.415 0.094 0.035 0.112 0.219 0.125

\
`# Shifted NIGM (r = 1)`\
`w_NIGM_r1`` ``<-`` `[`eNIGM`](https://kkfyuen.github.io/igm4ahp/reference/eNIGM.md)`(``A2``, r ``=`` ``1``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Normalized Priority Vector W (r = 1):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``w_NIGM_r1``, ``3``)``, ``"\n"``)`

    ## Normalized Priority Vector W (r = 1): 0.415 0.094 0.035 0.112 0.219 0.125

### 4.6. Step 6: Lagrangian-IGM (lagrangeIGM) for r = 0 and r = 1

We evaluate `lagrangeIGM` for both \\r = 0\\ (unshifted block \\C_2\\)
and \\r = 1\\ (shifted block \\G_2\\, as presented in *Yuen, 2024*,
Section 4.2):

\
`# LIGM with r = 0`\
`res_LIGM_r0_A2`` ``<-`` `[`lagrangeIGM`](https://kkfyuen.github.io/igm4ahp/reference/lagrangeIGM.md)`(``A2``, r ``=`` ``0``, detail ``=`` ``TRUE``)`\
`w_LIGM_r0_A2`` ``<-`` ``res_LIGM_r0_A2``$``W``[``1``:``6``]`\
`lambda_r0_A2`` ``<-`` ``res_LIGM_r0_A2``$``G_inv``[``7``, ``7``]`\
\
`# LIGM with r = 1 (matching manuscript Passage 98)`\
`res_LIGM_r1_A2`` ``<-`` `[`lagrangeIGM`](https://kkfyuen.github.io/igm4ahp/reference/lagrangeIGM.md)`(``A2``, r ``=`` ``1``, detail ``=`` ``TRUE``)`\
\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Augmented KKT Inverse Gram Matrix G_bar^(-1) (7x7, r = 1) ---"``)`

    ## [1] "--- Augmented KKT Inverse Gram Matrix G_bar^(-1) (7x7, r = 1) ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``res_LIGM_r1_A2``$``G_inv``, ``3``)``)`

    ##        [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]
    ## [1,]  0.061 -0.008 -0.004 -0.017 -0.022 -0.010  0.415
    ## [2,] -0.008  0.018  0.000 -0.003 -0.003 -0.004  0.094
    ## [3,] -0.004  0.000  0.006 -0.001 -0.001 -0.001  0.035
    ## [4,] -0.017 -0.003 -0.001  0.032 -0.009 -0.004  0.112
    ## [5,] -0.022 -0.003 -0.001 -0.009  0.042 -0.007  0.219
    ## [6,] -0.010 -0.004 -0.001 -0.004 -0.007  0.025  0.125
    ## [7,]  0.415  0.094  0.035  0.112  0.219  0.125 -1.633

\
`w_LIGM_r1_A2`` ``<-`` ``res_LIGM_r1_A2``$``W``[``1``:``6``]`\
`lambda_r1_A2`` ``<-`` ``res_LIGM_r1_A2``$``G_inv``[``7``, ``7``]`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"\nExtracted Priority Weights W (r = 1):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``w_LIGM_r1_A2``, ``3``)``, ``"\n"``)`

    ## 
    ## Extracted Priority Weights W (r = 1): 0.415 0.094 0.035 0.112 0.219 0.125

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Extracted Lagrange Multiplier lambda (r = 0):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``lambda_r0_A2``, ``3``)``, ``"\n"``)`

    ## Extracted Lagrange Multiplier lambda (r = 0): -0.633

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Extracted Lagrange Multiplier lambda (r = 1):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``lambda_r1_A2``, ``3``)``, ``"(Matches manuscript Passage 98)\n"``)`

    ## Extracted Lagrange Multiplier lambda (r = 1): -1.633 (Matches manuscript Passage 98)

### 4.7. Step 7: Weighted Least Squares Error (WLSE) Comparison

We evaluate the Weighted Least Squares Error \\f(W) = \sum\_{i,j} (w_i -
a\_{ij} w_j)^2\\ using `WLSE(w, A2)`:

\
`# WLSE for IGM closed-form weights`\
`wlse_PIGM_A2`` ``<-`` `[`WLSE`](https://kkfyuen.github.io/igm4ahp/reference/WLSE.md)`(``w_PIGM_A2``, ``A2``)`\
\
`# WLSE for Saaty EVM weights`\
`w_EVM_A2`` ``<-`` `[`SaatyEigen`](https://kkfyuen.github.io/igm4ahp/reference/SaatyEigen.md)`(``A2``)`\
`wlse_EVM_A2`` ``<-`` `[`WLSE`](https://kkfyuen.github.io/igm4ahp/reference/WLSE.md)`(``w_EVM_A2``, ``A2``)`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"WLSE (PIGM / IGM Closed-Form) :"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``wlse_PIGM_A2``, ``4``)``, ``"\n"``)`

    ## WLSE (PIGM / IGM Closed-Form) : 0.6335

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"WLSE (Saaty Eigenvalue Method):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``wlse_EVM_A2``, ``4``)``, ``"\n"``)`

    ## WLSE (Saaty Eigenvalue Method): 0.8518

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"WLS Error Advantage of IGM   :"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``wlse_EVM_A2`` ``-`` ``wlse_PIGM_A2``, ``4``)``, ``"\n"``)`

    ## WLS Error Advantage of IGM   : 0.2183

------------------------------------------------------------------------

## 5. Summary & Audit Findings of Manuscript Section 4

\
`summary_table`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`  Example ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Example 1 (4x4 Consistent)"``, ``"Example 2 (6x6 Inconsistent)"``)``,`\
`  Dimension ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"4 x 4"``, ``"6 x 6"``)``,`\
`  Consistency_CR ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``cr_A1``, ``3``)``, `[`round`](https://rdrr.io/r/base/Round.html)`(``cr_A2``, ``3``)``)``,`\
`  PIGM_Top_Rank ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"C1 (0.533)"``, ``"Learning (0.415)"``)``,`\
`  WLSE_Error ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0.0000``, `[`round`](https://rdrr.io/r/base/Round.html)`(``wlse_PIGM_A2``, ``4``)``)`\
`)`\
\
`knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(``summary_table``, caption ``=`` ``"Summary of Section 4 Numerical Demonstrations"``)`

| Example | Dimension | Consistency_CR | PIGM_Top_Rank | WLSE_Error |
|:---|:---|---:|:---|---:|
| Example 1 (4x4 Consistent) | 4 x 4 | 0.000 | C1 (0.533) | 0.0000 |
| Example 2 (6x6 Inconsistent) | 6 x 6 | 0.229 | Learning (0.415) | 0.6335 |

Summary of Section 4 Numerical Demonstrations {.table
style="width:100%;"}

#### Audit Notes on Manuscript Section 4:

1.  **Example 1 Gram Matrix Typo**: In Passage 92 of the docx
    manuscript, row 3 of matrix \\G\\ contained a stray leading `'0'`
    (`0-3.25-1.524.25-1.5`). The true row 3 is
    `[-3.250, -1.500, 24.250, -1.500]`.
2.  **Example 1 Unshifted Matrix Header Typo**: In Passage 93, the
    header for \\C\\ was written as `G=C=DTD=G-1=`. The correct
    definition is \\C = D_e^T D_e = G - E\\.
3.  **Example 2 KKT Multiplier Shift**: In Passage 98, the docx reports
    \\\lambda = -1.633\\ for \\r = 1\\. The package
    `lagrangeIGM(A, r=0)` returns \\\lambda = -0.6335\\ for \\r = 0\\.
    The relationship \\\lambda_r = \lambda_0 - r = -0.6335 - 1 =
    -1.6335\\ confirms both produce identical priority weight vectors
    \\W = \[0.415, 0.094, 0.035, 0.112, 0.219, 0.125\]\\.
