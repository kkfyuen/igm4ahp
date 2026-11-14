# Section 5: Monte Carlo Simulation 1 — Equivalence Verification of IGM Methods

## 1. Overview & Objectives

This vignette reproduces **Simulation 1** from **Section 5 (“Simulations
and Verifications”)** of *Yuen (2024)* (“Inverse Gram Matrix Methods for
Prioritization in Analytic Hierarchy Process: Explainability of Weighted
Least Squares Optimization Method”).

The full benchmark simulation dataset of 1,000,000 random Pairwise
Reciprocal Matrices (PRMs) evaluated across all four Inverse Gram Matrix
(IGM) formulations—taking 21.94 minutes (1,316.48 seconds) with 0.00%
discrepancy error recorded at 8 decimal places precision (\\\epsilon =
8\\)—is published as **Supplementary data 1** (*Yuen, 2024*) and is
available for direct download:
[Supplementary1.txt](https://kkfyuen.github.io/igm4ahp/articles/simData/Supplementary1.txt).

The primary objective of Simulation 1 is to verify the **exact numerical
equivalence** of the priority weight vectors derived by four IGM method
formulations across a massive dataset of randomly generated PRMs:

1.  **Pseudo-IGM (`PIGM` / `pseudoIGM`)**: Closed-form solution based on
    full design matrix pseudoinverse \\W = \frac{e (D^T D)^{-1}}{e (D^T
    D)^{-1} e^T}\\.
2.  **Normalized-IGM (`eNIGM(A, r = r)`)**: Element-form solution
    utilizing scalar shift parameter \\r \neq 0\\ (e.g., \\r \in
    \[-1000, 1000\] \setminus \\0\\\\).
3.  **Unshifted Lagrangian-IGM (`lagrangeIGM(A, r = 0)`)**: KKT system
    solution solving augmented block Gram matrix \\\bar{G}\\ without
    scalar shift.
4.  **Shifted Lagrangian-IGM (`lagrangeIGM(A, r = r)`)**: KKT system
    solution using parameterized shift block \\\bar{G}\_r\\.

Per the theoretical proofs in Theorems 1, 3, and 4 (*Yuen, 2024*), all
four formulations represent mathematically identical representations of
the unique global minimizer for the Weighted Least Squares (WLS)
objective. Simulation 1 confirms that floating-point matrix inversion
operations preserve this exact identity across random matrix orders \\n
\in \[3, 15\]\\ and scale ratings \\\tau \in \\1/9, \dots, 9\\\\.

------------------------------------------------------------------------

## 2. Package Setup

Load the **`igm4ahp`** package containing the simulation functions
(`runSimIGM`, `runSimWithTiming`), PRM generation utilities
(`randomPRR`, `getPRMFromRatings`), and closed-form IGM algorithms.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`igm4ahp`](https://kkfyuen.github.io/igm4ahp)`)`

------------------------------------------------------------------------

## 3. Step-by-Step Single Instance Illustration

Before executing the batch Monte Carlo loop, we demonstrate the
equivalence checking procedure on a single randomly generated PRM
instance using Algorithm 2 (*Yuen, 2024*).

### 3.1. Generating a Random PRM

We generate a random upper-triangular ratings vector for a random matrix
dimension \\n \in \[3, 15\]\\ using
**[`randomPRR()`](https://kkfyuen.github.io/igm4ahp/reference/randomPRR.md)**
and reconstruct the square PRM \\A\\ using
**[`getPRMFromRatings()`](https://kkfyuen.github.io/igm4ahp/reference/getPRMFromRatings.md)**.

\
`# Set seed for reproducible single-instance demonstration`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``42``)`\
\
`# Sample random order n and ratings vector rA from 9-point scale`\
`n_demo`` ``<-`` `[`sample`](https://rdrr.io/r/base/sample.html)`(``3``:``15``, ``1``)`\
`rA_demo`` ``<-`` `[`randomPRR`](https://kkfyuen.github.io/igm4ahp/reference/randomPRR.md)`(``n ``=`` ``n_demo``)`\
`A_demo`` ``<-`` `[`getPRMFromRatings`](https://kkfyuen.github.io/igm4ahp/reference/getPRMFromRatings.md)`(``rA_demo``)`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Random Matrix Dimension n:"``, ``n_demo``, ``"\n"``)`

    ## Random Matrix Dimension n: 3

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Upper-triangular Rating Count m:"``, `[`length`](https://rdrr.io/r/base/length.html)`(``rA_demo``)``, ``"\n\n"``)`

    ## Upper-triangular Rating Count m: 3

\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Generated Pairwise Reciprocal Matrix A (First 4x4 block) ---"``)`

    ## [1] "--- Generated Pairwise Reciprocal Matrix A (First 4x4 block) ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``A_demo``[``1``:`[`min`](https://rdrr.io/r/base/Extremes.html)`(``4``, ``n_demo``)``, ``1``:`[`min`](https://rdrr.io/r/base/Extremes.html)`(``4``, ``n_demo``)``]``, ``3``)``)`

    ##      [,1] [,2]  [,3]
    ## [1,]    1  0.2 0.111
    ## [2,]    5  1.0 2.000
    ## [3,]    9  0.5 1.000

### 3.2. Computing Priorities across IGM Variants

We evaluate priority weights using all four IGM function calls with a
random scalar shift \\r = 73\\:

\
`r_shift`` ``<-`` ``73`\
\
`# 1. Pseudo-IGM`\
`w_PIGM`` ``<-`` `[`PIGM`](https://kkfyuen.github.io/igm4ahp/reference/PIGM.md)`(``A_demo``)`\
\
`# 2. Element-form NIGM with shift r`\
`w_NIGM_r`` ``<-`` `[`eNIGM`](https://kkfyuen.github.io/igm4ahp/reference/eNIGM.md)`(``A_demo``, r ``=`` ``r_shift``)`\
\
`# 3. Lagrangian-IGM (r = 0, excluding last lambda element)`\
`res_LIGM_0`` ``<-`` `[`lagrangeIGM`](https://kkfyuen.github.io/igm4ahp/reference/lagrangeIGM.md)`(``A_demo``, r ``=`` ``0``)`\
`w_LIGM_0`` ``<-`` ``res_LIGM_0``[``1``:``n_demo``]`\
\
`# 4. Lagrangian-IGM (r = 73, excluding last lambda element)`\
`res_LIGM_r`` ``<-`` `[`lagrangeIGM`](https://kkfyuen.github.io/igm4ahp/reference/lagrangeIGM.md)`(``A_demo``, r ``=`` ``r_shift``)`\
`w_LIGM_r`` ``<-`` ``res_LIGM_r``[``1``:``n_demo``]`\
\
`# Combine into comparison table rounded to 8 decimal places (eps = 8)`\
`df_comp`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(`\
`  PIGM      ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(``w_PIGM``, ``8``)``,`\
`  eNIGM_r   ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(``w_NIGM_r``, ``8``)``,`\
`  LIGM_r0   ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(``w_LIGM_0``, ``8``)``,`\
`  LIGM_r73  ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(``w_LIGM_r``, ``8``)`\
`)`\
[`rownames`](https://rdrr.io/r/base/colnames.html)`(``df_comp``)`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"w"``, ``1``:``n_demo``)`\
\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Single Instance Priority Vector Equivalence (8 Decimal Places) ---"``)`

    ## [1] "--- Single Instance Priority Vector Equivalence (8 Decimal Places) ---"

\
[`print`](https://rdrr.io/r/base/print.html)`(``df_comp``)`

    ##          PIGM    eNIGM_r    LIGM_r0   LIGM_r73
    ## w1 0.05946621 0.05946621 0.05946621 0.05946621
    ## w2 0.58437587 0.58437587 0.58437587 0.58437587
    ## w3 0.35615792 0.35615792 0.35615792 0.35615792

\
`# Compute absolute discrepancy error er (Eq. 83; Yuen, 2024)`\
`er_single`` ``<-`` `[`sum`](https://rdrr.io/r/base/sum.html)`(`[`abs`](https://rdrr.io/r/base/MathFun.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``w_PIGM``, ``8``)`` ``-`` `[`round`](https://rdrr.io/r/base/Round.html)`(``w_NIGM_r``, ``8``)``)``)`` ``+`\
`             `[`sum`](https://rdrr.io/r/base/sum.html)`(`[`abs`](https://rdrr.io/r/base/MathFun.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``w_PIGM``, ``8``)`` ``-`` `[`round`](https://rdrr.io/r/base/Round.html)`(``w_LIGM_0``, ``8``)``)``)`` ``+`\
`             `[`sum`](https://rdrr.io/r/base/sum.html)`(`[`abs`](https://rdrr.io/r/base/MathFun.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``w_PIGM``, ``8``)`` ``-`` `[`round`](https://rdrr.io/r/base/Round.html)`(``w_LIGM_r``, ``8``)``)``)`\
`er_single`` ``<-`` `[`round`](https://rdrr.io/r/base/Round.html)`(``er_single``, ``6``)`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"\nTotal Discrepancy Error (er):"``, ``er_single``, ``"(Exact Equivalence Confirmed)\n"``)`

    ## 
    ## Total Discrepancy Error (er): 0 (Exact Equivalence Confirmed)

------------------------------------------------------------------------

## 4. Automated Monte Carlo Simulation via `runSimIGM()`

To verify robustness against numerical drift and matrix ill-conditioning
across diverse matrix sizes, Section 5 executes a Monte Carlo simulation
over random PRMs.

Here, we execute a demonstration simulation with **\\N = 10\\** sample
iterations using
**[`runSimIGM()`](https://kkfyuen.github.io/igm4ahp/reference/runSimIGM.md)**.
Results are saved to a CSV log inside the local **`./data`** directory
relative to the executing script location.

\
`# Ensure ./data directory exists in current working directory`\
`data_dir`` ``<-`` ``"./data"`\
`if`` ``(``!`[`dir.exists`](https://rdrr.io/r/base/files2.html)`(``data_dir``)``)`` ``{`\
`  `[`dir.create`](https://rdrr.io/r/base/files2.html)`(``data_dir``, recursive ``=`` ``TRUE``)`\
`}`\
\
`# Define CSV file path under ./data/`\
`csv_file`` ``<-`` `[`file.path`](https://rdrr.io/r/base/file.path.html)`(``data_dir``, ``"sim1_error_log.csv"``)`\
\
`# Execute Monte Carlo Simulation 1 for N = 10 iterations`\
`sim_error`` ``<-`` `[`runSimIGM`](https://kkfyuen.github.io/igm4ahp/reference/runSimIGM.md)`(``size ``=`` ``10``, seedNo ``=`` ``2000``, errorFile ``=`` ``csv_file``)`

    ## [1] "i:1"
    ## [1] "i:2"
    ## [1] "i:3"
    ## [1] "i:4"
    ## [1] "i:5"
    ## [1] "i:6"
    ## [1] "i:7"
    ## [1] "i:8"
    ## [1] "i:9"
    ## [1] "i:10"
    ## [1] "end"

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Simulation Error Status (FALSE = No errors detected):"``, ``sim_error``, ``"\n"``)`

    ## Simulation Error Status (FALSE = No errors detected): FALSE

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Simulation log saved to:"``, ``csv_file``, ``"\n"``)`

    ## Simulation log saved to: ./data/sim1_error_log.csv

### 4.1. Reading and Printing Simulation Results

We read the generated simulation log from `./data/sim1_error_log.csv`
and display the full table of results for all \\N = 10\\ iterations:

\
`# Read simulation CSV output from ./data directory`\
`sim_data`` ``<-`` `[`read.csv`](https://rdrr.io/r/utils/read.table.html)`(``csv_file``, sep ``=`` ``";"``, stringsAsFactors ``=`` ``FALSE``)`\
[`colnames`](https://rdrr.io/r/base/colnames.html)`(``sim_data``)`` ``<-`` `[`trimws`](https://rdrr.io/r/base/trimws.html)`(`[`colnames`](https://rdrr.io/r/base/colnames.html)`(``sim_data``)``)`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Total Simulation Iterations Logged:"``, `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``sim_data``)``, ``"\n"``)`

    ## Total Simulation Iterations Logged: 10

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Maximum Discrepancy Error Recorded:"``, `[`max`](https://rdrr.io/r/base/Extremes.html)`(``sim_data``$``Error``)``, ``"\n\n"``)`

    ## Maximum Discrepancy Error Recorded: 0

\
`# Display complete simulation table for N = 10`\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Complete Simulation Results (N = 10) ---"``)`

    ## [1] "--- Complete Simulation Results (N = 10) ---"

\
`knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(``sim_data``, caption ``=`` ``"Simulation 1 Logged Output (./data/sim1_error_log.csv)"``)`

| ID | Error | r | A | n | w | WLSE |
|---:|---:|---:|:---|---:|:---|---:|
| 1 | 0 | 259 | 1/2,7,6,7,1/5,1/6,1/4,1,1/3,1/5,6,3,9,1/9,1,3,1/7,1/9,1/3,7,1/8 | 7 | 0.1815,0.2077,0.1014,0.0745,0.0877,0.1941,0.1532 | 8.6145 |
| 2 | 0 | -19 | 9,1/7,7,2,9,2,8,1,1/5,1/4,2,1/7,9,1/9,1/6,1/8,1/8,7,3,5,1/4,1/9,4,1/3,1/8,1,7,7,1/4,4,6,2,3,1/2,1/3,1/7 | 9 | 0.1519,0.0368,0.0734,0.0523,0.2765,0.0538,0.0485,0.0718,0.2351 | 5.3040 |
| 3 | 0 | -758 | 1/4,1/5,1/5,1/9,1/4,1/5 | 4 | 0.0773,0.0548,0.1805,0.6874 | 0.5871 |
| 4 | 0 | -3 | 1/4,1/6,1/5,9,1/3,7,7,1/2,1/5,3,1/6,1/7,1/4,2,7,1/4,9,1/3,1/2,1/5,1/6,1/5,8,1/5,1,5,1/8,1/8,1/6,1/5,1/6,3,1/9,1/6,1/5,1/6 | 9 | 0.1273,0.0909,0.125,0.1721,0.0406,0.0714,0.1067,0.1129,0.1531 | 7.8688 |
| 5 | 0 | -904 | 1,2,1/2,1/7,1/3,4,1/9,1/4,1/9,5 | 5 | 0.1384,0.0995,0.0909,0.3581,0.313 | 4.5056 |
| 6 | 0 | -359 | 1/7,1/7,1/8,1/6,1/2,2,7,1/2,1/3,1,6,1/3,5,9,9,2,5,4,1/3,6,1/7,6,1/2,8,4,5,3,2 | 8 | 0.0379,0.1418,0.5326,0.0578,0.0844,0.0579,0.0512,0.0364 | 1.5874 |
| 7 | 0 | 546 | 1/9,1/7,1/6,1/3,1/3,1/8,1,1/2,7,1/9,6,7,1/4,1/7,1/7,3,7,1/5,1/7,9,2,8,1/8,4,1,1/7,1/3,5,7,6,9,4,7,1/8,2,7,1/4,5,5,1/9,1/3,1/3,6,8,9,1/9,1/5,1/2,4,8,5,1/5,8,5,1/3,1/2,5,1/7,8,4,2,1/6,1,4,4,1/9,1/7,4,7,1/4,4,1/6,1/9,9,6,1/2,1/3,3,2,3,4,9,9,7,5,1/2,1/4,1/8,1/9,4,2,1/7,1/5,1/7,3,2,2,1/8,1/5,1/2,1/9,1/7,8,5,1/3 | 15 | 0.0381,0.0599,0.1555,0.0596,0.0768,0.0604,0.0634,0.071,0.0423,0.0404,0.1027,0.041,0.079,0.0548,0.0553 | 9.8067 |
| 8 | 0 | -836 | 1/8,1/9,1/5,1/8,9,3,8,2,5,1/6,1/2,8,1/3,3,4,6,1/8,2,6,1/8,9,1/5,1/2,1/3,1/4,2,1/4,3,1/9,9,1/2,9,1/9,7,9,7 | 9 | 0.062,0.2186,0.0914,0.1318,0.1261,0.0571,0.1965,0.0701,0.0465 | 8.8554 |
| 9 | 0 | -131 | 4,8,1/3,1/7,1/5,1/3,6,8,1/2,1/8,6,8,7,1/6,5,8,1/7,1/5,4,1/6,1/8,7,7,2,1/2,3,1/9,9,1/9,6,3,1,1/2,1/4,2,1/7 | 9 | 0.1421,0.1178,0.0633,0.1791,0.0656,0.1132,0.0735,0.1471,0.0985 | 9.6026 |
| 10 | 0 | -272 | 1,2,1/8,2,1/9,7,4,1/9,1/2,1/8,1/9,1,9,4,1/6,1/4,9,1/6,8,1/4,1,1,1/4,1/4,1/9,1/2,1/3,1/8,1/9,7,1/6,2,1/2,1/9,1/3,1,1,1/6,5,1/3,4,1/2,7,4,1/8 | 10 | 0.0482,0.06,0.1008,0.054,0.0525,0.0808,0.104,0.3453,0.0921,0.0624 | 6.7143 |

Simulation 1 Logged Output (./data/sim1_error_log.csv) {.table}

------------------------------------------------------------------------

## 5. Timed Execution & Performance Benchmarking via `runSimWithTiming()`

Because closed-form IGM methods require only direct matrix inversion
(\\\mathcal{O}(n^3)\\ operations), batch execution is exceptionally
fast. Package wrapper function
**[`runSimWithTiming()`](https://kkfyuen.github.io/igm4ahp/reference/runSimWithTiming.md)**
records execution timestamps and logs results to `./data/`.

\
`# Define output file paths under ./data/`\
`rdata_timing_file`` ``<-`` `[`file.path`](https://rdrr.io/r/base/file.path.html)`(``data_dir``, ``"sim1_timing_data.RData"``)`\
`csv_timing_file``   ``<-`` `[`file.path`](https://rdrr.io/r/base/file.path.html)`(``data_dir``, ``"sim1_timing_log.csv"``)`\
\
`# Execute timed simulation for N = 10 iterations demo`\
`timing_results`` ``<-`` `[`runSimWithTiming`](https://kkfyuen.github.io/igm4ahp/reference/runSimWithTiming.md)`(`\
`  size ``=`` ``10``,`\
`  seedNo ``=`` ``2024``,`\
`  fileIGMData ``=`` ``rdata_timing_file``,`\
`  errorFile ``=`` ``csv_timing_file`\
`)`

    ## [1] "i:1"
    ## [1] "i:2"
    ## [1] "i:3"
    ## [1] "i:4"
    ## [1] "i:5"
    ## [1] "i:6"
    ## [1] "i:7"
    ## [1] "i:8"
    ## [1] "i:9"
    ## [1] "i:10"
    ## [1] "end"

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"--- Simulation 1 Execution Timing (N = 10) ---\n"``)`

    ## --- Simulation 1 Execution Timing (N = 10) ---

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Start Time       :"``, `[`format`](https://rdrr.io/r/base/format.html)`(``timing_results``$``startTime``)``, ``"\n"``)`

    ## Start Time       : 2026-11-14 21:46:24

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"End Time         :"``, `[`format`](https://rdrr.io/r/base/format.html)`(``timing_results``$``endTime``)``, ``"\n"``)`

    ## End Time         : 2026-11-14 21:46:24

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Elapsed Duration :"``, `[`round`](https://rdrr.io/r/base/Round.html)`(`[`as.numeric`](https://rdrr.io/r/base/numeric.html)`(``timing_results``$``timeSlot``, units ``=`` ``"secs"``)``, ``4``)``, ``"seconds\n"``)`

    ## Elapsed Duration : 0.0091 seconds

\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Discrepancy Flag :"``, ``timing_results``$``error``, ``"\n\n"``)`

    ## Discrepancy Flag : FALSE

\
`# Read and print timing log results from ./data directory`\
`timing_data`` ``<-`` `[`read.csv`](https://rdrr.io/r/utils/read.table.html)`(``csv_timing_file``, sep ``=`` ``";"``, stringsAsFactors ``=`` ``FALSE``)`\
[`colnames`](https://rdrr.io/r/base/colnames.html)`(``timing_data``)`` ``<-`` `[`trimws`](https://rdrr.io/r/base/trimws.html)`(`[`colnames`](https://rdrr.io/r/base/colnames.html)`(``timing_data``)``)`\
\
[`print`](https://rdrr.io/r/base/print.html)`(``"--- Timed Simulation Log Results (./data/sim1_timing_log.csv) ---"``)`

    ## [1] "--- Timed Simulation Log Results (./data/sim1_timing_log.csv) ---"

\
`knitr``::`[`kable`](https://rdrr.io/pkg/knitr/man/kable.html)`(``timing_data``, caption ``=`` ``"Timed Simulation 1 Logged Output"``)`

| ID | Error | r | A | n | w | WLSE |
|---:|---:|---:|:---|---:|:---|---:|
| 1 | 0 | 23 | 1/5,5,9,5,1,1/8 | 4 | 0.2032,0.6832,0.0548,0.0587 | 1.3282 |
| 2 | 0 | -488 | 3,8,6,1/8,4,1/9,3,1/7,1,6,9,8,3,2,1/6,1/3,8,1/7,1/8,9,9,1/7,1/9,1/4,1,7,5,8,3,1/8,9,4,1/5,8,1/4,1/7,1/9,3,1,1/7,3,6,8,1/6,6,9,1/2,1,8,3,1/8,4,9,1/8,1,8,5,1/4,1/6,9,7,8,5,1/3,1/7,1/6,1/3,3,6,3,1,7,1/8,1/2,2,1/8,9,6,1/8,3,1/9,1/8,1/7,5,1/5,1/5,1/8,9,1/6,8,1,1/8,5,2,1/9,1/9,1/9,1/2,1/2,8,2,1/5,1/5,3,1/3 | 15 | 0.0805,0.0523,0.0556,0.0801,0.0652,0.1146,0.0524,0.0304,0.0573,0.0524,0.0468,0.1144,0.0597,0.0915,0.0468 | 11.4184 |
| 3 | 0 | 49 | 6,1,1/6,4,1/6,1/7,1/5,9,6,1/5 | 5 | 0.0971,0.0481,0.6574,0.0811,0.1163 | 1.3845 |
| 4 | 0 | 95 | 1/2,9,5,1/4,1/9,4,1,2,1/3,1/9,1,4,9,1/9,1/5,1/5,7,1,1/9,5,4,4,1/3,2,3,3,1/8,8,1/5,1/3,2,6,4,8,1/9,3,3,1/6,1,1/8,1/2,1/3,1/3,1/9,6 | 10 | 0.1069,0.0641,0.0757,0.0965,0.1233,0.0804,0.0828,0.0673,0.1678,0.1352 | 8.9017 |
| 5 | 0 | 330 | 7,4,1/7,7,1/6,7,9,8,1/7,3,6,1/2,1/9,1/7,8 | 6 | 0.1857,0.1739,0.2072,0.0804,0.1194,0.2335 | 11.7394 |
| 6 | 0 | 194 | 1/9,4,5,3,1/3,1/9,1/5,1/8,7,7,1/6,2,6,7,1/4,1,4,1/7,9,5,6,1/9,7,8,1/7,5,3,8,3,8,9,1/2,1/3,5,1/6,5,4,1/6,4,1/4,1/2,1/2,1/7,2,2,1/3,4,1/8,6,8,1/3,7,4,1,1/6,1/3,1/4,8,1,9,1/2,1/6,9,3,1/9,1,6,9,1/6,6,4,1/8,1/9,1/9,1/7,4,1/3,1/4,2,1/7,1/3,1/5,4,9,5,1/2,2,6,9,5,3 | 14 | 0.0492,0.074,0.1112,0.0744,0.0875,0.0824,0.094,0.0332,0.0444,0.1063,0.0652,0.0694,0.0593,0.0495 | 9.6184 |
| 7 | 0 | -4 | 1/5,8,1/3,1,3,8,3,1/8,7,2,1/7,9,1/2,1/9,1/6,4,1/7,4,1/2,1,1/8,1/4,1/2,1/5,6,1/5,7,3,1/3,1/2,1/2,6,1/8,1/3,2,2,1/9,1/7,1/9,1/9,1/6,3,1/7,7,5,1/3,1/3,6,5,2,1/9,1/5,1/7,2,1/5,9,1/3,8,1/4,1/5,1/7,1/2,1/7,8,7,1/3 | 12 | 0.0867,0.0579,0.0586,0.0624,0.057,0.0891,0.0655,0.1883,0.0622,0.0888,0.0933,0.0901 | 9.2510 |
| 8 | 0 | 696 | 1/4,8,1/4,1/9,1/3,2,4,1/2,1,1/3,1/4,7,6,1/7,1/8,5,2,2,1/9,1,8,1,3,7,1/7,1/9,6,1/3,3,4,1/9,1/8,1/5,1/2,5,1/8,9,1/8,1/3,1/8,1/2,4,8,8,1/6,1/6,4,7,1/9,1,1/8,2,1/5,2,1/4,9,1/8,6,1/5,9,2,1/4,8,1/5,8,5 | 12 | 0.0815,0.0701,0.0563,0.0437,0.1138,0.0614,0.1666,0.0972,0.0882,0.0982,0.0824,0.0406 | 9.4327 |
| 9 | 0 | 763 | 5,8,2,5,5,7 | 4 | 0.6963,0.1746,0.0817,0.0474 | 0.6113 |
| 10 | 0 | 443 | 1/3,9,1/4,1/6,2,1/4 | 4 | 0.2175,0.1308,0.0774,0.5743 | 2.5113 |

Timed Simulation 1 Logged Output {.table}

------------------------------------------------------------------------

## 6. Summary of Published Benchmark Results (Yuen, 2024)

The findings from Section 5 Simulation 1 of *Yuen (2024)* are summarized
below:

| Metric / Parameter | Paper Benchmark (Simulation 1) | Vignette Demo Configuration |
|:---|:---|:---|
| **Total Sample Size (\\N\\)** | \\1,000,000\\ random PRMs | \\N = 10\\ demo iterations |
| **Log Directory Path** | Local file storage / `Supplementary data 1` | `./data/sim1_error_log.csv` |
| **Matrix Dimension (\\n\\)** | Random \\n \in \[3, 15\]\\ | Validated via [`randomPRR()`](https://kkfyuen.github.io/igm4ahp/reference/randomPRR.md) |
| **Shift Parameter (\\r\\)** | Random \\r \in \[-1000, 1000\] \setminus \\0\\\\ | Passed to `eNIGM` & `lagrangeIGM` |
| **Precision Threshold (\\\epsilon\\)** | \\8\\ decimal places (\\\epsilon = 8\\) | Evaluated via `round(w, 8)` |
| **Discrepancy Error Rate** | **\\0 / 1,000,000\\ (0.00%)** | Verified exact zero error |
| **Total Execution Time** | \\21.94\\ minutes (\\1,316.5\\ seconds) | Direct \\\mathcal{O}(n^3)\\ matrix inversion |

#### **Key Conclusions**:

1.  **Mathematical Identicalness**: `PIGM`, `eNIGM`, and `lagrangeIGM`
    produce identical priority weight vectors across \\100\\\\ of test
    cases.
2.  **Invariance to Regularization Shift \\r\\**: The scalar parameter
    \\r \neq 0\\ cancels out completely during normalization in `eNIGM`
    and `lagrangeIGM`, having zero effect on priority weights, rankings,
    or consistency indices.
3.  **Computational Efficiency**: Direct \\\mathcal{O}(n^3)\\ matrix
    inversion processes over 45,000 matrix instances per minute,
    establishing high performance for large-scale decision models.
