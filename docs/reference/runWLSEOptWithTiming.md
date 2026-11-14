# Run and Benchmark Non-Linear WLS Optimization Simulation

Wraps the
[`runWLSEOpt`](https://kkfyuen.github.io/igm4ahp/reference/runWLSEOpt.md)
simulation function with automated execution timing, recording the start
time, end time, and total elapsed duration, and returning/saving a
consolidated list containing simulation outcome metrics and timing
metadata.

## Usage

``` r
runWLSEOptWithTiming(
  size = 10000,
  seedNo = 2000,
  filewlseOptData = "wlseOptData.Rdata",
  sigDigits = 4,
  errorFile = "wlseOptData.txt"
)
```

## Arguments

- size:

  Integer. Total number of random PRM simulation instances to evaluate
  (default: `10000`).

- seedNo:

  Integer. Seed number for pseudo-random number generator initialization
  (default: `2000`).

- filewlseOptData:

  Character string. Output RData file path to save timing and result
  metadata (default: `"wlseOptData.Rdata"`).

- sigDigits:

  Integer. Number of significant decimal places used for rounding
  priority weight vectors during equality verification (default: `4`).

- errorFile:

  Character string. Output text file path for logging simulation results
  and discrepancy metrics (default: `"wlseOptData.txt"`).

## Value

A named `list` containing:

- error:

  Logical status flag returned by
  [`runWLSEOpt`](https://kkfyuen.github.io/igm4ahp/reference/runWLSEOpt.md)
  (`FALSE` if no errors found, `TRUE` if discrepancy detected).

- startTime:

  `POSIXct` timestamp recorded when execution started.

- endTime:

  `POSIXct` timestamp recorded when execution completed.

- timeSlot:

  `difftime` object representing total elapsed execution duration.

## References

Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization in
Analytic Hierarchy Process: Explainability of Weighted Least Squares
Optimization Method. *arXiv preprint arXiv:2401.01190*.
<https://doi.org/10.48550/arXiv.2401.01190>

## Author

Kevin Kam Fung Yuen

## Examples

``` r
if (FALSE) { # \dontrun{
  # Run timed WLS optimization benchmark simulation
  tmp <- tempfile(fileext = ".Rdata")
  res <- runWLSEOptWithTiming(
    size = 100,
    seedNo = 2000,
    filewlseOptData = tmp
  )
  print(res$timeSlot)
} # }
```
