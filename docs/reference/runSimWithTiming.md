# Run and Benchmark Inverse Gram Matrix Simulation

Wraps the
[`runSimIGM`](https://kkfyuen.github.io/igm4ahp/reference/runSimIGM.md)
simulation function with automated execution timing, recording the start
time, end time, and total elapsed duration, and returning/saving a
consolidated list containing simulation outcome metrics and timing
metadata.

## Usage

``` r
runSimWithTiming(
  size = 1e+06,
  seedNo = 2000,
  fileIGMData = "igmData.Rdata",
  errorFile = "igmData.csv"
)
```

## Arguments

- size:

  Integer. Total number of random PRM simulation iterations to execute
  (default: `1000000`).

- seedNo:

  Integer. Seed number for pseudo-random number generator initialization
  (default: `2000`).

- fileIGMData:

  Character string. Output RData file path to save the timing and result
  metadata (default: `"igmData.Rdata"`).

- errorFile:

  Character string. CSV file path for logging simulation results and
  discrepancies (default: `"igmData.csv"`).

## Value

A named `list` containing:

- error:

  Logical status flag returned by
  [`runSimIGM`](https://kkfyuen.github.io/igm4ahp/reference/runSimIGM.md)
  (`FALSE` if no errors found, `TRUE` if discrepancy detected).

- startTime:

  `POSIXct` timestamp recorded when execution started.

- endTime:

  `POSIXct` timestamp recorded when execution completed.

- timeSlot:

  `difftime` object representing the total elapsed execution time.

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
  # Run timed simulation and access timing metrics from result list
  res <- runSimWithTiming(size = 10, seedNo = 2000, fileIGMData = tempfile(fileext = ".Rdata"))
  print(res$timeSlot)
} # }
```
