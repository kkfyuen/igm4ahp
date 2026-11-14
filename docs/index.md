# igm4ahp: Inverse Gram Matrix Methods for Analytic Hierarchy Process

igm4ahp is an R package implementing exact, closed-form Inverse Gram
Matrix (IGM) prioritization methods for Pairwise Reciprocal Matrices
(PRMs) in the Analytic Hierarchy Process (AHP).

Based on the theoretical framework established in Yuen (2024), the
package derives exact global minimizers for the Weighted Least Squares
(WLS) optimization problem in \\O(n^3)\\ algebraic operations—completely
bypassing the need for non-linear iterative optimization solvers.

## Key Features

- **Exact Closed-Form Prioritization**: Solves global Weighted Least
  Squares (WLS) priority vectors directly without iteration.
- **Three IGM Algorithms**:
  - [`PIGM()`](https://kkfyuen.github.io/igm4ahp/reference/PIGM.md) /
    [`pseudoIGM()`](https://kkfyuen.github.io/igm4ahp/reference/pseudoIGM.md):
    Pseudo-Inverse Gram Matrix formulation via design matrix \\D\\.
  - [`eNIGM()`](https://kkfyuen.github.io/igm4ahp/reference/eNIGM.md):
    Element-form Normalized Inverse Gram Matrix computation.
  - [`lagrangeIGM()`](https://kkfyuen.github.io/igm4ahp/reference/lagrangeIGM.md):
    Augmented Lagrangian-form Inverse Gram Matrix solving KKT block
    matrices.
- **Saaty Benchmark Utilities**:
  - [`SaatyEigen()`](https://kkfyuen.github.io/igm4ahp/reference/SaatyEigen.md):
    Principal eigenvector prioritization (Saaty’s EVM).
  - [`ConsistencyIndex()`](https://kkfyuen.github.io/igm4ahp/reference/ConsistencyIndex.md)
    &
    [`ConsistencyRatio()`](https://kkfyuen.github.io/igm4ahp/reference/ConsistencyRatio.md):
    Automated consistency evaluations (Saaty, 1980).
  - [`WLSE()`](https://kkfyuen.github.io/igm4ahp/reference/WLSE.md):
    Objective function calculation for Weighted Least Squares Error.
- **Interactive Web Documentation Portal**: Explore live package
  vignettes at <https://kkfyuen.github.io/igm4ahp>.

## Installation

You can install the development version of igm4ahp directly from GitHub
using `remotes` or `devtools`:

\
`if`` ``(``!`[`requireNamespace`](https://rdrr.io/r/base/ns-load.html)`(``"remotes"``, quietly ``=`` ``TRUE``)``)`` ``{`\
`  `[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"remotes"``)`\
`}`\
`remotes``::`[`install_github`](https://remotes.r-lib.org/reference/install_github.html)`(``"kkfyuen/igm4ahp"``, build_vignettes ``=`` ``TRUE``)`

## Quick Start Example

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`igm4ahp`](https://kkfyuen.github.io/igm4ahp)`)`\
\
`# Define a 3x3 Pairwise Reciprocal Matrix (PRM)`\
`A`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(`\
`  ``1``,   ``2``, ``5``,`\
`  ``1``/``2``, ``1``, ``3``,`\
`  ``1``/``5``, ``1``/``3``, ``1`\
`)``, nrow ``=`` ``3``, byrow ``=`` ``TRUE``)`\
\
`# 1. Evaluate Matrix Consistency`\
`cr`` ``<-`` `[`ConsistencyRatio`](https://kkfyuen.github.io/igm4ahp/reference/ConsistencyRatio.md)`(``A``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Consistency Ratio (CR):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``cr``, ``4``)``, ``"\n"``)`\
\
`# 2. Derive Priority Weights via Pseudo-IGM (PIGM)`\
`w_pigm`` ``<-`` `[`PIGM`](https://kkfyuen.github.io/igm4ahp/reference/PIGM.md)`(``A``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"PIGM Priority Weights:\n"``)`\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``w_pigm``, ``4``)``)`\
\
`# 3. Derive Priority Weights via Element-form Normalized IGM (eNIGM)`\
`w_enigm`` ``<-`` `[`eNIGM`](https://kkfyuen.github.io/igm4ahp/reference/eNIGM.md)`(``A``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"eNIGM Priority Weights:\n"``)`\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``w_enigm``, ``4``)``)`\
\
`# 4. Compare with Saaty's Eigenvalue Method (EVM)`\
`w_evm`` ``<-`` `[`SaatyEigen`](https://kkfyuen.github.io/igm4ahp/reference/SaatyEigen.md)`(``A``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Saaty EVM Priority Weights:\n"``)`\
[`print`](https://rdrr.io/r/base/print.html)`(`[`round`](https://rdrr.io/r/base/Round.html)`(``w_evm``, ``4``)``)`\
\
`# 5. Calculate Weighted Least Squares Error (WLSE)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"PIGM WLS Error:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(`[`WLSE`](https://kkfyuen.github.io/igm4ahp/reference/WLSE.md)`(``w_pigm``, ``A``)``, ``6``)``, ``"\n"``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"EVM WLS Error: "``, `[`round`](https://rdrr.io/r/base/Round.html)`(`[`WLSE`](https://kkfyuen.github.io/igm4ahp/reference/WLSE.md)`(``w_evm``, ``A``)``, ``6``)``, ``"\n"``)`

## Vignettes and Documentation Portal

View full online documentation and step-by-step numerical demonstrations
at the interactive portal:

- Documentation Website: <https://kkfyuen.github.io/igm4ahp>

Or access compiled vignettes locally within R:

\
`# Browse all available package vignettes`\
[`vignette`](https://rdrr.io/r/utils/vignette.html)`(``package ``=`` ``"igm4ahp"``)`\
\
`# Open Section 4 numerical examples vignette`\
[`vignette`](https://rdrr.io/r/utils/vignette.html)`(``"section4_igm_examples"``, package ``=`` ``"igm4ahp"``)`

## Citation

If you use `igm4ahp` in your research, please cite both the software
release and the foundational research paper:

> Yuen, K. K. F. (2026). *igm4ahp: Inverse Gram Matrix Methods for
> Analytic Hierarchy Process* (Version v0.1.0) \[Computer software\].
> Zenodo. <https://doi.org/10.5281/zenodo.24223557>

> Yuen, K. K. F. (2024). Inverse Gram Matrix Methods for Prioritization
> in Analytic Hierarchy Process: Explainability of Weighted Least
> Squares Optimization Method. *arXiv preprint arXiv:2401.01190*.
> <https://doi.org/10.48550/arXiv.2401.01190>

## License

This package is licensed under the GPL (\>= 3) License.
