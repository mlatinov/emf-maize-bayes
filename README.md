## Reproducing the analysis

### Prerequisites

- **R** (≥ 4.6.0) — [install from CRAN](https://cran.r-project.org/)
- **CmdStan** (2.39.0) — see [installation instructions](https://mc-stan.org/docs/cmdstan-guide/installation.html)
- **renv** (R package, will be installed automatically)

### Setup

Clone the repository and open it in R or RStudio:

```bash
git clone https://github.com/mlatinov/emf-maize-bayes.git
cd emf-maize-bayes
```

Open R inside the project directory. The `.Rprofile` will automatically activate `renv`. Then restore the exact package versions used:

```r
install.packages("renv")  # if not already installed
renv::restore()
```

This installs the same R package versions used to produce the results in the paper. First-time installation may take 20–40 minutes depending on your system, as some packages compile from source.

### CmdStan installation

If you don't already have CmdStan installed:

```r
library(cmdstanr)
install_cmdstan(version = "2.39.0")
```

This installs CmdStan in `~/.cmdstan` by default. Verify the installation:

```r
cmdstanr::cmdstan_version()
```

### Running the analysis

To inspect the pipeline graph and see all available targets:

```r
tar_visnetwork()
```

The analysis is a `targets` pipeline. To rebuild everything from raw data through final figures:

```r
library(targets)
tar_make("target_name")
```

### Troubleshooting

- **CmdStan not found**: Run `cmdstanr::set_cmdstan_path("/path/to/cmdstan")` or reinstall via `install_cmdstan()`.
- **Package version conflicts**: `renv::restore()` should handle these. If problems persist, delete the `renv/library/` directory and re-run `renv::restore()`.
- **Slow first-time runs**: Initial model compilation in Stan adds 30–60 seconds per model. Subsequent runs use the cached compiled binaries.

### Questions

For questions about the code, data, or methodology, open an issue at: https://github.com/mlatinov/emf-maize-bayes/issues