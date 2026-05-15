---
name: "Cellenics Pipeline Codebase Guide"
description: "Instructions for understanding and developing in the Cellenics pipeline R package with emphasis on code style conventions and testing practices."
---

# Cellenics Pipeline Agent Instructions

## Project Overview

This is an R package (`pipeline-runner`) that executes dependency-managed tasks for the Cellenics single-cell analysis platform. Tasks are defined in `pipeline-runner/R/` and submitted as AWS Step Function activities.

Key directories:
- `pipeline-runner/R/` — Task definitions and utility functions
- `pipeline-runner/tests/testthat/` — Test suite (testthat framework, Edition 3)
- `pipeline-runner/data-raw/` — Data generation scripts
- `.github/agents/` — Custom AI agent definitions
- `local-runner/` — Local development environment with Docker

## Build & Test Commands

All commands use the Makefile at the repo root:

```bash
make install      # Install R and Node.js dependencies
make build        # Build Docker image (biomage-pipeline-runner)
make test         # Run all tests in Docker
make test-file FILE=test-name.R  # Run specific test file
make snapshot     # Regenerate renv.lock from installed packages
```

Tests are executed inside Docker:
```bash
docker run --entrypoint /bin/bash biomage-pipeline-runner \
  -c "R -e 'testthat::test_local()'"
```

**Development note**: `renv` manages dependencies. Install development packages (styler, roxygen2, etc.) with `renv::install()` — they are excluded from lockfile by config.

## Code Style Conventions

**All new code and edits must follow these rules.** Use `styler::style_file()` after writing code, then validate with `lintr::lint()`.

### Line Length & Formatting

- **Maximum 80 characters per line** (hard limit)
- Break long function calls at opening/closing brackets:
  ```r
  config <- list(
    name = input$experimentName,
    samples = input$sampleIds,
    organism = input$organism,
    input = list(type = input$input$type),
    sampleOptions = input$sampleOptions
  )
  ```
- Break function signatures across lines at commas:
  ```r
  copy_processed_rds <- function(
    from_experiment_id,
    to_experiment_id,
    sample_ids_map,
    pipeline_config
  ) {
  ```

### Naming & Variables

- Use **snake_case** for all variable and function names (not camelCase)
- Avoid abbreviations; prefer clarity: `sample_id` not `sid`
- Internal helper functions start with `.` or a descriptive prefix

### Pipes & Operators

- **Use native pipe `|>`** (R 4.1+), not `%>%`
- Avoid the magrittr package where possible

### Strings & Quotes

- Use **double quotes `""`** exclusively, never single quotes `''`
- Example: `message("Processing sample")` not `message('Processing sample')`

### Comments

- **Comments above code**, not beside
- Lowercase comments; do not capitalize
- Example:
  ```r
  # check if sample exists
  if (!exists("sample_id")) {
    stop("sample_id is missing")
  }
  ```
- Not: `if (!exists("sample_id")) { # Check if sample exists`
- Not: `if (!exists("sample_id")) { # Check if sample exists`

### Package Calls

- Use **fully qualified calls** `package::function()` instead of `library(package)`
- Example: `dplyr::left_join()`, `Matrix::colSums()`, `Seurat::CreateSeuratObject()`

### Documentation

- Document all exported functions with **roxygen2** docstrings (see [roxygen2 vignette](https://roxygen2.r-lib.org/))
- Example:
  ```r
  #' Download user files from S3
  #'
  #' @param input The input object from the request
  #' @param pipeline_config result of \code{load_config}
  #' @param input_dir Path where S3 object will be stored
  #'
  #' @return list with 'output' slot containing config
  #' @export
  #'
  download_user_files <- function(input, pipeline_config, ...) {
    ...
  }
  ```

## Testing & Coverage

- Write tests in `tests/testthat/test-*.R` files
- Use `describe()` blocks to group related tests
- **Target >80% code coverage** for critical paths and error handling
- Test both happy path and error conditions
- Use `mockery::mock()` or `testthat::with_mock()` to isolate dependencies
- Fixtures and test data go in `tests/testthat/mock_data/` or helper files
- See [Test Manager agent](/.github/agents/test-manager.agent.md) for testing workflow

## Git Workflow

- Write **small, focused commit messages** (preferably one line)
- **Always sign commits** with `-s` flag:
  ```bash
  git commit -s -m "add tests for scDblFinder_simple"
  ```
- Reference issue numbers in commit bodies when applicable

## Workflow for Code Changes

1. **Edit code** in `pipeline-runner/R/`
2. **Run styler**: `styler::style_file("R/myfile.R")`
3. **Run linter**: `lintr::lint("R/myfile.R")` and fix issues
4. **Write/update tests** in `tests/testthat/`
5. **Run tests locally**: `make test-file FILE=test-myfeature.R` or `make test`
6. **Commit with signoff**: `git commit -s -m "description"`

## Docker & Local Development

- Development happens in VS Code dev container (see `.devcontainer/`)
- Build with `make build` before testing — ensures Docker image is current
- Tests run in Docker to match production environment
- R version must match Dockerfile (4.2.0) for Bioconductor packages

## Dependency Management

- Add packages with `install.packages("pkg_name")`
- Update lockfile: `renv::snapshot(prompt = FALSE)`
- Commit `renv.lock` changes
- Development packages (`styler`, `roxygen2`, etc.) install with `renv::install()` — they're in DESCRIPTION but excluded from lockfile

## Key Files & Examples

- [gem2s-1-download_user_files.R](pipeline-runner/R/gem2s-1-download_user_files.R) — File I/O and S3 operations
- [test-gem2s-1-download_user_files.R](pipeline-runner/tests/testthat/test-gem2s-1-download_user_files.R) — S3 mocking patterns
- [test-qc-5-filter_doublets.R](pipeline-runner/tests/testthat/test-qc-5-filter_doublets.R) — Seurat object testing
- [helper-functions.R](pipeline-runner/tests/testthat/helper-functions.R) — Shared test utilities

## Common Pitfalls to Avoid

- ❌ Don't use `library(pkg)` → ✅ Use `pkg::function()`
- ❌ Don't write long lines (>80 chars) → ✅ Break at brackets/commas
- ❌ Don't use `%>%` → ✅ Use `|>` or explicit assignment
- ❌ Don't use single quotes for strings → ✅ Use double quotes `""`
- ❌ Don't skip linting → ✅ Run `lintr::lint()` after editing
- ❌ Don't commit without `-s` → ✅ Always add signoff flag

## When You Need Help

- Review existing code in `pipeline-runner/R/` for pattern examples
- Check test files for mocking and fixture patterns
- Consult [Test Manager agent](/.github/agents/test-manager.agent.md) for test coverage questions
- Refer to [renv documentation](https://rstudio.github.io/renv/) for dependency issues
