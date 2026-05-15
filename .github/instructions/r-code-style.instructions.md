---
name: "R Code Style & Linting"
description: "Use when: writing or editing R code in pipeline-runner/R/. Apply styler and lintr, enforce 80-char line limit, snake_case naming, double quotes, native pipes, and roxygen2 documentation."
applyTo: "pipeline-runner/R/**/*.R"
---

# R Code Style Instructions

When working with R files in `pipeline-runner/R/`, apply these rules to every edit:

## Required Workflow

After **every code change**, run these tools:

1. **Fix style issues automatically**:
   ```r
   styler::style_file("path/to/file.R")
   ```

2. **Check for remaining linting issues**:
   ```r
   lintr::lint("path/to/file.R")
   ```

3. **Fix any linting warnings** manually (see [Common Linting Issues](#common-linting-issues) below)

## Hard Rules (Non-Negotiable)

### Line Length: 80 Characters Maximum

Every line must be ≤ 80 characters. This includes comments and strings.

**If a line exceeds 80 chars, break it:**

```r
# ❌ TOO LONG (92 chars)
translated_cell_sets$cellSets <- unflatten_cell_sets(translate_cell_sets(original_cell_sets$cellSets, sample_ids_map))

# ✅ CORRECT
translated_cell_sets$cellSets <- unflatten_cell_sets(
  translate_cell_sets(original_cell_sets$cellSets, sample_ids_map)
)
```

### Naming Convention: snake_case Only

All variable and function names must use `snake_case`. **Never use camelCase.**

```r
# ❌ WRONG
fromExperimentId <- input$fromExperimentId
copyProcessedRds <- function(...) { }

# ✅ CORRECT
from_experiment_id <- input$fromExperimentId
copy_processed_rds <- function(...) { }
```

### Strings: Double Quotes Only

Use `""` for all strings, never `''`.

```r
# ❌ WRONG
message('Processing sample')
filepath <- './data/sample.txt'

# ✅ CORRECT
message("Processing sample")
filepath <- "./data/sample.txt"
```

### Pipe: Use |> Not %>%

Use the native R pipe `|>`, not magrittr's `%>%`.

```r
# ❌ WRONG
data %>% dplyr::filter(x > 5) %>% dplyr::select(name)

# ✅ CORRECT
data |> dplyr::filter(x > 5) |> dplyr::select(name)
```

### Package Calls: Always Fully Qualified

Use `package::function()` instead of `library(package)` + bare function calls.

```r
# ❌ WRONG
library(dplyr)
left_join(df1, df2)
library(Matrix)
colSums(matrix)

# ✅ CORRECT
dplyr::left_join(df1, df2)
Matrix::colSums(matrix)
```

### Comments: Above Code, Lowercase

- Comments go **above** the code they describe, not beside
- Start comments in **lowercase**
- Use `#` with a space after it

```r
# ❌ WRONG
if (!exists("sample_id")) { # Check if sample exists
  stop("Missing sample_id")  # ID is required
}

# ✅ CORRECT
# check if sample exists
if (!exists("sample_id")) {
  # id is required
  stop("Missing sample_id")
}
```

## Formatting Rules

### Function Definitions

Break at opening parenthesis if the signature is long:

```r
# ✅ Good
copy_processed_rds <- function(
  from_experiment_id,
  to_experiment_id,
  sample_ids_map,
  pipeline_config
) {
  # body
}
```

### Function Calls & Lists

Break at commas and closing brackets:

```r
# ✅ Good
config <- list(
  name = input$experimentName,
  samples = input$sampleIds,
  organism = input$organism,
  input = list(type = input$input$type),
  sampleOptions = input$sampleOptions
)

# ✅ Good
download_and_store(
  bucket = originals_bucket,
  key = s3_path,
  file_path = local_fpath,
  s3 = s3
)
```

### Operators & Spacing

- Space around `<-`, `=`, `==`, `!=`
- No space inside parentheses: `function(x)` not `function( x )`
- Space after commas: `list(a, b)` not `list(a,b)`

```r
# ✅ Correct spacing
x <- 5
result <- function(a, b) {
  a + b
}
```

## Documentation: roxygen2 Style

**All exported functions** must have roxygen2 documentation. Place above the function:

```r
#' Download user files from S3
#'
#' Retrieves user-uploaded files from the specified S3 bucket
#' and stores them locally.
#'
#' @param input The input object from the API request
#' @param pipeline_config Configuration list with S3 credentials
#' @param input_dir Local directory path for downloads (default: "/input")
#' @param prev_out Previous output passed through the pipeline
#'
#' @return List with slots:
#'   \item{data}{Empty list}
#'   \item{output}{Configuration for next pipeline step}
#'
#' @examples
#' \dontrun{
#'   result <- download_user_files(input, config)
#' }
#'
#' @export
#'
download_user_files <- function(
  input,
  pipeline_config,
  prev_out = list(),
  input_dir = "/input"
) {
  # implementation
}
```

**Internal helper functions** (not exported) can skip the `@export` tag and use a simpler format:

```r
# validate input has required fields
validate_input <- function(input) {
  # implementation
}
```

## Common Linting Issues & Fixes

### Issue: Line Too Long

**Fix**: Break at logical points (brackets, commas, operators)

```r
# lintr: Line too long (95 characters)
result <- this_is_a_very_long_function_name(argument1, argument2, argument3)

# Fixed:
result <- this_is_a_very_long_function_name(
  argument1,
  argument2,
  argument3
)
```

### Issue: Unexpected Symbol `}` or `{`

**Fix**: Check bracket placement and spacing. Correctly formatted:

```r
if (condition) {
  # code here
} else {
  # code here
}
```

### Issue: Use `<-` Not `=` for Assignment

**Fix**: Use `<-` for variable assignment (except in function arguments)

```r
# ❌ WRONG
x = 5

# ✅ CORRECT
x <- 5

# Exception: function arguments use =
function(x = 5) { }
```

### Issue: Unused Variables

**Remove** unused variables or add `# nolint: unused_import` if intentional.

```r
# ❌ Unused
library(pkg)  # Never used later

# ✅ Used
pkg::function()
```

### Issue: `<-` Not Found

**Fix**: Ensure assignments use `<-` properly

```r
# ❌ WRONG (missing left side)
<- some_value

# ✅ CORRECT
my_var <- some_value
```

## Workflow Example

You've edited `pipeline-runner/R/my-feature.R`:

```bash
# 1. Apply automatic style fixes
Rscript -e "styler::style_file('pipeline-runner/R/my-feature.R')"

# 2. Check for remaining issues
Rscript -e "lintr::lint('pipeline-runner/R/my-feature.R')"

# Output shows issues to fix manually
# Filename:Line:Col: linter: message

# 3. Fix any reported issues in the editor

# 4. Verify all tests pass
make test-file FILE=test-my-feature.R

# 5. Commit with signoff
git commit -s -m "add feature to my-feature.R"
```

## When to Ignore Rules

Use `# nolint` comments **sparingly** for legitimate exceptions:

```r
# nolint: line_length_linter
very_long_variable_name_that_cannot_be_shortened <- important_calculation()

# nolint: object_name_linter
INPUT <- read_config()  # Intentional screaming snake_case for constants
```

**But**: Always ask yourself first—can the line be broken or the name simplified?

## See Also

- [AGENTS.md](../AGENTS.md) — Full project guide
- [Test Manager agent](../agents/test-manager.agent.md) — Testing workflow
- [styler documentation](https://styler.r-lib.org/)
- [lintr documentation](https://github.com/r-lib/lintr)
- [roxygen2 documentation](https://roxygen2.r-lib.org/)
