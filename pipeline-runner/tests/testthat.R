library(testthat)
library(pipeline)
print(getwd())

load('../pipeline-runner/R/sysdata.rda')

test_check("pipeline")
