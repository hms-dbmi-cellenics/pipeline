library(testthat)
library(pipeline)
print(getwd())
load('R/sysdata.rda')

test_check("pipeline")
