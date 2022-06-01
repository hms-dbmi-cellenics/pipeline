library(testthat)
library(pipeline)
print(getwd())
print(list.files("./",recursive=TRUE))
print(list.files("../",recursive=TRUE))
print(list.files("/src",recursive=TRUE))
load('../pipeline-runner/R/sysdata.rda')

test_check("pipeline")
