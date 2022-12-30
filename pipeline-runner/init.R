# if Seurat not attached can cause errors when accessing metadata
library(Seurat)
library(zeallot)
library(tryCatchLog)
library(magrittr)

# increase maxSize from the default of 500MB to 32GB
options(future.globals.maxSize = 32 * 1024 * 1024^2)

# show line numbers for tryCatchLog
options(keep.source.pkgs = TRUE)

for (f in list.files("R", ".R$", full.names = TRUE)) {
  source(f, keep.source = TRUE)
}
load("R/sysdata.rda") # constants

init()
