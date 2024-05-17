# if Seurat not attached can cause errors when accessing metadata
library(Seurat)
library(zeallot)
library(tryCatchLog)
library(magrittr)
library(uuid)

# v5 is the default but making explicit
options(Seurat.object.assay.version = "v5")

# increase maxSize from the default of 500MB to 32GB
options(future.globals.maxSize = 32 * 1024 * 1024^2)

# show line numbers for tryCatchLog
options(keep.source.pkgs = TRUE,
        keep.source = TRUE)

# time stamp used for directory to store log/dump files in event of error
debug_timestamp <- format(Sys.time(), format = "%Y-%m-%d_at_%H-%M-%OS3")

for (f in list.files("R", ".R$", full.names = TRUE)) {
  source(f)
}
load("R/sysdata.rda") # constants

init()
