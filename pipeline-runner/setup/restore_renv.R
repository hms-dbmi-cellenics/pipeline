lockfile <- renv:::renv_lockfile_read('renv.lock')
library <- Sys.getenv('RENV_LIB')

# some Bioconductor and GitHub packages aren't being installed
# we can see which ones
current <- renv:::snapshot(
  project = '.',
  library = library,
  lockfile = NULL,
  type = 'all'
)

# skip 'crossgrade' packages = same version installed but lockfile record names or info differ
diff <- renv:::renv_lockfile_diff_packages(current, lockfile)
skip <- names(diff)[diff == 'crossgrade']

# restore all except crossgrade
renv::restore(lockfile='renv.lock', library = library, exclude = skip, clean = TRUE)

# deleting renv cache
root <- renv::paths$root()
unlink(root, recursive = TRUE)

# strip debug from shared libraries
# see http://dirk.eddelbuettel.com/blog/2017/08/20/#010_stripping_shared_libraries
cmd <- sprintf("strip --strip-debug %s/*/libs/*.so", library)
system(cmd)
