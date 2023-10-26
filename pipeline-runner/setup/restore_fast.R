options(renv.config.pak.enabled = TRUE)
library <- Sys.getenv('RENV_LIB')

renv::install("R.utils")
renv::restore(lockfile = 'renv.lock', library = library)
