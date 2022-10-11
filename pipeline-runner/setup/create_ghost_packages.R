# functions that need an alternative to

# can replace theme_cowplot for Seurat DotPlot as styling not needed
theme_cowplot <- function(...) ggplot2::theme_classic()

# pbapply used by Seurat to add progress bar
# can just replace with normal apply family
pblapply <- function(X, FUN, ...) lapply(X, FUN, ...)
pbapply <- function(X, FUN, ...) apply(X, FUN, ...)
pbsapply <- function(X, FUN, ...) sapply(X, FUN, ...)

# used by DropletUtils::EmptyDrops
# can replace with built-in pseudo random number generator (slower)
generateSeedVectors <- function(nseeds, nwords = 2L) {
  res <- sample(-.Machine$integer.max:.Machine$integer.max,
                 nwords*nseeds,
                 replace = TRUE)

  res <- split(res, ceiling(seq_along(res)/nwords))
  names(res) <- NULL
  return(res)
}

alt_fnames <- c('theme_cowplot', 'pblapply', 'pbapply', 'pbsapply', 'generateSeedVectors')


install_ghost_package <- function(package) {


  library <- Sys.getenv('RENV_LIB', unset = .libPaths()[1])
  installed <- requireNamespace(package)
  if (!installed) return(NULL)

  library(package, character.only = TRUE)
  fun_names <- ls(paste0('package:', package))

  for (fun_name in fun_names) {
    fun.exists <- fun_name %in% alt_fnames

    if (!fun.exists) {
      ghost_fun <- function(...) {
        stop('removed due to licensing conflicts.')
      }

      assign(fun_name, ghost_fun)
    }
  }

  package.skeleton(package, list = fun_names, environment = environment())

  rd_files <- list.files(file.path(package, 'man'),'[.]Rd$', full.names = TRUE)
  for (rd_file in rd_files) {
    add_title(rd_file)
    escape_rd_comments(rd_file)
  }

  change_package_version(package, version = packageVersion(package))

  remove.packages(package, lib = library)
  renv::install(file.path('.', package), repos = NULL, type = 'source', library = library)
  unlink(package, recursive = TRUE)

}

# change default package version so that it meets required version
change_package_version <- function(package, version) {

  desc_path <- file.path(package, 'DESCRIPTION')
  desc <- readLines(desc_path)

  desc[desc == "Version: 1.0"] <- paste("Version:", version)
  writeLines(desc, desc_path)
}

# add title to Rd files created by package.skeleton
add_title <- function(rd_file) {
  rd <- readLines(rd_file)
  title.line <- which(rd == "\\title{")+1
  rd[title.line] <- 'pretend title'

  writeLines(rd, rd_file)
}

escape_rd_comments <- function(rd_file) {
  rd <- readLines(rd_file)
  name.line <- grep('name\\{', rd)
  rd[name.line] <- gsub('%', '\\\\%', rd[name.line])

  alias.line <- grep('alias\\{', rd)
  rd[alias.line] <- gsub('%', '\\\\%', rd[alias.line])

  writeLines(rd, rd_file)
}

# packages that are not used by Cellenics and have potentially problematic
# licenses (AGPL, GPL-2 only, BSD-4-Clause, NAIST-2003)

ghost_list <- c(
  'dqrng',
  'RhpcBLASctl',
  'cowplot',
  'ROCR',
  'gplots',
  'gtools',
  'pbapply',
  'rtracklayer',
  'rjson',
  'ggridges',
  'lme4',
  'minqa',
  'pheatmap',
  'pscl',
  'slam',
  'spdep',
  'units',
  'sf',
  'stringi'
)


for (package in ghost_list) {
  install_ghost_package(package)
}

# used by check_package_licenses.R
saveRDS(ghost_list, 'exclude.rds')

