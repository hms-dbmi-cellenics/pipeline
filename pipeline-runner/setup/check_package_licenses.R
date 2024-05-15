get_licenses <- function(package_dir = '.') {

  lock <- renv:::renv_lockfile_load(package_dir)
  licenses <- list()

  for (pkg_info in lock$Packages) {
    pkg_name <- pkg_info$Package

    licenses$`module name` <- c(licenses$`module name`, pkg_name)
    licenses$license <- c(licenses$license, packageDescription(pkg_name, fields="License"))
  }

  licenses <- as.data.frame(licenses, check.names = FALSE)
  return(licenses)
}


licenses <- get_licenses()
licenses_to_avoid <- c('AGPL-3', 'AGPL-3 | file LICENSE', 'GPL-2', 'GPL-2 | file LICENSE', 'GPL')
problematic <- licenses[licenses$license %in% licenses_to_avoid, 'module name']

# manually checked as compatible with GPL-3
safe <- c('speedglm', 'codetools', 'formatR', 'gtable', 'snow', 'stringr', 'mime', 'highr', 'knitr')

# packages that were removed as they are not needed by Cellenics
exclude <- readRDS('exclude.rds')

# packages that are requesting update to >= GPL-2
requesting <- c('Rserve', 'proxy')

installed.problematic <- problematic[!problematic %in% c(safe, exclude, requesting)]

if (length(installed.problematic)) {
  stop(
    'Investigate licenses/necessity of following packages: ',
    paste(installed.problematic, collapse = ', ')
  )
}


