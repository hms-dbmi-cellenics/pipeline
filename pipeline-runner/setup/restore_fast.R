renv::restore(packages = 'pkgdepends', prompt = FALSE)
lockfile <- renv:::renv_lockfile_read('renv.lock')

# get records to install
records <- as.list(lockfile$Packages)

# exclude packages causing issues (let renv install later)
records <- records[!names(records) %in% c('SeuratWrappers')]

# convert into specs compatible with pkgdepends, and install

# functions from renv 1.0.0 that format the package specifications for pkgbuild
`%||%` <- function(a, b) {
  if (!is.null(a)) {
    return(a)
  } else {
    return(b)
  }
}

map_chr <- function (x, f, ...) {
  f <- match.fun(f)
  vapply(x, f, ..., FUN.VALUE = character(1))
}

renv_record_format_remote <- function (record) {
  remotes <- c("RemoteUsername", "RemoteRepo")
  if (all(remotes %in% names(record)))
    return(renv_record_format_short_remote(record))
  paste(record$Package, record$Version, sep = "@")
}

renv_record_format_short_remote <- function (record) {
  text <- paste(record$RemoteUsername, record$RemoteRepo, sep = "/")
  subdir <- record$RemoteSubdir %||% ""
  if (nzchar(subdir))
    text <- paste(text, subdir, sep = ":")
  if (!is.null(record$RemoteRef)) {
    ref <- record$RemoteRef
    if (!identical(ref, "master"))
      text <- paste(text, record$RemoteRef, sep = "@")
  }
  else if (!is.null(record$RemoteSha)) {
    sha <- substring(record$RemoteSha, 1L, 8L)
    text <- paste(text, sha, sep = "@")
  }
  text
}

remotes <- map_chr(records, renv_record_format_remote)

# fix explicit URLs (e.g. Matrix.utils and spatstat.core -- see below)
urls <- lapply(records, `[[`, 'RemoteUrl')
has.url <- !sapply(urls, is.null)
remotes[has.url] <- paste0('url::', unlist(urls))


# fix github refs
sources <- sapply(records, `[[`, 'Source')
is.github <- sources == 'GitHub'
remotes[is.github] <- gsub('@HEAD', '', remotes[is.github])
is.github.no.ver <- !grepl('@', remotes) & is.github
shas <- sapply(records, `[[`, 'RemoteSha')
remotes[is.github.no.ver] <- paste0(remotes[is.github.no.ver], '@', shas[is.github.no.ver])

# try all at once
# invalid version specification 'NA' was from Matrix.utils and spatstat.core
# too resolve: search for packages removed from CRAN
# avail <- available.packages()
# not.avail <- names(records)[!names(records) %in% row.names(avail)]
# look up not.avail on CRAN to see if was removed
# e.g.:
# renv::install('https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz')
# renv::snapshot()

prop.install <- pkgdepends::new_pkg_installation_proposal(
  remotes, config = list(library = Sys.getenv('RENV_LIB'), dependencies = NA), policy = 'lazy')

# will be cached for dl
prop.install$download()
prop.install$install()

# TODO: figure out why pkgdepends couldn't install some packages

# cleanup downloaded packages
cache_info <- pkgcache::pkg_cache_summary()
unlink(cache_info$cachepath, force = TRUE, recursive = TRUE)
