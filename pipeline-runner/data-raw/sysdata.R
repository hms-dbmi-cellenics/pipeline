# constants used in GEM2S
gem2s <- list(
  max.edrops.fdr = 0.001,
  max.empty.counts = 100,
  max.empty.drops = 50
)

# path where dump/log files are saved
# mounted as a volume outside container to local-runner/debug
DEBUG_PATH <- "/debug"

usethis::use_data(gem2s, DEBUG_PATH, internal = TRUE, overwrite = TRUE)
