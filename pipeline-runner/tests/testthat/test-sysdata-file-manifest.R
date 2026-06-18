# Guards the camelCase-key gotcha (plan learning 5): the API camelCases the
# files map via knex (xenium_cell_boundaries -> xeniumCellBoundaries), so the
# download manifest must key file_names off the camelCased file-type names.
# Miss this and gem2s downloads nothing for the file.

test_that("xenium download manifest lists the three camelCased file types", {
  file_types <- unlist(file_types_by_technology[["xenium"]])

  expect_setequal(
    file_types,
    c("xeniumCellFeatureMatrix", "xeniumCells", "xeniumCellBoundaries")
  )
})

test_that("xenium file types resolve to the on-disk names ReadXenium expects", {
  file_types <- unlist(file_types_by_technology[["xenium"]])

  # every file type has a file_names entry (no NULL -> nothing downloaded)
  missing <- file_types[vapply(file_types, function(ft) is.null(file_names[[ft]]), logical(1))]
  expect_length(missing, 0)

  on_disk <- vapply(file_types, function(ft) file_names[[ft]], character(1))
  expect_setequal(
    on_disk,
    c("cell_feature_matrix.h5", "cells.parquet", "cell_boundaries.parquet")
  )
})
