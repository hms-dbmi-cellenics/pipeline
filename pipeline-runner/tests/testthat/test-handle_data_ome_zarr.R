read_zarray <- function(path, level = "0") {
  jsonlite::fromJSON(file.path(path, level, ".zarray"))
}

test_that("rgb_image_to_ome_zarr stores uint8 chunks compressed with blosc/zstd", {
  path <- file.path(tempdir(), "test-image.ome.zarr")
  on.exit(unlink(path, recursive = TRUE))
  image_arr <- array(runif(64 * 48 * 3), dim = c(64, 48, 3))

  rgb_image_to_ome_zarr(image_arr, path, "test")

  for (level in as.character(0:4)) {
    zarray <- read_zarray(path, level)
    expect_equal(zarray$dtype, "|u1")
    expect_equal(zarray$compressor$id, "blosc")
    expect_equal(zarray$compressor$cname, "zstd")
  }
  expect_equal(read_zarray(path)$shape, c(3, 64, 48))
})

test_that("rgb_image_to_ome_zarr stores all 3 channels of a tile in one chunk", {
  path <- file.path(tempdir(), "test-image-chunks.ome.zarr")
  on.exit(unlink(path, recursive = TRUE))
  image_arr <- array(runif(64 * 48 * 3), dim = c(64, 48, 3))

  rgb_image_to_ome_zarr(image_arr, path, "test")

  # (c, y, x) like the array, so a chunk is a 512x512 tile with every channel
  for (level in as.character(0:4)) {
    expect_equal(read_zarray(path, level)$chunks, c(3, 512, 512))
  }
})

test_that("polygons_to_bitmask_ome_zarr stores float32 labels compressed with blosc/zstd", {
  path <- file.path(tempdir(), "test-segmentation.ome.zarr")
  on.exit(unlink(path, recursive = TRUE))
  # one 10x10 square cell with cells_id 4 (stored 1-indexed as label 5)
  polygons <- data.frame(
    cells_id = c(4, 4, 4, 4),
    x = c(5, 15, 15, 5),
    y = c(5, 5, 15, 15)
  )

  polygons_to_bitmask_ome_zarr(polygons, height = 32, width = 32, output_path = path)

  for (level in as.character(0:4)) {
    zarray <- read_zarray(path, level)
    expect_equal(zarray$dtype, "<f4")
    expect_equal(zarray$compressor$cname, "zstd")
  }

  zarr <- reticulate::import("zarr")
  labels <- reticulate::py_to_r(zarr$open(file.path(path, "0"), mode = "r")[["__getitem__"]](
    reticulate::tuple(0L)
  ))
  expect_setequal(unique(as.vector(labels)), c(0, 5))
})
