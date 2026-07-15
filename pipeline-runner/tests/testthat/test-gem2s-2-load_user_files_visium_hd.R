# Build a minimal VisiumV2 object with a tissue image and cell polygons that
# occupy only a sub-region of the image, mirroring the Space Ranger case where
# the hires H&E is much larger than the CytAssist capture area.
make_test_segmentations <- function(image_width = 600, image_height = 500,
                                    scale_factor = 0.1,
                                    poly_x = c(50, 90), poly_y = c(120, 200)) {
  # two square cells inside the capture sub-region, in full-resolution coords
  square <- function(cell, cx, cy, r = 5) {
    data.frame(
      cell = cell,
      x = c(cx - r, cx + r, cx + r, cx - r, cx - r),
      y = c(cy - r, cy - r, cy + r, cy + r, cy - r)
    )
  }
  seg_df <- rbind(
    square("c1", poly_x[1], poly_y[1]),
    square("c2", poly_x[2], poly_y[2])
  )
  centroid_df <- data.frame(
    cell = c("c1", "c2"),
    x = poly_x,
    y = poly_y
  )

  centroids <- SeuratObject::CreateCentroids(centroid_df)
  fov <- SeuratObject::CreateFOV(
    coords = list(centroids = centroids),
    type = "centroids",
    key = "test_"
  )

  scale_factors <- list(
    spot = 1, fiducial = 1, hires = scale_factor, lowres = scale_factor / 10
  )
  class(scale_factors) <- "scalefactors"

  # a VisiumV2 extends FOV with the image + scale factors that the crop touches
  # (image = full extent, downscaled by the hires scale factor)
  vis <- methods::new(
    "VisiumV2",
    molecules = fov@molecules,
    boundaries = fov@boundaries,
    assay = fov@assay,
    key = fov@key,
    image = array(0.5, dim = c(image_height, image_width, 3)),
    scale.factors = scale_factors
  )

  # add the polygon boundaries as compact segmentations (x/y/cell in @sf.data),
  # matching the form simplify_segmentations() produces
  vis[["segmentations"]] <- SeuratObject::CreateSegmentation(seg_df, compact = TRUE)
  vis[["simplified.segmentations"]] <-
    SeuratObject::CreateSegmentation(seg_df, compact = TRUE)
  vis
}

test_that("crop_to_capture_area crops the image to the cell bounding box", {
  seg <- make_test_segmentations()
  cropped <- crop_to_capture_area(seg, margin_frac = 0)

  # bbox in full-res coords: x [45, 95], y [115, 205]; * 0.1 hires scale with
  # floor(lo) .. ceiling(hi): cols floor(4.5)=4 .. ceiling(9.5)=10  -> width 7,
  # rows floor(11.5)=11 .. ceiling(20.5)=21 -> height 11
  dims <- dim(cropped@image)
  expect_lt(dims[2], dim(seg@image)[2]) # width shrank
  expect_lt(dims[1], dim(seg@image)[1]) # height shrank
  expect_equal(dims[2], 7)
  expect_equal(dims[1], 11)
})

test_that("crop_to_capture_area shifts coordinates to the new origin", {
  seg <- make_test_segmentations()
  cropped <- crop_to_capture_area(seg, margin_frac = 0)

  poly <- cropped@boundaries[["segmentations"]]@sf.data
  cent <- cropped@boundaries[["centroids"]]@coords
  sf <- cropped@scale.factors$hires

  # polygons started at x>=45, y>=115; after cropping they sit within ~1px of
  # the new top-left origin (checked in hires pixels: coord * scale factor)
  expect_gte(min(poly[, "x"]), 0)
  expect_gte(min(poly[, "y"]), 0)
  expect_lt(min(poly[, "x"]) * sf, 2)
  expect_lt(min(poly[, "y"]) * sf, 2)

  # coord * hires still lands inside the cropped image (rows = y, cols = x)
  expect_lte(max(poly[, "x"]) * sf, dim(cropped@image)[2])
  expect_lte(max(poly[, "y"]) * sf, dim(cropped@image)[1])
  expect_lte(max(cent[, "x"]) * sf, dim(cropped@image)[2])
})

test_that("crop_to_capture_area preserves the hires scale factor", {
  seg <- make_test_segmentations()
  cropped <- crop_to_capture_area(seg, margin_frac = 0)
  expect_equal(cropped@scale.factors$hires, seg@scale.factors$hires)
})

test_that("crop_to_capture_area is a no-op when cells already fill the image", {
  # cells span the whole 600x500 image (x up to 6000, y up to 5000 at sf 0.1)
  seg <- make_test_segmentations(
    poly_x = c(100, 5900), poly_y = c(100, 4900)
  )
  cropped <- crop_to_capture_area(seg)
  expect_equal(dim(cropped@image), dim(seg@image))
  expect_identical(
    cropped@boundaries[["centroids"]]@coords,
    seg@boundaries[["centroids"]]@coords
  )
})
