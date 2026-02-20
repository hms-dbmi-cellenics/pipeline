test_that("extract_subset_user_metadata extracts track names correctly from simple keys", {
  # Test: "Group-T1" should extract track="Group", value="T1"
  subset_cellsets <- data.table::as.data.table(
    list(
      key = c("Group-T1", "Group-T2"),
      name = c("T1", "T2"),
      type = c("metadata", "metadata"),
      cell_id = c(1, 2)
    )
  )

  result <- extract_subset_user_metadata(subset_cellsets)

  expect_true("Group" %in% names(result))
  expect_equal(result$Group, c("T1", "T2"))
})


test_that("extract_subset_user_metadata handles track names with hyphens", {
  # Test: "Cell-Type-Activated" should extract track="Cell-Type", value="Activated"
  subset_cellsets <- data.table::as.data.table(
    list(
      key = c("Cell-Type-Activated", "Cell-Type-Resting"),
      name = c("Activated", "Resting"),
      type = c("metadata", "metadata"),
      cell_id = c(1, 2)
    )
  )

  result <- extract_subset_user_metadata(subset_cellsets)

  expect_true("Cell-Type" %in% names(result))
  expect_equal(result$`Cell-Type`, c("Activated", "Resting"))
})


test_that("extract_subset_user_metadata handles multiple metadata tracks", {
  # Test: Multiple independent tracks
  subset_cellsets <- data.table::as.data.table(
    list(
      key = c("Group-T1", "Group-T2", "CellType-A", "CellType-B"),
      name = c("T1", "T2", "A", "B"),
      type = c("metadata", "metadata", "metadata", "metadata"),
      cell_id = c(1, 2, 3, 4)
    )
  )

  result <- extract_subset_user_metadata(subset_cellsets)

  expect_equal(length(result), 2)
  expect_true("Group" %in% names(result))
  expect_true("CellType" %in% names(result))
  expect_equal(result$Group, c("T1", "T2"))
  expect_equal(result$CellType, c("A", "B"))
})


test_that("sync_metadata_from_cellsets updates seurat metadata correctly", {
  # Create a mock Seurat object
  counts <- matrix(1, nrow = 10, ncol = 20)
  scdata <- Seurat::CreateSeuratObject(counts = counts)
  scdata$cells_id <- 1:20
  scdata$Group <- rep(c("1", "2"), 10)  # Old values

  # Create cellsets with renamed values
  cellsets <- data.table::as.data.table(
    list(
      key = c("Group-T1", "Group-T1", "Group-T2", "Group-T2"),
      name = c("T1", "T1", "T2", "T2"),  # New values
      type = c("metadata", "metadata", "metadata", "metadata"),
      cell_id = c(1, 2, 11, 12)  # Subset of cells
    )
  )

  result <- sync_metadata_from_cellsets(scdata, cellsets)

  # Check that synced cells have new values
  expect_equal(result@meta.data$Group[1], "T1")
  expect_equal(result@meta.data$Group[2], "T1")
  expect_equal(result@meta.data$Group[11], "T2")
  expect_equal(result@meta.data$Group[12], "T2")
})


test_that("sync_metadata_from_cellsets handles empty metadata cellsets", {
  counts <- matrix(1, nrow = 10, ncol = 20)
  scdata <- Seurat::CreateSeuratObject(counts = counts)
  scdata$cells_id <- 1:20

  # Empty cellsets
  cellsets <- data.table::as.data.table(
    list(
      key = character(),
      name = character(),
      type = character(),
      cell_id = integer()
    )
  )

  result <- sync_metadata_from_cellsets(scdata, cellsets)

  expect_equal(result, scdata)
})


test_that("sync_metadata_from_cellsets handles multiple tracks with renamed values", {
  counts <- matrix(1, nrow = 10, ncol = 10)
  scdata <- Seurat::CreateSeuratObject(counts = counts)
  scdata$cells_id <- 1:10
  scdata$Group <- rep(c("1", "2"), 5)
  scdata$CellType <- rep(c("A", "B"), 5)

  # Cellsets with both group and celltype renamed
  cellsets <- data.table::as.data.table(
    list(
      key = c("Group-T1", "Group-T2", "CellType-TypeX", "CellType-TypeY"),
      name = c("T1", "T2", "TypeX", "TypeY"),
      type = c("metadata", "metadata", "metadata", "metadata"),
      cell_id = c(1, 2, 3, 4)
    )
  )

  result <- sync_metadata_from_cellsets(scdata, cellsets)

  expect_equal(result@meta.data$Group[1], "T1")
  expect_equal(result@meta.data$Group[2], "T2")
  expect_equal(result@meta.data$CellType[3], "TypeX")
  expect_equal(result@meta.data$CellType[4], "TypeY")
})
