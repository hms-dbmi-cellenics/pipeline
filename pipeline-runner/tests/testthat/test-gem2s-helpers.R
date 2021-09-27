test_that("check_input fails if sample number differs from metadata length", {
  input <- list(
      sampleNames = list('WT1', 'WT2'),
      metadata = list('Group' = list('WT'))
  )

  expect_error(check_input(input))
})


test_that("check_input passes if sample number is the same as metadata length", {
    input <- list(
        sampleNames = list('WT1', 'WT2'),
        metadata = list('Group' = list('WT', 'WT'))
    )

    expect_silent(check_input(input))
})


test_that("check_input passes with multiple metadata columns", {
    input <- list(
        sampleNames = list('WT1', 'WT2'),
        metadata = list('Group' = list('WT', 'WT'),
                        'Borg' = list('YES', 'NO'))
    )

    expect_silent(check_input(input))
})

