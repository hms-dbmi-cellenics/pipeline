test_that("update_sets_through_api builds a correct request", {})

with_fake_http(
  test_that("update_sets_through_api sends patch request", {
    expect_PATCH(update_sets_through_api(list(), "api_url", "experiment_id", "cell_set_key", "auth"))
  })
)
