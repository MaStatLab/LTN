test_that("makeDataParser works", {
  dietcovariates = c("BF",
                     "Solid_Food",
                     "Eggs",
                     "Fish",
                     "Soy_Prod",
                     "Rye",
                     "Barley",
                     "Buckwheat_Millet")
  data = makeDataParser(as.formula(
    paste0(
      '~sero_cont+log(Age_at_Collection)+Country+Gender+Case_Control+',
      paste(dietcovariates, collapse = '+'),
      '+(1|Subject_ID)'
    )
  ),
  ps_app,
  F)
  app_data = application_data_ut
  app_data$tree$edge.length = NULL
  data$tree$edge.length = NULL
  expect_identical(data, app_data)
})
