test_that("multiplication works", {
  dietcovariates=c("BF","Solid_Food","Eggs","Fish","Soy_Prod","Rye","Barley","Buckwheat_Millet")
  data=make_data(ps_app,
                 'sero_cont',F,
                 as.formula(paste0('~log(Age_at_Collection)+Country+Gender+Case_Control+',paste(dietcovariates,collapse = '+'))),
                 'Subject_ID')
  app_data=application_data_ut
  app_data$tree$edge.length=NULL
  data$tree$edge.length=NULL
  expect_identical(data,app_data)
})
