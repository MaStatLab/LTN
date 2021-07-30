test_that("tax_tree() works", {
  tree1=tax_tree(taxtab_ut)
  expect_identical(tree1, tree_ut)
  expect_equal(ape::is.binary(tree1),T)
})
