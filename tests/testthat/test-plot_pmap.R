test_that("consistent sign in plot_pmap()", {
  for (i in 1:5){
    pmap=stats::rbeta(9,1,3)
    tree=ape::rtree(10)
    alpha=stats::rnorm(9)
    e=plot_pmap(pmap,tree,'main',alpha=alpha,label_nodes = which(pmap>0.5))
    expect_identical(e,NULL)
    pmap=rep(0,9)
    e=plot_pmap(pmap,tree,'main',alpha=alpha,label_nodes = which(pmap>0.5))
    expect_identical(e,NULL)
  }
})
