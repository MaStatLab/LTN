test_that("gibbs_ltn consistent", {
  cnt=matrix(rpois(10,10),2,5)
  tree=ape::rtree(5)
  colnames(cnt)=1:5
  tree$tip.label=1:5
  data=phyloseq::phyloseq(phyloseq::otu_table(cnt,taxa_are_rows = F),phyloseq::phy_tree(tree))
  g1=gibbs_ltn(data=data,niter=5)
  yyl=apply(cnt,1, function(x) {count2y(x,tree)})
  Y=do.call(rbind,lapply(yyl, function(x){x[['Y']]}))
  YL=do.call(rbind,lapply(yyl, function(x){x[['YL']]}))
  g2=gibbs_ltn(Y=Y,YL=YL,niter=5)
  expect_identical(g1,g2)
})
