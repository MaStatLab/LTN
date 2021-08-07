
#' Gibbs sampler for LTN
#' @param data phyloseq object containing OTU table and phylogenetic tree
#' @param Y (numeric) matrix of y_i(A). Y and YL will be used if the phyloseq object is not provided
#' @param YL (numeric) matrix of y_i(A_l)
#' @param niter number of Gibbs iterations
#' @param SEED random seed used in initializing the sampler
#' @param lambda shrinkage parameter in graphical lasso prior
#' @export
gibbs_ltn=function(data=NULL,Y=NULL,YL=NULL,lambda=10,niter,SEED=1){
  if (!is.null(data)){
    # compute y(A) and y(Al)
    tree=phyloseq::phy_tree(data)
    cnt=phyloseq::otu_table(data)
    if (phyloseq::taxa_are_rows(data)){
      cnt=t(cnt)
    }
    cnt=cnt[,tree$tip.label]
    yyl=apply(cnt,1, function(x) {count2y(x,tree)})
    Y=t(do.call(rbind,lapply(yyl, function(x){x[['Y']]})))
    YL=t(do.call(rbind,lapply(yyl, function(x){x[['YL']]})))
  }
  return(gibbs_glasso(niter=niter,YL=YL,Y=Y,r=1,s=0.01,SEED=SEED,lambda=lambda))
}
