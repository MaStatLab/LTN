#' Gibbs sampler for LTN-based mixed-effects model
#' @description a wrapper for `gibbs_crossgroup`
#' @param formula a lme4 type formula where the first fixed effect is the group indicator
#' @param data a phyloseq object
#' @param test_baseline the group to be used as baseline in the linear model
#' @param niter number of Gibbs iterations
#' @param reffcov specification of the precision matrix of random effects. 1: diagonal, 2: sparse with graphical lasso prior
#' @param SEED random seed for initializing the parameters
#' @param gprior_m hyperparameter in prior of beta with default value 100
#' @param pnull prior probability of the joint null hypothesis with default value 0.5
#' @param lambda shrinkage parameter in the graphical lasso prior with default value 10
#' @param save_alpha_only whether to only return posterior samples of alpha(A)
#' @export
ltnme = function(formula,
                 data,
                 test_baseline,
                 niter,
                 reffcov = 2,
                 SEED = 1,
                 save_alpha_only = F,
                 gprior_m = 100,
                 pnull = 0.5,
                 lambda = 10,
                verbose = F) {
  input_data=makeDataParser(formula,data,test_baseline)
  cnt=input_data$cnt
  Xtest=input_data$Xtest
  Xadjust=input_data$Xadjust
  grouplabel=input_data$grouplabel
  g=input_data$g
  tree=input_data$tree
  if(is.null(Xadjust)){
    adjust=F
  }else{
    adjust=T
  }
  yyl=seqtab2y(cnt,tree)
  Y=yyl$Y
  YL=yyl$YL
  N=nrow(Y)
  p=ncol(Y)
  K=p+1
  gibbs=gibbs_crossgroup(N=N,
                         p=p,
                         g = g,
                         r=0,
                         YL=YL,
                         Y=Y,
                         Xtest=Xtest,
                         Xadjust = Xadjust,
                         grouplabel = grouplabel,
                         c0 = 1,
                         d0 = 0.001,
                         c1 = 1,
                         d1 = 0.001,
                         nu = 3,
                         niter=niter,
                         adjust = adjust,
                         reff = T,
                         reffcov = reffcov,
                         SEED = SEED,
                         save_alpha_only = save_alpha_only,
                         gprior_m = gprior_m,
                         pnull = pnull,
                         lambda_fixed=lambda,
                        verbose = verbose)
  BETA1=lapply(gibbs$BETA1,function(x){as.vector(x)})
  BETAMAT1=do.call(rbind,BETA1)
  pmap=matrix(apply(BETAMAT1[(niter/2):niter,],2,function(x){sum(x!=0)})/(niter/2+1),nrow = 1)
  pjap=1-sum(rowSums(BETAMAT1[(niter/2):niter,]!=0)==0)/(niter/2+1)
  alpha_mean=apply(BETAMAT1[(niter/2):niter,],2,mean)
  gibbs$ALPHA=gibbs$BETA1
  gibbs$BETA=gibbs$BETA2
  rmpara=c('BETA1','Z','ACC1','ACC2','A1','A2','TAU','PHI')
  gibbs[rmpara]=NULL
  return(list(gibbs_samples=gibbs,PJAP=pjap,PMAP=pmap,alpha_mean=alpha_mean))
}


