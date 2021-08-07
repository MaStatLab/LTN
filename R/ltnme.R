#' LTN-based mixed-effects model
#' @param reffcov 1: diagonal, 2: sparse
#' @param pnull prior probability of the null hypothesis
#' @param pi_only whether to return posterior samples of only pi and alpha
#' @export
ltnme = function(formula,
                 data,
                 test_baseline,
                 niter,
                 reffcov = 2,
                 SEED = 1,
                 pi_only = F,
                 gprior_m = 100,
                 pnull = 0.5,
                 lambda = 10) {
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
  return(gibbs_crossgroup(N=N,
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
                              pi_only = pi_only,
                              gprior_m = gprior_m,
                              pnull = pnull,
                              lambda_fixed=lambda))
}


