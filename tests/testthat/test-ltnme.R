test_that("wrapper of mixed-effects model works", {
  result1=ltnme(formula=~sero_cont + log(Age_at_Collection) + Country + Gender + Case_Control +
                  BF + Solid_Food + Eggs + Fish + Soy_Prod + Rye + Barley +
                  Buckwheat_Millet + (1 | Subject_ID),
                   data=ps_app,
                   test_baseline=F,
                   niter=5,
                   reffcov = 2,
                   SEED = 1,
                   save_alpha_only = T,
                   gprior_m = 100,
                   pnull = 0.5,
                   lambda = 10)
  input_data=application_data_ut
  r=0
  niter=5
  model_index=2
  lambda_fixed=10
  gprior_m=100
  cnt=input_data$cnt
  tree=input_data$tree
  yyl=seqtab2y(cnt,tree)
  Y=yyl$Y
  YL=yyl$YL
  N=nrow(Y)
  p=ncol(Y)
  K=p+1
  g=input_data$g
  Xtest=input_data$Xtest
  Xadjust=input_data$Xadjust
  grouplabel=input_data$grouplabel
  SSS=1
  result2=gibbs_crossgroup(N=N,p=p,g=g,r=r,YL=YL,Y=Y,Xtest=Xtest,Xadjust=Xadjust,grouplabel=grouplabel,c0=1,d0=0.001,c1=1,d1=0.001,niter=5,adjust=T,reff=T,SEED=SSS,save_alpha_only=F,gprior_m=gprior_m,reffcov=model_index,lambda_fixed=lambda_fixed)
  expect_identical(result1$gibbs_samples$ALPHA,result2$BETA1)
})
