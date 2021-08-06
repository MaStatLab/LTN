
# INPUT:
# niter -- number of Gibbs iterations
# Y,YL -- p x N
# r,s: hyperparameters for shrinkage parameter lambda. lambda ~ Ga(r,s). The default choice is from Wang (2012).
# SEED: random seed
# lambda: 'hp' -- the Gamma hyperprior will be used; or a real number used as a fixed lambda, and r and s will be ignored
# OUTPUT:
# OMEGA: list of posterior samples of precision matrix Omega (pre-order of nodes)
# MU: p * niter, each column is a posterior sample of mu (pre-order of nodes)
# LAM: posterior samples of lambda
gibbs_glasso=function(niter,YL,Y,r=1,s=0.01,SEED=1,lambda='hp'){
  #initialization
  set.seed(SEED)
  LAM=rep(0,niter)
  if (lambda=='hp'){
    lambda=stats::rgamma(1,r,s)
  }
  LAM[1]=lambda
  nsim=ncol(Y)
  kappa=NULL
  for (i in 1:nsim){
    kappa=cbind(kappa,YL[,i]-Y[,i]/2)
  }
  p=nrow(Y)
  K=p+1
  PHI=diag(rep(1,p))
  OMEGA=list()
  set.seed(SEED)
  OMEGA[[1]]=stats::rWishart(1,p+2,PHI)[,,1]
  omega=OMEGA[[1]]
  TAU=list()
  tau=matrix(0,nrow=p,ncol=p)
  set.seed(SEED)
  for (l in 2:p){
    for (k in 1:(l-1)){
      mu_prime=sqrt(lambda^2/omega[k,l]^2)
      ukl=statmod::rinvgauss(1, mu_prime, lambda^2)
      tau[k,l]=1/ukl
      tau[l,k]=1/ukl
    }
  }
  TAU[[1]]=tau
  PSI=list()
  set.seed(SEED)
  psi=matrix(0,ncol = ncol(Y),nrow=nrow(Y))
  for (i in 1:nrow(Y)){
    for (j in 1:ncol(Y)){
      if (Y[i,j]==0){
        psi[i,j]=0
      }
      else{
        tmp=YL[i,j]/Y[i,j]
        if (tmp==0){
          psi[i,j]=-5
        }
        else{
          if (tmp==1){
            psi[i,j]=5
          }
          else{
            psi[i,j]=stats::qlogis(tmp)
          }
        }
      }
    }
  }
  PSI[[1]]=psi
  Lam=diag(rep(5,p))
  Z=list()
  z=matrix(0,ncol=nsim,nrow=p)
  set.seed(SEED)
  for (i in 1:nsim){
    for (j in 1:p){
      if (Y[j,i]!=0){
        z[j,i]=BayesLogit::rpg(1,Y[j,i]+0.0001,0)
        while (is.na(z[j,i])){
          z[j,i]=BayesLogit::rpg(1,Y[j,i]+0.0001,0)
          warning('NA in PG variable')
        }
      }
    }
  }
  Z[[1]]=z
  set.seed(SEED)
  MU=matrix(0,ncol = niter,nrow=p)
  Lam=diag(rep(5,p))
  #gibbs sampling
  for (it in 2:niter){
    # print(it)
    lambda=LAM[it-1]
    omega=OMEGA[[it-1]]
    psi=PSI[[it-1]]
    tau=TAU[[it-1]]
    mu=MU[,it-1]
    #update z
    z=matrix(0,ncol=nsim,nrow=p)
    for (i in 1:nsim){
      for (j in 1:p){
        if (Y[j,i]!=0){ #only update nodes with non-zero y_i(A)
          z[j,i]=BayesLogit::rpg(1,as.numeric(Y[j,i]),psi[j,i])
        }
      }
    }
    Z[[it]]=z
    #update psi
    psi=matrix(0,ncol = nsim,nrow=p)
    for (i in 1:nsim){
      cov_mat=chol2inv(chol(omega+diag(z[,i])))
      mean_vec=cov_mat%*%(omega%*%mu+kappa[,i])
      psi[,i]=mvtnorm::rmvnorm(1,mean_vec,cov_mat)
    }
    PSI[[it]]=psi
    #update mu
    cov_mat=chol2inv(chol(chol2inv(chol(Lam))+nsim*omega))
    mean_vec=cov_mat%*%omega%*%apply(psi,1,sum)
    MU[,it]=mvtnorm::rmvnorm(1,mean_vec,cov_mat)
    mu=MU[,it]
    #update omega and tau together (blocked gibbs)
    psi_diff=psi-matrix(rep(mu,nsim),ncol=nsim,byrow=F)
    S=psi_diff%*%t(psi_diff)
    #step1
    for (k in 1:p){
      #(b)
      s22=S[k,k]
      gam=stats::rgamma(1,nsim/2+1,(s22+lambda)/2)
      Dtau=diag(tau[-k,k])
      omega11=omega[-k,-k]
      C=chol2inv(chol((s22+lambda)*chol2inv(chol(omega11))+chol2inv(chol(Dtau))))
      s21=S[-k,k]
      beta_mean=-C%*%s21
      be=mvtnorm::rmvnorm(1,beta_mean,C)
      #(c)
      omega[-k,k]=be
      omega[k,-k]=t(be)
      omega[k,k]=gam+be%*%chol2inv(chol(omega[-k,-k]))%*%t(be)
    }
    # step2
    for (l in 2:p){
      for (k in 1:(l-1)){
        mu_prime=sqrt(lambda^2/omega[k,l]^2)
        ukl=statmod::rinvgauss(1, mu_prime, lambda^2)
        tau[k,l]=1/ukl
        tau[l,k]=1/ukl
      }
    }
    OMEGA[[it]]=omega
    TAU[[it]]=tau
    if (lambda=='hp'){
      # update lambda
      lambda=stats::rgamma(1,r+p*(p+1)/2,s+sum(abs(omega))/2)
    }
    LAM[it]=lambda
  }
  return(list(OMEGA=OMEGA,LAM=LAM,MU=MU))
}
