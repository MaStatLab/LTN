
#' clr covariance of compositions generated from DTM
#' @export
clrcov_dtm_sim_log=function(nsim,tree,theta,tau,SSS=1,savesamp=F,dir=NULL){
  set.seed(SSS)
  p=length(theta)
  K=p+1
  para=rbind(theta,tau)
  logTHETA=apply(para,2,function(x){rbeta_log(nsim,x[1]*x[2],(1-x[1])*x[2])})
  print(sum(is.infinite(logTHETA)))
  print(sum(is.na(logTHETA)))
  rt=data.tree::as.Node(tree)
  maxLen=max(do.call(c,(lapply(ape::nodepath(tree),length))))
  #nodepath traversal of leaves: pre-order
  nodemat=do.call(rbind,lapply(ape::nodepath(tree),function(x){if(length(x)<maxLen){c(x,rep(rt$name,maxLen-length(x)))}else{x}}))
  nodemat=apply(nodemat, 1:2, function(x){paste0('n',x)})
  rt$Set(name=1:K,traversal = 'pre-order',filterFun = data.tree::isLeaf)
  rt$Do(function(x) {
    x$leftchild=names(x$children[1])
    x$rightchild=names(x$children[2])
  }, traversal = "pre-order",
  filterFun = data.tree::isNotLeaf)
  lc=stats::na.omit(rt$Get('leftchild'))
  rc=stats::na.omit(rt$Get('rightchild'))
  nam=c(paste0('n',lc),paste0('n',rc),paste0('n',rt$name))
  LEAFP=apply(logTHETA,1,function(x){theta2p_log(x,nam,nodemat)})
  gamma_mat=t(apply(LEAFP,2,function(x){x-(1/K)*(sum(x))}))
  if (savesamp){
    saveRDS(list(logtheta=logTHETA,leafp=LEAFP),dir)}
  return(stats::cov(gamma_mat))
}

theta2p_log=function(theta_vec,nam,nodemat){
  nodep=c(theta_vec,VGAM::log1mexp(-theta_vec),0)
  names(nodep)=nam
  leafp=apply(nodemat,1,function(x){sum(nodep[x])})
  return(leafp)
}

#' clr covariance of compositions generated from DTM
#' @export
clrcov_sim=function(mu,sigma,tree,iter){
  rt=data.tree::as.Node(tree)
  K=length(tree$tip.label)
  maxLen=max(do.call(c,(lapply(ape::nodepath(tree),length))))
  #nodepath traversal of leaves: pre-order
  nodemat=do.call(rbind,lapply(ape::nodepath(tree),function(x){if(length(x)<maxLen){c(x,rep(rt$name,maxLen-length(x)))}else{x}}))
  nodemat=apply(nodemat, 1:2, function(x){paste0('n',x)})
  rt$Set(name=1:K,traversal = 'pre-order',filterFun = data.tree::isLeaf)
  rt$Do(function(x) {
    x$leftchild=names(x$children[1])
    x$rightchild=names(x$children[2])
  }, traversal = "pre-order",
  filterFun = data.tree::isNotLeaf)
  lc=stats::na.omit(rt$Get('leftchild'))
  rc=stats::na.omit(rt$Get('rightchild'))
  nam=c(paste0('n',lc),paste0('n',rc),paste0('n',rt$name))
  set.seed(1)
  psi=mvtnorm::rmvnorm(iter,mu,sigma)
  LEAFP=apply(psi,1,function(x){psi2p(x,nam,nodemat)})
  gamma_mat=t(apply(LEAFP,2,function(x){log(x)-(1/K)*(sum(log(x)))}))
  return(stats::cov(gamma_mat))
}

psi2p=function(psi_vec,nam,nodemat){
  nodep=c(stats::plogis(psi_vec),1-stats::plogis(psi_vec),1)
  names(nodep)=nam
  leafp=apply(nodemat,1,function(x){prod(nodep[x])})
  return(leafp)
}

#' matrix loss
#' @export
matloss=function(diff){
  frob=sqrt(sum(diff^2))
  l1=max(colSums(abs(diff)))
  linf=max(abs(diff))
  spec=norm(diff,'2')
  return(c(frob,l1,linf,spec))
}

#' add 0.5 pseudo count
#' @export
add_pseudo=function(mat){
  mat[which(mat==0)]=0.5
  return(mat)
}


clrcov_sample=function(cnt){
  clrX=t(apply(cnt,1,function(x){log(x)-sum(log(x))/ncol(cnt)}))
  return(stats::cov(clrX))
}

node_binom=function(x) {
  leftchild=names(x$children[1])
  rightchild=names(x$children[2])
  x$yl <- stats::rbinom(1,x$y,x$theta_A)
  xl=x[[leftchild]]
  xl$y <- x$yl
  xr=x[[rightchild]]
  xr$y<-x$y-x$yl
}

#' generate data from DTM
#' @param theta DTM mean, preorder
#' @param tau DTM dispersion parameter, preorder
#' @export
dtm_sim=function(nsim,tree,theta,tau,total){
  p=length(theta)
  Y=matrix(0,nsim,p)
  YL=Y
  cnt=matrix(0,nsim,p+1)
  rt=data.tree::as.Node(tree)
  for (i in 1:nsim){
    theta_i=NULL
    for (j in 1:p){
      theta_i=c(theta_i,stats::rbeta(1,theta[j]*tau[j],(1-theta[j])*tau[j]))
    }
    rt$y=total
    rt$Set(theta_A=theta_i,filterFun = data.tree::isNotLeaf)
    rt$Do(node_binom, traversal = "pre-order",filterFun = data.tree::isNotLeaf)
    Y[i,]=rt$Get('y', traversal = "pre-order",filterFun = data.tree::isNotLeaf)
    YL[i,]=rt$Get('yl', traversal = "pre-order",filterFun = data.tree::isNotLeaf)
    cnt[i,]=rt$Get('y', traversal = "pre-order",filterFun = data.tree::isLeaf)
  }
  return(list(Y=Y,YL=YL,cnt=cnt))
}

leaf2nodeprob=function(leafprob,tree){
  rt=data.tree::as.Node(tree)
  rt$Set(prob=leafprob,filterFun=data.tree::isLeaf,traversal = "pre-order")
  rt$Do(function(x) {
    leftchild=names(x$children[1])
    rightchild=names(x$children[2])
    xl=x[[leftchild]]
    xr=x[[rightchild]]
    x$prob <- xl$prob+xr$prob
    x$leftprob <- xl$prob
  }, traversal = "post-order",filterFun = data.tree::isNotLeaf)
  node_prob=rt$Get('prob', traversal = "pre-order", filterFun=data.tree::isNotLeaf)
  left_prob=rt$Get('leftprob', traversal = "pre-order", filterFun=data.tree::isNotLeaf)
  return(list(node_prob=node_prob,left_prob=left_prob))
}

gibbs_glasso_cnt=function(niter,cnt,tree,r=1,s=0.01,SSS=1){
  #initialization
  set.seed(SSS)
  LAM=rep(0,niter)
  lambda=stats::rgamma(1,r,s)
  LAM[1]=lambda
  yyl=apply(cnt,1, function(x) {count2y(x,tree)})
  Y=t(do.call(rbind,lapply(yyl, function(x){x[['Y']]})))
  YL=t(do.call(rbind,lapply(yyl, function(x){x[['YL']]})))
  nsim=ncol(Y)
  kappa=NULL
  for (i in 1:nsim){
    kappa=cbind(kappa,YL[,i]-Y[,i]/2)
  }
  p=nrow(Y)
  PHI=diag(rep(1,p))
  OMEGA=list()
  set.seed(SSS)
  OMEGA[[1]]=stats::rWishart(1,p+2,PHI)[,,1]
  omega=OMEGA[[1]]
  TAU=list()
  tau=matrix(0,nrow=p,ncol=p)
  set.seed(SSS)
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
  set.seed(SSS)
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
  set.seed(SSS)
  for (i in 1:nsim){
    for (j in 1:p){
      if (Y[j,i]!=0){
        z[j,i]=BayesLogit::rpg(1,Y[j,i],0)
      }
    }
  }
  Z[[1]]=z
  set.seed(SSS)
  MU=matrix(0,ncol = niter,nrow=p)
  Lam=diag(rep(5,p))
  #gibbs sampling
  for (it in 2:niter){
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
          z[j,i]=BayesLogit::rpg(1,Y[j,i],psi[j,i])
        }
      }
    }
    Z[[it]]=z
    #update psi
    psi=matrix(0,ncol = nsim,nrow=p)
    for (i in 1:nsim){
      cov_mat=solve(omega+diag(z[,i]))
      mean_vec=cov_mat%*%(omega%*%mu+kappa[,i])
      psi[,i]=mvtnorm::rmvnorm(1,mean_vec,cov_mat)
    }
    PSI[[it]]=psi
    #update mu
    cov_mat=solve(solve(Lam)+nsim*omega)
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
      C=solve((s22+lambda)*solve(omega11)+solve(Dtau))
      s21=S[-k,k]
      beta_mean=-C%*%s21
      be=mvtnorm::rmvnorm(1,beta_mean,C)
      #(c)
      omega[-k,k]=be
      omega[k,-k]=t(be)
      omega[k,k]=gam+be%*%solve(omega[-k,-k])%*%t(be)
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
    # update lambda
    lambda=stats::rgamma(1,r+p*(p+1)/2,s+sum(abs(omega))/2)
    LAM[it]=lambda
    print(it)
  }
  return(list(OMEGA=OMEGA,LAM=LAM,MU=MU))
}

#' clr covariance of compositions generated from LTN(mu,sigma)
#' @export
clrcov_sim_log=function(mu,sigma,tree,iter,savesamp=F,dir=NULL){
  rt=data.tree::as.Node(tree)
  K=length(tree$tip.label)
  maxLen=max(do.call(c,(lapply(ape::nodepath(tree),length))))
  #nodepath traversal of leaves: pre-order
  nodemat=do.call(rbind,lapply(ape::nodepath(tree),function(x){if(length(x)<maxLen){c(x,rep(rt$name,maxLen-length(x)))}else{x}}))
  nodemat=apply(nodemat, 1:2, function(x){paste0('n',x)})
  rt$Set(name=1:K,traversal = 'pre-order',filterFun = data.tree::isLeaf)
  rt$Do(function(x) {
    x$leftchild=names(x$children[1])
    x$rightchild=names(x$children[2])
  }, traversal = "pre-order",
  filterFun = data.tree::isNotLeaf)
  lc=stats::na.omit(rt$Get('leftchild'))
  rc=stats::na.omit(rt$Get('rightchild'))
  nam=c(paste0('n',lc),paste0('n',rc),paste0('n',rt$name))
  set.seed(1)
  psi=mvtnorm::rmvnorm(iter,mu,sigma)
  LEAFP=apply(psi,1,function(x){psi2p_log(x,nam,nodemat)})
  LEAFP=LEAFP[,colSums(is.infinite(LEAFP))==0]
  if (savesamp){
    saveRDS(list(psi=psi,leafp=LEAFP),dir)
  }
  gamma_mat=t(apply(LEAFP,2,function(x){x-(1/K)*sum(x)}))
  return(stats::cov(gamma_mat))
}
rbeta_log=function(nmc,a,b){
  x=stats::rbeta(nmc,a,b)
  y=x[sapply(x,function(t){t!=0 && t!=1})]
  while (length(y)<nmc){
    y1=stats::rbeta(1,a,b)
    if (y1!=0 && y1!=1){
      y=c(y,y1)
    }
  }
  return(log(y))
}

#' transform log-odds to compositions
#' @export
psi2p_log=function(psi_vec,nam,nodemat){
  nodep=c(psi_vec-log(1+exp(psi_vec)),log(1-stats::plogis(psi_vec)),0)
  names(nodep)=nam
  leafp=apply(nodemat,1,function(x){sum(nodep[x])})
  return(leafp)
}



#' transform OTU counts to node counts (without nodemat, not very effective)
#' @param seqtab OTU table
#' @param tree -- phylogenetic tree
#' @param preorder -- whether the OTUs are ordered according to the preorder traversal of the tree
#' @return Y -- {y(A): A is internal node}, nodes are ordered according to the preorder traversal of the tree
#' @return YL -- {y(A_l): A is internal node}, nodes are ordered according to the preorder traversal of the tree
#' @export
seqtab2y=function(seqtab,tree,preorder=F){
  K=length(tree$tip.label)
  tree$node.label=as.character((K+1):(2*K-1)) # must have non-empty node labels, otherwise 'data.tree::as.Node' does not work
  rt=data.tree::as.Node(tree)
  if (!preorder){
    otu_order=rt$Get('name',filterFun=data.tree::isLeaf,traversal = "pre-order")
    seqtab_ordered=seqtab[,otu_order]
  } else {
    seqtab_ordered=seqtab
  }
  yyl=apply(seqtab_ordered,1, function(x) {count2y(x,tree)})
  Y=do.call(rbind,lapply(yyl, function(x){x[['Y']]}))
  YL=do.call(rbind,lapply(yyl, function(x){x[['YL']]}))
  return(list(Y=Y,YL=YL))
}

#' count vector to y and yl
#' @param cnt a row of OTU table (preorder)
#' @param tree phylogenetic tree
count2y=function(cnt,tree){
  K=length(tree$tip.label)
  tree$node.label=as.character((K+1):(2*K-1)) # must have non-empty node labels, otherwise 'data.tree::as.Node' does not work
  rt=data.tree::as.Node(tree)
  rt$Set(y=cnt,filterFun=data.tree::isLeaf,traversal = "pre-order")
  rt$Do(function(x) {
    leftchild=names(x$children[1])
    rightchild=names(x$children[2])
    xl=x[[leftchild]]
    xr=x[[rightchild]]
    x$yl <- xl$y
    x$y<-xl$y+xr$y
  }, traversal = "post-order",filterFun = data.tree::isNotLeaf)
  Y=rt$Get('y', traversal = "pre-order", filterFun=data.tree::isNotLeaf)
  YL=rt$Get('yl', traversal = "pre-order", filterFun=data.tree::isNotLeaf)
  return(list(Y=Y,YL=YL))
}


#' MoM estimate of Dirichlet Multinomial
#' X_i|q_i ~ Multi(N_i,q_i)
#' q_i ~ Dir(v*pi)
#' theta=1/(1+v)
#' For beta-binomials in DTM, K=2
#' @param X n samples * K categories counts
#' @export
dm_mom=function(X){
  n=nrow(X)
  K=ncol(X)
  pihat=colSums(X)/sum(X)
  Nc=(sum(X)-sum(rowSums(X)^2)/sum(X))/(n-1)
  pi_ij=t(apply(X, 1, function(x){x/sum(x)}))
  Sj=sapply(1:K,function(j){sum((pi_ij[,j]-pihat[j])^2*rowSums(X))/(n-1)})
  Gj=matrix(rowSums(X),nrow=1)%*%(pi_ij*(1-pi_ij))/sum(rowSums(X)-1)
  thetahat=sum(Sj-Gj)/(sum(Sj+(Nc-1)*Gj))
  vhat=1/thetahat-1
  return(list(pihat=pihat,thetahat=thetahat,vhat=vhat))
}
#' @export
dtm_sim=function(nsim,tree,theta,tau,total){
  p=length(theta)
  Y=matrix(0,nsim,p)
  YL=Y
  cnt=matrix(0,nsim,p+1)
  rt=data.tree::as.Node(tree)
  for (i in 1:nsim){
    theta_i=NULL
    for (j in 1:p){
      theta_i=c(theta_i,stats::rbeta(1,theta[j]*tau[j],(1-theta[j])*tau[j]))
    }
    rt$y=total
    rt$Set(theta_A=theta_i,filterFun = data.tree::isNotLeaf)
    rt$Do(node_binom, traversal = "pre-order",filterFun = data.tree::isNotLeaf)
    Y[i,]=rt$Get('y', traversal = "pre-order",filterFun = data.tree::isNotLeaf)
    YL[i,]=rt$Get('yl', traversal = "pre-order",filterFun = data.tree::isNotLeaf)
    cnt[i,]=rt$Get('y', traversal = "pre-order",filterFun = data.tree::isLeaf)
  }
  return(list(Y=Y,YL=YL,cnt=cnt))
}
#' @export
node_binom=function(x) {
  leftchild=names(x$children[1])
  rightchild=names(x$children[2])
  x$yl <- stats::rbinom(1,x$y,x$theta_A)
  xl=x[[leftchild]]
  xl$y <- x$yl
  xr=x[[rightchild]]
  xr$y<-x$y-x$yl
}

#' MoM estimate of DTM at each node
#' @export
dtm_mom=function(cnt,tree,add_pseudo=T){
  if (add_pseudo){
    cnt[cnt==0]=0.5
  }
  yyl=seqtab2y(cnt,tree,F)
  Y=yyl$Y
  YL=yyl$YL
  p=ncol(Y)
  dmmat=lapply(1:p, function(j){nodeX(j,Y,YL)})
  est=lapply(dmmat,dm_mom)
  thetahat=do.call(c,lapply(est, function(x){x[[1]][1]}))
  tauhat=do.call(c,lapply(est, function(x){x[[3]]}))
  return(list(thetahat=thetahat,tauhat=tauhat))
}

nodeX=function(j,Y,YL){
  yl=YL[,j]
  yr=Y[,j]-YL[,j]
  return(cbind(yl,yr))
}



#' evaluate clrcov when ILR(p) ~ MVN(mu,sig)
#' @param tree phylogenetic tree, has to have node labels
#' @param mu mean(preorder)
#' @param sig covariance (preorder)
#' @export
clrcov_ilr=function(tree,mu,sig){
  p=length(mu)
  K=p+1
  ilrE=diag(p)
  colnames(ilrE)=tree$node.label
  E=philr::philrInv(ilrE,tree)
  H=t(apply(E, 1, function(x){philr::clrp(x,rep(1,K))})) # philr::clrp=clr when p=1
  clrcov_true=t(H)%*%sig%*%H
  if (sum(is.na(clrcov_true))+sum(is.infinite(clrcov_true))==0){
    return(clrcov_true)
  }else{
    print('NA or inf in clrcov!')
    return(NA)
  }
}


########################
#application
##############
#' build phylogenetic tree based on taxonomy table
#' @param taxtab taxonomy table
#' @return a rooted full binary tree (some node labels are arbitrarily assigned)
#' @export
tax_tree=function(taxtab){
  taxpath=matrix(rep('ROOT',nrow(taxtab)),ncol=1)
  taxtab_otu=cbind(taxtab,rownames(taxtab))
  taxranks=colnames(taxtab)
  taxranks=c(taxranks,'OTUlevel')
  colnames(taxtab_otu)=taxranks
  set.seed(1)
  for (i in 1:length(taxranks)){
    taxmat=cbind(taxpath,taxtab_otu[,i])
    nodes=unique(taxpath)
    taxpath=apply(taxmat,1,function(x){paste(x,collapse = '/')})
    nchild=sapply(nodes, function(x){length(unique(taxmat[taxmat[,1]==x,2]))})
    children=lapply(nodes, function(x){unique(taxmat[taxmat[,1]==x,2])})
    for (j in 1:length(nchild)){
      if (nchild[j]==1){
        taxpath[taxmat[,1]==nodes[j]]=nodes[j]
      } else if (nchild[j]>2){
        randphylo=ape::rtree(n=nchild[j],tip.label = children[[j]],rooted = T)
        randNode=data.tree::as.Node(randphylo)
        randpath=ape::nodepath(randphylo) # pre-order, [ NOT the order of children[[j]] !!! ]
        leaf_preorder=randNode$Get('name',traversal = 'pre-order',filterFun = data.tree::isLeaf)
        leaf_preorder=gsub(' ','_',leaf_preorder)
        randpath_named=lapply(1:nchild[j], function(x){paste0(paste(paste0(taxranks[i],randpath[[x]][-c(1,(length(randpath[[x]])))]),collapse = '/'),'/',leaf_preorder[x])})
        randpath_named=lapply(randpath_named,function(x){gsub(paste0(taxranks[i],'/'),'',x)})
        for (k in 1:nchild[j]){
          taxpath[taxpath==paste0(nodes[j],'/',leaf_preorder[k])]=paste0(taxmat[taxpath==paste0(nodes[j],'/',leaf_preorder[k]),1],'/',randpath_named[[k]])
        }
      }
    }
  }
  pathdf=data.frame(pathString=taxpath)
  pathlist=apply(pathdf,1,function(x){strsplit(x,split='/')})
  pathdf_otu=data.frame(pathString=do.call(c,lapply(1:nrow(taxtab),function(x){paste(c(pathlist[[x]][[1]][-length(pathlist[[x]][[1]])],rownames(pathdf)[x]),collapse = '/')})))
  tree=data.tree::as.Node(pathdf_otu)
  phylotree=ape::as.phylo(tree)
  return(phylotree)
}



#' make data for application
#' @export
make_data=function(ps,test_var,test_baseline,formula_covariates,sub_var){
  cnt=phyloseq::otu_table(ps) # taxa are cols
  K=ncol(cnt)
  tree=phyloseq::phy_tree(ps)
  tree$node.label=as.character((K+1):(2*K-1))
  cnt=cnt[,tree$tip.label] # re-order OTUs
  sampdat=phyloseq::sample_data(ps)
  sub_id=as.numeric(as.factor(as.matrix(sampdat[,sub_var])))
  G=max(sub_id)
  Xtest=matrix(as.numeric(sampdat[,test_var]!=test_baseline),ncol=1)
  Xadjust=stats::model.matrix(formula_covariates,data.frame(sampdat))
  return(list(cnt=cnt,Xtest=Xtest,Xadjust=Xadjust,grouplabel=sub_id,g=G,tree=tree))
}





