
#' clr covariance of compositions generated from DTM
#' @param nsim number of Monte Carlo samples
#' @param tree phylogenetic tree
#' @param theta,tau parameters of the DT distribution
#' @param SSS random seed
#' @param savesamp whether to save the Monte Carlo samples
#' @param dir if `savesamp`, the directory to save the Monte Carlo samples
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

#' clr covariance of compositions generated from LTN
#' @param mu,sigma paramters of LTN
#' @param  tree phylogenetic tree
#' @param iter number of Monte Carlo samples
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
#' @description The Frobenius, L1, L_infinity, spectral norm of the matrix `diff`
#' @export
matloss=function(diff){
  frob=sqrt(sum(diff^2))
  l1=max(colSums(abs(diff)))
  linf=max(abs(diff))
  spec=norm(diff,'2')
  return(c(frob,l1,linf,spec))
}

#' add 0.5 pseudo count
#' @description add 0.5 pseudo count to 0 counts
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
#' @param nsim number of samples
#' @param tree DTM tree
#' @param theta DTM mean
#' @param tau DTM dispersion parameter
#' @param total the total counts
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
#' @description helper function
#' @export
psi2p_log=function(psi_vec,nam,nodemat){
  nodep=c(psi_vec-log(1+exp(psi_vec)),log(1-stats::plogis(psi_vec)),0)
  names(nodep)=nam
  leafp=apply(nodemat,1,function(x){sum(nodep[x])})
  return(leafp)
}



#' transform OTU counts to node counts
#' @param seqtab OTU table
#' @param tree -- phylogenetic tree
#' @param preorder -- whether the OTUs are ordered according to the preorder traversal of the tree
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

#' transform count vector to y(A) and y(Al)
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
#' @description Method-of-moments estimates of the DM model
#' @description X_i|q_i ~ Multi(N_i,q_i), q_i ~ Dir(v*pi), theta=1/(1+v)
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
#' sample from DTM
#' @param nsim number of samples
#' @param tree DT tree
#' @param theta DT mean
#' @param tau DT dispersion parameter
#' @param total total counts in each sample
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
#' binomial experiment at a node
#' @description helper function
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
#' @param cnt OTU counts
#' @param tree DT tree
#' @param add_pseudo whether to add 0.5 pseudo counts to 0
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

#' make data for two-group comparison
#' @description helper function
#' @param ps phyloseq object
#' @param test_var the variable that defines the two contrasting groups
#' @param test_baseline the baseline group in the mixed effects model
#' @param formula_covariates formula of other fixed effects to be adjusted
#' @param sub_var random effects
#' @export
make_data=function(ps,test_var,test_baseline,formula_covariates,sub_var){
  cnt=phyloseq::otu_table(ps)
  if (phyloseq::taxa_are_rows(ps)){
    cnt=t(cnt)
  }
  K=ncol(cnt)
  tree=phyloseq::phy_tree(ps)
  tree$node.label=as.character((K+1):(2*K-1))
  cnt=cnt[,tree$tip.label] # re-order OTUs
  sampdat=phyloseq::sample_data(ps)
  sub_id=as.numeric(as.factor(as.matrix(sampdat[,sub_var])))
  G=max(sub_id)
  Xtest=matrix(as.numeric(sampdat[,test_var]!=test_baseline),ncol=1)
  if(!is.null(formula_covariates)){
      Xadjust=stats::model.matrix(formula_covariates,data.frame(sampdat))
  }
  else{
    Xadjust=NULL
  }
  return(list(cnt=cnt,Xtest=Xtest,Xadjust=Xadjust,grouplabel=sub_id,g=G,tree=tree))
}

#' plot PMAPs along the tree
#' @param pmap vector of PMAPs
#' @param tree phylogenetic tree
#' @param main.text title
#' @param alpha (optional) posterior mean of alpha
#' @param label labels of internal nodes
#' @param label_nodes a set of k nodes to be labelled by A1, ..., Ak
#' @export
plot_pmap=function(pmap,tree,main.text,alpha=NULL,label=NULL,label_nodes=NULL,tip_label=NULL){
  col_Pal=grDevices::colorRampPalette(c('white', 'red'))
  graphics::layout(t(1:2), widths=c(96,4))
  graphics::par(mar=rep(0.5, 4), oma=c(1,0.5,2,2), las=1)
  K=length(tree$tip.label)
  node_col = col_Pal(500)[as.numeric(cut(c((pmap),0,1), breaks = 500)) ]
  if (is.null(tip_label)){
    tree$tip.label=rep('',K)
  }else{
    tree$tip.label=tip_label
  }
  graphics::plot(tree,main='',show.node.label=FALSE,direction="downwards",show.tip.label = TRUE,cex.main=1,use.edge.length=F,align.tip.label=T,node.depth=2)
  graphics::mtext(main.text,side=3,cex=2.5)
  ape::nodelabels(bg=node_col,frame='none',cex=3,pch=21,col ='black')
  if (!is.null(label)){
    ape::nodelabels(text=label,frame='none',cex=1.2)
  }
  if (!is.null(label_nodes) & is.null(label)){
    label=rep('',K-1)
    label[label_nodes]=paste0('A',1:length(label_nodes))
    ape::nodelabels(text=label,frame='none',cex=1.2)
  }
  if (is.null(alpha)){
      alpha=rep(0,K-1)
  }
    edges=data.frame(tree$edge)
    edges$leftSign=''
    edges$rightSign=''
    leftBranch=sapply((K+1):(2*K-1),function(n){min(which(edges[,1]==n))})
    rightBranch=sapply((K+1):(2*K-1),function(n){max(which(edges[,1]==n))})
    edges[leftBranch,'leftSign']=c('-','+')[as.numeric(alpha>0)+1]
    edges[rightBranch,'rightSign']=c('-','+')[as.numeric(alpha<0)+1]
    for (j in seq_along(alpha)){
      if (alpha[j]==0){
        edges[leftBranch[j],'leftSign']=''
        edges[rightBranch[j],'rightSign']=''
      }
    }
    ape::edgelabels(text=edges$leftSign,frame='none',cex=1.2,adj = 1)
    ape::edgelabels(text=edges$rightSign,frame='none',cex=1.2,adj = -0.2)
  legend_image <- grDevices::as.raster(matrix(col_Pal(500), ncol=1))
  graphics::image(z=t(1:500), col=legend_image, axes=FALSE)
  graphics::mtext('PMAP',side=3,cex=1)
  graphics::axis(side=4,cex.axis=0.8,tick=T)
  # check sign
  lr=sapply((K+1):(2*K-1),function(n){which(edges[,1]==n)})
  try(if (!all(apply(rbind(lr,alpha), 2, function(x) {
    (edges[x[1], 'leftSign'] != edges[x[2], 'rightSign'] &
      edges[x[1], 'leftSign'] != '' &
      edges[x[2], 'rightSign'] != '')|
      (edges[x[1], 'leftSign'] == edges[x[2], 'rightSign'] &
       edges[x[1], 'leftSign'] ==''&
       x[3]==0)
  })))
    stop(paste0("inconsistent left and right sign", which(!(
      apply(rbind(lr,alpha), 2, function(x) {
        (edges[x[1], 'leftSign'] != edges[x[2], 'rightSign'] &
           edges[x[1], 'leftSign'] != '' &
           edges[x[2], 'rightSign'] != '')|
          (edges[x[1], 'leftSign'] == edges[x[2], 'rightSign'] &
             edges[x[1], 'leftSign'] ==''&
             x[3]==0)
      })
    )), '\n')))
  #left and right sign being exactly opposite
  alpha_sign=c('-','+')[as.numeric(alpha>0)+1]
  alpha_sign[alpha==0]=''
  try(if (!all(sapply(1:length(alpha), function(x) {
    alpha_sign[x] == edges[which(edges[, 1] == x + K)[1], 'leftSign']
  })))
    stop(paste0('inconsistent left and node sign', which((
      sapply(1:length(alpha), function(x) {
        alpha_sign[x] != edges[which(edges[, 1] == x + K)[1], 'leftSign']
      })
    )),'\n')))
  # left sign same as alpha sign
}

#' time series plot
#' @param data data.frame
#' @param facet_by subject-level binary variable
#' @export
plot_timepoints=function(data,x,y,color,shape,facet_by,shape_values=NULL,color_values=NULL,color_gradientn=NULL,xlab,ylab,color_labs,shape_labs){
  g=ggplot2::ggplot(data,ggplot2::aes_string(x=x,y=y,color=color,shape=shape))+ggplot2::geom_point()+ggplot2::theme_bw()+ggplot2::facet_wrap(stats::as.formula(paste0('~',facet_by)),nrow=2,scales='free_y',strip.position='left')+ggplot2::xlab(xlab)+ggplot2::ylab(ylab)+ggplot2::labs(color=color_labs,shape=shape_labs)
  if(!is.null(shape_values)){
    g=g+ggplot2::scale_shape_manual(values=shape_values)
  }
  if(!is.null(color_values)){
    g=g+ggplot2::scale_color_manual(values=color_values)
  }
  if (!is.null(color_gradientn)){
    g=g+ggplot2::scale_color_gradientn(colours=color_gradientn)
  }
  gt = ggplot2::ggplot_gtable(ggplot2::ggplot_build(g))
  n_ind=rowSums(table(as.matrix(data[,facet_by]),as.matrix(data[,y]))!=0)
  gt$heights[11]=n_ind[2]/n_ind[1]*gt$heights[7]
  grid::grid.draw(gt)
}

#' convert phyloseq object into the format of input of the LTN-based mixed-effects model sampler
#' @param formula a one-sided lme4 type formula, where the first variable is the group being tested
#' @param data a phyloseq object with an OTU table and a tree
#' @param test_baseline which level of the groups to be used as baseline
#' @export
makeDataParser=function(formula,data,test_baseline){
  term=attributes(stats::terms(formula))$term.labels
  isRandomEffect=grepl("\\|",term)
  test_var=term[1]
  if(sum(!isRandomEffect)==1){
    formula_covariates=NULL
  }else{
      formula_covariates=stats::as.formula(paste0('~',paste0(term[!isRandomEffect][-1],collapse = '+')))
  }
  if (sum(isRandomEffect)>1){
    warning('More than one random effects. Used the first one.')
  }
  sub_var=strsplit(term[isRandomEffect][1],split='\\| ')[[1]][2]
  return(make_data(ps=data,test_var = test_var,test_baseline = test_baseline,formula_covariates = formula_covariates,sub_var = sub_var))
}




