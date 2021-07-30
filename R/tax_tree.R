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
