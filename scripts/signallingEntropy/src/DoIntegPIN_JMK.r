### DoIntegPIN.R
### Author: Andrew Teschendorff
### Date: 31 Dec 2013.
### This function takes as input
### (i) adj.m: an adjacency matrix representing a PPI network, with rownames/colnames of the matrix annotated to Entrez gene IDs.
### (ii) exp.m: an expression data matrix with rows labeling the genes (also annotated to entrez gene IDs.
### The output of the function is a list with following entries:
### a: the maximally connected network upon integrating the input PPI with the gene expression data matrix.
### e: the corresponding gene expression data matrix (same number of rows as $a).

DoIntegPIN <- function(adj.m,exp.m){
  require(igraph);
  
  ### find genes in network with gene expression profiles
  commonEID.v <- intersect(rownames(adj.m),rownames(exp.m));

  match(commonEID.v,rownames(exp.m)) -> map.idx;
  expPIN.m <- exp.m[map.idx,];
  match(commonEID.v,rownames(adj.m)) -> map1.idx;
  adjEXP.m <- adj.m[map1.idx,map1.idx];

  gr.o <- graph_from_adjacency_matrix(adjEXP.m,mode="undirected");
  comp.l <- components(gr.o)

  cd.v <- summary(factor(comp.l$member),maxsum=gorder(gr.o));   ### enlarge the maxsum
  mcID <- as.numeric(names(cd.v)[which.max(cd.v)]);
  maxc.idx <- which(comp.l$member==mcID);
  adjMC.m <- adjEXP.m[maxc.idx,maxc.idx];
  expMC.m <- expPIN.m[maxc.idx,];
 
  return(list(a=adjMC.m,e=expMC.m));
}






### DoIntegPIN.full
### Author: Johann M. Kraus just removed code from Andrew Teschendorff
### Date: 13 Jun 2023
### This function takes as input
### (i) adj.m: an adjacency matrix representing a PPI network, with rownames/colnames of the matrix annotated to Entrez gene IDs.
### (ii) exp.m: an expression data matrix with rows labeling the genes (also annotated to entrez gene IDs.
### The output of the function is a list with following entries:
### a: the full network upon integrating the input PPI with the gene expression data matrix.
###    this is the difference to Teschendorff, as we want all edges in the PPI not only the maximally connected network
### e: the corresponding gene expression data matrix (same number of rows as $a).

DoIntegPIN.full <- function(adj.m,exp.m){
  
  ### find genes in network with gene expression profiles
  commonEID.v <- intersect(rownames(adj.m),rownames(exp.m))
  match(commonEID.v,rownames(exp.m)) -> map.idx
  expPIN.m <- exp.m[map.idx,]
  match(commonEID.v,rownames(adj.m)) -> map1.idx
  adjEXP.m <- adj.m[map1.idx,map1.idx]
  return(list(a=adjEXP.m,e=expPIN.m))
}

### CompStochMatrix: This function computes the stochastic matrix of a sample specified by a gene expression profile exp.v in a network specified by an adjacency matrix adj.m. 
### INPUT:
### adj.m: an adjacency matrix representing a connected PPI network, with rownames/colnames of the matrix annotated to Entrez gene IDs.
### exp.v: a genome-wide expression vector of a sample with names labeling the genes (also annotated to entrez gene IDs) of same length as rows of adj.m.
### the stochastic vector of a node is independent of the node's gene expression value
CompStochMatrix <- function(exp.v,adj.m){
  ### compute stochastic matrix
  p.m <- matrix(0,nrow=length(exp.v),ncol=length(exp.v));
  rownames(p.m) <- rownames(adj.m);
  colnames(p.m) <- rownames(adj.m);
  
  for(g in 1:nrow(adj.m)){
    nn.idx <- which(adj.m[g,]==1);
    term2.v <- exp.v[nn.idx]/sum(exp.v[nn.idx]);
    p.m[g,nn.idx] <- term2.v;
  }  
  return(p.m)
}

