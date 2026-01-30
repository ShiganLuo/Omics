### CompSR_AS.R

### CompS: Computes local entropy of a node with stochastic vector p.v
### correspondig to lines 7-12 from CompSR.R (Teschendorff)
CompS <- function(p.v){
  tmp.idx <- which(p.v>0);
  S <-  - sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )
  return(S);
}

### CompNS: 
### compute the normalised local entropy of a node with stochastic vector p.v
### corresponding to lines 14-24 from CompSR.R (Teschendorff)
CompNS <- function(p.v){
  tmp.idx <- which(p.v>0);
  if(length(tmp.idx)>1){
    S <-  - sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )/log(length(tmp.idx));
  }else { ### one degree nodes have zero entropy, avoid singularity.
    S <- 0;
  }
  return(S);
}


### CompMaxSR: Computes the maximum entropy rate of a network with adjacency matrix adj.m (assumed connected).
CompMaxSR <- function(adj.m){
  
  require(igraph);
  # find right eigenvector of adjacency matrix
  fa <- function(x,extra=NULL) {
    as.vector(adj.m %*% x)
  }
  ap.o <- arpack(fa, options=list(n=nrow(adj.m),nev=1,which="LM",maxiter=10000),sym=TRUE);
  v <- ap.o$vectors;
  lambda <- ap.o$values;
  maxSR <- log(lambda); ### maximum entropy
  return(maxSR);
}


### CompSR.ne
### compute the non-equilibrium entropies for the full networks
### corresponding to line 117 from CompSR.R (Teschendorff)
### INPUT: p.m -> proability matrix pre-calculated and saved for one sample
### OUTPUT: list containing the global entropy (as average of the normalised locals) and the normalised locals  
CompSR.ne <- function(p.m){
  NS.v <- apply(p.m,1,CompNS)
  return(list(sr=mean(NS.v),ns=NS.v))
}

### CompSR
### MAIN function for computing the global entropy normalised by the maximum
### comparable to lines 102-123 from CompSR.R (Teschendorff)
### INPUT:
### p.m: probability matrix combining already the gene expression profile exp.v and the adjacency matrix adj.m
### exp.v: a genome-wide expression vector of a sample with names labeling the genes (also annotated to entrez gene IDs) of same length as rows of adj.m.
### maxSR: optionally, the maximum entropy rate of the network.

### OUTPUT: a list of four objects
### sr: the entropy rate of the sample (normalised between 0 and 1 if maxSR was provided).
### inv: the stationary distribution of the sample.
### s: the unnormlised local entropies of the sample.
### ns: the normalised local entropies of the sample.
CompSR <- function(p.m,maxSR=NULL){
  
  require(igraph);

  fp <- function(x,extra=NULL) {
    as.vector(p.m %*% x)
  }
  fpt <- function(x,extra=NULL) {
    as.vector(t(p.m) %*% x)
  }
  ap.o <- arpack(fpt, options=list(n=nrow(p.m),nev=1,which="LM",maxiter=10000),sym=FALSE);
  invP.v <- abs(as.numeric(ap.o$vectors));
  invP.v <- invP.v/sum(invP.v);
  S.v <- apply(p.m,1,CompS);
  SR <- sum(invP.v*S.v);
  if(is.null(maxSR)==FALSE){## if provided then normalise relative to maxSR
    SR <- SR/maxSR;
  }
  
  return(list(sr=SR,s=S.v,inv=invP.v));
}
