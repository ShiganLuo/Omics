### time to simulate behaviour for changes on PIN
### reduce PIN
reducePIN <- function(adj.m, perc, sn=123){
  set.seed(sn)
  allEdges <- which(lower.tri(adj.m) & adj.m==1)
  rmEdges <- sample(allEdges,perc*length(allEdges))
  adj.m[rmEdges] <- 0
  
  return(list(adj.m=adj.m*t(adj.m),rmEdges=rmEdges))
}

### enlarge PIN
enlargePIN <- function(adj.m, perc, sn=123){
  set.seed(sn)
  allEdges <- which(lower.tri(adj.m) & adj.m==1)
  nonEdges <- which(lower.tri(adj.m) & adj.m==0)
  addEdges <- sample(nonEdges,perc*length(allEdges))
  adj.m[addEdges] <- 1
  
  adj.m <- 1 - adj.m
  adj.m <- adj.m*t(adj.m)
  
  return(list(adj.m=(1 - adj.m),addEdges=addEdges))
}
 
### rewire PIN
rewirePIN <- function(adj.m, perc, sn=123){
  set.seed(sn)
  
  allEdges <- which(lower.tri(adj.m) & adj.m==1, arr.ind = TRUE)
  rwEdges <- sample(allEdges,round(perc*nrow(allEdges)))     # edges to rewire
  
  for(i in rwEdges){
    sourceNode <- allEdges[i, 1]
    nonNeighs <-  which(adj.m[sourceNode,1:sourceNode-1]==0)
    targetNode <- sample(nonNeighs,1)
    
    # rewire on whole matrix
    adj.m[sourceNode, allEdges[i, 2]] <- 0
    adj.m[allEdges[i, 2], sourceNode] <- 0
    adj.m[sourceNode, targetNode] <- 1
    adj.m[targetNode, sourceNode] <- 1
  }
  
  return(list(adj.m=adj.m,rwEdges=rwEdges))
}

### flip PIN
flipPIN <- function(adj.m, perc, sn=123){
  set.seed(sn)
  rmEdges <- which(lower.tri(adj.m) & adj.m==1)
  rmEdges <- sample(rmEdges,perc*length(rmEdges))     # edges to be removed
  
  addEdges <- which(lower.tri(adj.m) & adj.m==0)
  addEdges <- sample(addEdges,perc*length(addEdges))  # edges to be added
  
  # first remove 
  adj.m[rmEdges] <- 0
  adj.m <- adj.m*t(adj.m)
  
  # then add
  adj.m[addEdges] <- 1
  
  adj.m <- 1 - adj.m
  adj.m <- adj.m*t(adj.m)
  
  return(list(adj.m=(1 - adj.m),flEdges=c(rmEdges,addEdges)))
}


