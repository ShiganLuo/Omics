### STRING correction
### Input: p.m: probability matrix, thresh: a threshold, ppiString: the list of String interactions each with a specific score 
### Output: filtered adjacency matrix by the specified threshold
correctSTRING <- function(p.m, thresh, ppiString){
  library(tidyr)
  library(dplyr)
  
  # keep only String interactions > thresh
  filteredString <- ppiString %>% filter(combined_score > thresh) %>%
    select(protein1, protein2)
  
  # keep the common interactions in p.m and STRING
  adj.m <- xtabs(~., unique(filteredString))
  
  p.m <- p.m[rownames(p.m) %in% rownames(adj.m),
             colnames(p.m) %in% colnames(adj.m)]
  
  adj.m <- adj.m[rownames(p.m),colnames(p.m)]
  
  return(adj.m)  
}


### topological correction
### Input: p.m probability matrix, thresh: a threshold, corrmethod: topological method between "jaccard", "dice", "invlogweighted"
### Output: adjacency matrix filtered by the specified threshold and method
correctTopology <- function(p.m, thresh, corrmethod){
  library(igraph)
  
  g <- graph_from_adjacency_matrix(ifelse(p.m>0,1,0))
  
  neighG = similarity(g, vids=V(g), mode = c("all"),loops=F, method=corrmethod)
  if (thresh != 0){
    neighG[neighG<thresh]=0
    neighG[neighG>0]=1
  }
  adj.m = as_adjacency_matrix(g)
  adj.m = as.matrix(adj.m) * as.matrix(neighG)

  return(adj.m)
}


### semantic correction
### Input: p.m probability matrix, thresh: a threshold and
### neighGPath: semantic neighborhood matrix corresponding to a correction method among "Resnik", "Lin", "Rel", "Jiang", "Wang"
### generated from GOSemSim with "org.Hs.eg.db" OrgDb and "biological process" (BP) ontology
### Output: adjacency matrix filtered by the specified threshold and method
correctSemantic <- function(p.m, thresh, neighGPath){
  load(neighGPath)      # neighG
  cat("neighG loaded \n")
  
  ### apply semantic correction
  if (thresh != 0){
    neighG[neighG<thresh]=0
    neighG[neighG>0]=1
  }
  
  ### get the common genes between neigh and p.m
  commonGenes = sort(as.numeric(intersect(rownames(p.m), rownames(neighG))))
  
  # filter neigh for common genes
  neighG <- neighG[which(rownames(neighG) %in% commonGenes),]
  neighG <- neighG[,which(rownames(neighG) %in% commonGenes)]  
  
  # filter probability matrix for common genes
  p.m <- p.m[which(rownames(p.m) %in% commonGenes),]
  p.m <- p.m[,which(rownames(p.m) %in% commonGenes)]
  
  neighG <- neighG[sort(rownames(neighG)),]
  neighG <- neighG[,sort(colnames(neighG))]  
  
  p.m <- p.m[sort(rownames(p.m)),]
  p.m <- p.m[,sort(colnames(p.m))]  
  
  adj.m = ifelse(p.m>0,1,0) * neighG
  
  return(adj.m)
}




