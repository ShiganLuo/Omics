box::use(
  igraph[graph_from_adjacency_matrix, components],
  ../util/LogUtil[log_msg]
)
#' Integrate PPI Network with Gene Expression Data
#'
#' This function aligns a Protein-Protein Interaction (PPI) network with gene 
#' expression data, filters for common genes, and extracts the Largest Connected 
#' Component (LCC).
#'
#' @param adj_m A square adjacency matrix (sparse or dense) representing the PPI 
#'   network. Rownames and colnames must be Entrez Gene IDs.
#' @param exp_m A matrix or data frame of gene expression data. Rownames must 
#'   be Entrez Gene IDs.
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item \strong{a}: The adjacency matrix of the maximally connected component (LCC).
#'   \item \strong{e}: The expression matrix subsetted to match the genes in \code{a}.
#' }
#'
#' @note 
#' Original Author: Andrew Teschendorff (31 Dec 2013).rewrite by shigan luo
#' Optimized version for large sparse matrices and modern snake_case style.
#'
#' @importFrom igraph graph_from_adjacency_matrix components
#' @export
do_integrate_pin <- function(adj_m, exp_m) {
  # 1. 寻找交集并初步过滤对齐
  common_ids <- intersect(rownames(adj_m), rownames(exp_m))

  adj_filtered <- adj_m[common_ids, common_ids, drop = FALSE]
  exp_filtered <- exp_m[common_ids, , drop = FALSE]

  # 2. 构建图并寻找连通分量
  graph_obj <- graph_from_adjacency_matrix(adj_filtered, mode = "undirected")
  comp_info <- components(graph_obj)

  # 3. 识别最大连通块 (LCC)
  max_comp_idx <- which.max(comp_info$csize)
  main_nodes_mask <- comp_info$membership == max_comp_idx

  removed_nodes <- nrow(adj_m) - sum(main_nodes_mask)
  log_msg("INFO", paste("Total number of nodes removed (including unmatched and disconnected nodes):", removed_nodes))

  # 5. 隐式返回结果列表
  list(
    a = adj_filtered[main_nodes_mask, main_nodes_mask, drop = FALSE],
    e = exp_filtered[main_nodes_mask, , drop = FALSE]
  )
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
#' @export
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
#' @export
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


#' Compute Stochastic Transition Matrix via Mass-Action Principle
#'
#' @description 
#' This function computes a stochastic transition matrix for a given sample by 
#' integrating its gene expression profile with a protein-protein interaction (PPI) 
#' network. The transition probability from node i to node j is proportional to 
#' the expression of node j, following the mass-action principle.
#'
#' @param exp_v A numeric vector of gene expression values for a single sample. 
#' Names must correspond to the Entrez IDs in \code{adj_m}.
#' @param adj_m A square adjacency matrix (binary or weighted) representing the 
#' network topology. Row and column names must be Entrez IDs.
#'
#' @return A non-negative stochastic matrix (rows sum to 1) representing the 
#' transition probabilities between genes/proteins in the network.
#' 
#' @details 
#' The calculation uses optimized matrix operations: 
#' \eqn{P_{ij} = \frac{A_{ij} \cdot E_j}{\sum_{k} A_{ik} \cdot E_k}}. 
#' Note that the transition probability \eqn{P_{ij}} is independent of the source 
#' node's own expression \eqn{E_i}, focusing instead on the availability of the 
#' target interactors.
#'
#' @examples
#' # p_m <- compute_stochastic_matrix(sample_exp, ppi_adj)
#'
#' @export
compute_stochastic_matrix <- function(exp_v, adj_m) {

  # 1. Multiply the adjacency matrix by the expression vector
  # Utilizing R's broadcasting: each column j of adj_m is multiplied by exp_v[j]
  # We transpose to ensure the orientation aligns with row-wise stochasticity
  weighted_m <- t(adj_m * exp_v)

  # 2. Calculate row sums (total weight of neighbors for each node)
  row_sums <- rowSums(weighted_m)

  # 3. Normalization: Divide each row by its sum
  # Handling potential division by zero for isolated nodes
  stoch_m <- weighted_m / ifelse(row_sums == 0, 1, row_sums)

  stoch_m
}
