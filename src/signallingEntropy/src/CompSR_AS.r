### CompSR_AS.R
libs <- c("Matrix", "igraph")
suppressPackageStartupMessages(suppressWarnings(lapply(libs, library, character.only = TRUE)))
### CompS: Computes local entropy of a node with stochastic vector p.v
### correspondig to lines 7-12 from CompSR.R (Teschendorff)
#' @export
CompS <- function(p.v){
  tmp.idx <- which(p.v>0);
  S <-  - sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )
  return(S);
}

### CompNS: 
### compute the normalised local entropy of a node with stochastic vector p.v
### corresponding to lines 14-24 from CompSR.R (Teschendorff)
#' @export
CompNS <- function(p.v){
  tmp.idx <- which(p.v>0);
  if(length(tmp.idx)>1){
    S <-  - sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )/log(length(tmp.idx));
  }else {
    S <- 0;
  }
  return(S);
}


### CompMaxSR: Computes the maximum entropy rate of a network with adjacency matrix adj.m (assumed connected).
#' @export
CompMaxSR <- function(adj.m){
  # find right eigenvector of adjacency matrix
  fa <- function(x,extra=NULL) {
    as.vector(adj.m %*% x)
  }
  ap_o <- arpack(fa, options=list(n=nrow(adj.m),nev=1,which="LM",maxiter=10000),sym=TRUE);
  v <- ap_o$vectors;
  lambda <- ap_o$values;
  maxSR <- log(lambda);
  return(maxSR);
}


### CompSR.ne
### compute the non-equilibrium entropies for the full networks
### corresponding to line 117 from CompSR.R (Teschendorff)
### INPUT: p_m -> proability matrix pre-calculated and saved for one sample
### OUTPUT: list containing the global entropy (as average of the normalised locals) and the normalised locals
#' @export
CompSR.ne <- function(p_m){
  NS.v <- apply(p_m,1,CompNS)
  return(list(sr=mean(NS.v),ns=NS.v))
}


### OUTPUT: a list of four objects
### sr: the entropy rate of the sample (normalised between 0 and 1 if maxSR was provided).
### inv: the stationary distribution of the sample.
### s: the unnormlised local entropies of the sample.
### ns: the normalised local entropies of the sample.
#' @export
CompSR <- function(p_m,maxSR=NULL){
  
  fp <- function(x,extra=NULL) {
    as.vector(p_m %*% x)
  }
  fpt <- function(x,extra=NULL) {
    as.vector(t(p_m) %*% x)
  }
  ap_o <- arpack(fpt, options=list(n=nrow(p_m),nev=1,which="LM",maxiter=10000),sym=FALSE);
  invP_v <- abs(as.numeric(ap_o$vectors));
  invP_v <- invP_v/sum(invP_v);
  S.v <- apply(p_m,1,CompS);
  SR <- sum(invP_v*S.v);
  if(is.null(maxSR)==FALSE){## if provided then normalise relative to maxSR
    SR <- SR/maxSR;
  }
  
  return(list(sr=SR,s=S.v,inv=invP_v));
}

#' Compute Signaling Entropy Rate (SR)
#'
#' @description
#' This function calculates the Signaling Entropy Rate (SR) of a network given its 
#' stochastic transition matrix. SR serves as a core metric for quantifying the 
#' complexity of the network's state space and the efficiency of signal transduction.
#'
#' @details
#' \strong{Mathematical Procedure:}
#' \enumerate{
#'   \item \strong{Stationary Distribution (\eqn{\pi}):} The function employs the ARPACK 
#'         algorithm to solve for the left principal eigenvector of the transition 
#'         matrix \eqn{P}, satisfying the equilibrium condition \eqn{\pi P = \pi}.
#'   \item \strong{Local Node Entropy (\eqn{S_i}):} Computes the Shannon entropy for 
#'         each row of the matrix: \eqn{S_i = -\sum_j p_{ij} \log(p_{ij})}. 
#'         This represents the local uncertainty of signaling outward from node \eqn{i}.
#'   \item \strong{Global Signaling Entropy (SR):} Calculates the weighted average of 
#'         local entropies across the steady-state distribution: \eqn{SR = \sum_i \pi_i S_i}.
#' }
#'
#' @param p_m A numeric stochastic transition matrix (row-normalized). 
#' Supports standard \code{matrix} or sparse \code{dgCMatrix} formats.
#' @param max_sr Optional numeric value representing the maximum theoretical entropy 
#' for normalization. If provided, the function returns \eqn{SR / max\_sr}.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \strong{sr}: The calculated Signaling Entropy Rate (scalar).
#'   \item \strong{s}: A numeric vector of local entropies for each node.
#'   \item \strong{inv}: The stationary distribution vector (\eqn{\pi}), also 
#'         known as the invariant measure of the network.
#' }
#'
#'
#'
#' @importFrom igraph arpack
#' @importFrom Matrix t rowSums
#' @export
comp_sr_optimized <- function(p_m, max_sr = NULL) {
  if (!inherits(p_m, "Matrix") && !is.matrix(p_m)) {
    p_m <- methods::as(p_m, "dgCMatrix")
  }

  pm_t <- Matrix::t(p_m)
  n <- nrow(pm_t)

  ap_o <- igraph::arpack(
    function(x, extra = NULL) as.numeric(pm_t %*% x),
    options = list(
      n = n,
      nev = 1,
      ncv = max(10, min(n, 20)),
      which = "LM",
      maxiter = 20000,
      tol = 1e-10
    ),
    sym = FALSE
  )

  inv_v <- abs(as.numeric(ap_o$vectors))
  inv_v <- inv_v / sum(inv_v)

  if (inherits(p_m, "sparseMatrix")) {
    log_pm <- p_m
    log_pm@x <- log(log_pm@x)
    s_v <- -Matrix::rowSums(p_m * log_pm, na.rm = TRUE)
  } else {
    log_pm <- p_m
    idx <- p_m > 0
    log_pm[idx] <- log(p_m[idx])
    log_pm[!idx] <- 0
    s_v <- -Matrix::rowSums(p_m * log_pm, na.rm = TRUE)
  }

  sr_val <- sum(inv_v * s_v)

  if (!is.null(max_sr)) {
    sr_val <- sr_val / max_sr
  }

  list(
    sr = sr_val,
    s = s_v,
    inv = inv_v
  )
}