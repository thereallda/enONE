#' Remove unwanted variation using control genes
#' 
#' @description Adapted from \code{RUVSeq} (https://doi.org/10.1038%2Fnbt.2931)
#' 
#' @param object A counts matrix.
#' @param log Whether to perform log2-transformation with 1 offset on data matrix, 
#' default: TRUE. 
#' @param k The number of factors of unwanted variation to be estimated from the data.
#' @param tolerance Tolerance in the selection of the number of positive singular 
#' values, i.e., a singular value must be larger than tolerance to be considered positive.
#' @param control.idx ID of control genes.
#' @param drop The number of singular values to drop in the estimation of 
#' unwanted variation, default drop the first singular value that represent the 
#' difference between enrichment and input
#'
#' @return list contain a matrix of unwanted factors (W) and corrected counts matrix (normalizedCounts).
#' @export
#'
enRUVg <- function(object, log=TRUE, k=1, tolerance=1e-8, control.idx, drop=1) {
  if(log) {
    Y <- t(log2(object+1))
  } else {
    Y <- t(object)
  }
  Ycenter <- apply(Y, 2, function(x) scale(x, center = TRUE, scale = FALSE))
  svdWa <- svd(Ycenter[, control.idx])
  first <- 1+drop
  # d	a vector containing the singular values of x, of length min(n, p), sorted decreasingly.
  k <- min(k, max(which(svdWa$d > tolerance)))
  # u	a matrix whose columns contain the left singular vectors of x, present if nu > 0. Dimension c(n, nu).
  W <- svdWa$u[, (first:k), drop = FALSE]
  colnames(W) <- paste("W", seq(1, ncol(W)), sep="_")
  
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
  correctedY <- Y - W %*% alpha
  return(list(W = W, normalizedCounts = t(correctedY), alpha = alpha))
}

#' Remove unwanted variation using control genes within replicates
#' 
#' @description Adapted from \code{RUVSeq} (https://doi.org/10.1038%2Fnbt.2931)
#'
#' @param object A counts matrix.
#' @param log Whether to perform log2-transformation with 1 offset on data matrix, 
#' default: TRUE. 
#' @param k The number of factors of unwanted variation to be estimated from the data.
#' @param tolerance Tolerance in the selection of the number of positive singular 
#' values, i.e., a singular value must be larger than tolerance to be considered positive.
#' @param control.idx ID of control genes.
#' @param drop The number of singular values to drop in the estimation of 
#' unwanted variation, default drop the first singular value that represent the 
#' difference between enrichment and input
#' @param sc.idx A numeric matrix specifying the replicate samples for which to 
#' compute the count differences used to estimate the factors of unwanted variation 
#' (see details).
#' 
#' @details Each row of sc.idx should correspond to a set of replicate samples. 
#' The number of columns is the size of the largest set of replicates; 
#' rows for smaller sets are padded with -1 values. 
#' 
#' For example, if the sets of replicate samples are (1,2,3),(4,5),(6,7,8), then scIdx should be
#' 1 2 3 
#' 4 5 -1
#' 6 7 8
#' @return list contain a matrix of unwanted factors (W) and corrected counts matrix (normalizedCounts).
#' @export
#'
enRUVs <- function(object, log=TRUE, k=2, tolerance=1e-8, control.idx, sc.idx, drop=1) {
  if(log) {
    Y <- t(log2(object+1))
  } else {
    Y <- t(object)
  }
  sc.idx <- sc.idx[rowSums(sc.idx > 0) >= 2, , drop = FALSE]
  Yctls <- matrix(0, prod(dim(sc.idx)), ncol(Y))
  m <- nrow(Y)
  n <- ncol(Y)
  c <- 0
  for (ii in 1:nrow(sc.idx)) {
    for (jj in 1:(ncol(sc.idx))) {
      if (sc.idx[ii, jj] == -1)
        next
      c <- c + 1
      Yctls[c, ] <- Y[sc.idx[ii, jj], , drop = FALSE] -
        colMeans(Y[sc.idx[ii, (sc.idx[ii, ] > 0)], , drop = FALSE])
    }
  }
  Yctls <- Yctls[rowSums(Yctls) != 0, ]
  Y <- rbind(Y, Yctls)
  sctl <- (m + 1):(m + nrow(Yctls))
  svdRes <- svd(Y[sctl, ], nu = 0, nv = k)
  first <- 1+drop
  k <- min(k, max(which(svdRes$d > tolerance)))
  
  a <- diag(as.vector(svdRes$d[(first:k)]), ncol=k-drop, nrow=k-drop) %*% t(as.matrix(svdRes$v[, (first:k), drop = FALSE]))
  colnames(a) <- colnames(Y)
  # Estimate the unwanted factors W by OLS regression of Z
  W <- Y[, control.idx] %*% t(solve(a[, control.idx, drop = FALSE] %*% t(a[, control.idx, drop = FALSE]), a[, control.idx, drop = FALSE]))
  #Wa <- W %*% a
  correctedY <- Y[1:m, ] - W[1:m, ] %*% a
  
  W <- as.matrix(W[1:m,])
  colnames(W) <- paste("W", seq(1, ncol(W)), sep="_")
  return(list(W = W, normalizedCounts = t(correctedY), alpha = a))
}
