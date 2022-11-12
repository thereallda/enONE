###################################################
#		Esimate sequencing depth
###################################################

#' Esimate sequencing depth
#'
#' @description Functions from PossionSeq (https://doi.org/10.1093/biostatistics/kxr031)
#' 
#' @param n The data matrix. The rows are counts for a gene, and the columns are counts from an experiment.
#' @param iter Number of iterations used. Default value: 5. The default value is usually a good choice.
#' @param ct.sum if the total number of reads of a gene across all experiments <= ct.sum, this gene will not be considered for estimating sequencing depth. Default value: 5.
#' @param ct.mean if the mean number of reads of a gene across all experiments <= ct.mean, this gene will not be considered for estimating sequencing depth. Default value: 0.5.
#'
#' @return estimated sequencing depth. a vector. their product is 1.
#' @export
#' 
#' @examples
#' mat <- matrix(rnbinom(300, mu=100, size=1), ncol=5)
#' PS.Est.Depth(mat)
#' 
PS.Est.Depth <- function(n, iter=5, ct.sum=5, ct.mean=0.5) {
  n <- .PS.Filter(dat=list(n=n), ct.sum=ct.sum, ct.mean=ct.mean)$n
  
  seq.depth <- .Est.Depth(n=n, iter=iter)$cmeans
  
  seq.depth <- exp(log(seq.depth) - mean(log(seq.depth)))
  
  return(seq.depth)
}

###################################################
#		Esimate Cmeans using PossionSeq methods
###################################################
#' @importFrom stats quantile
.Est.Depth <- function(n, iter=5) {
  SMALL.VAL <- 1e-8
  cmeans <- colSums(n) / sum(n)
  keep <- NULL
  
  for (i in 1 : iter)
  {
    n0 <- rowSums(n) %*% t(cmeans)
    prop <- rowSums((n - n0) ^ 2 / (n0 + SMALL.VAL))
    qs <- stats::quantile(prop, c(0.25, 0.75))
    keep <- (prop >= qs[1]) & (prop <= qs[2])
    
    cmeans <- colMeans(n[keep, ])
    cmeans <- cmeans / sum(cmeans)
  }
  
  return(list(cmeans=cmeans, keep=keep))
}


############################################################
#		Filter genes with too small counts
############################################################
.PS.Filter <- function(dat, ct.sum=5, ct.mean=0.5){
  if (is.null(dat$gname))
  {
    dat$gname <- 1 : nrow(dat$n)
  }
  
  keep <- (rowMeans(dat$n) > ct.mean) & (rowSums(dat$n) > ct.sum)
  # cat(length(keep) - sum(keep), "genes has been filtered because they contains too small number of reads across the experiments.")
  
  dat$n <- dat$n[keep, ]
  dat$gname <- dat$gname[keep]
  
  return(dat)
}