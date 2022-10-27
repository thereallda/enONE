#' Applies normalization on sequencing data
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is the number of samples and p is the number of features. 
#' @param scaling.method Vector of normalization methods that are applied to the data.
#'   Available methods are: \code{c("TC", "UQ", "TMM", "DESeq", "PossionSeq")}. 
#'   Select one or multiple methods. By default all normalization methods will be applied.
#' @param ruv.norm Whether to perform RUV normalization. 
#' @param ruv.k The number of factors of unwanted variation to be estimated from the data.
#' @param ruv.drop The number of singular values to drop in the estimation of 
#'   unwanted variation, default drop the first singular value that represent the 
#'   difference between enrichment and input. 
#' @param control.idx Vector of the negative control genes for RUV normalization. 
#' @param sc.idx A numeric matrix specifying the replicate samples for which to 
#'   compute the count differences used to estimate the factors of unwanted variation.
#' @param enrich.idx Matrix with two rows indicating the column index of 
#'   enrichment and input samples in the raw/normalized count data matrix. 
#'   The first row is the column index of input and the second row is the 
#'   column index of enrichment samples.
#' @param spike.in.prefix A character specify the prefix of spike-in id. 
#'
#' @return List of objects containing normalized data and associated normalization factors. 
#' @export
ApplyNormalization <- function(data, 
                               scaling.method = c("TC", "UQ", "TMM", "DESeq", "PossionSeq"),
                               ruv.norm = TRUE, 
                               ruv.k = 1, 
                               ruv.drop = 0, 
                               control.idx = NULL,
                               sc.idx = NULL,
                               enrich.idx = NULL,
                               spike.in.prefix = NULL) {
  
  scaling.method <- match.arg(scaling.method,
                      choices = c("TC", "UQ", "TMM", "DESeq", "PossionSeq"),
                      several.ok = TRUE)
  
  # scaling
  if (is.null(spike.in.prefix)) {
    data.scaled <- lapply(scaling.method, function(i) {
      norm.f <- get(paste0("norm", i)) 
      data.norm <- norm.f(data)
    })
    names(data.scaled) <- scaling.method
    data.norm <- data.scaled
    data.raw <- list(dataNorm = data, normFactor = rep(1, ncol(data)))
    data.norm[["Raw"]] <- data.raw
  } else {
    data.spike.in <- data[grep(spike.in.prefix, rownames(data)),]
    data.non.spike.in <- data[grep(spike.in.prefix, rownames(data), invert = TRUE),]
    
    data.scaled <- lapply(scaling.method, function(i) {
      
      norm.f <- get(paste0("norm", i)) 
      data.norm.spike.in <- norm.f(data.spike.in)
      data.norm.non.spike.in <- norm.f(data.non.spike.in)
      # return non-spike-in and negative control counts
      dataNorm <- rbind(data.norm.non.spike.in$dataNorm, data.norm.spike.in$dataNorm[control.idx,])
      
      return(list(
        dataNorm = dataNorm,
        normFactor = data.norm.non.spike.in$normFactor
      ))
      
    })
    names(data.scaled) <- scaling.method
    data.norm <- data.scaled
    data.raw <- list(dataNorm = data[c(grep(spike.in.prefix, rownames(data), invert=TRUE,value=TRUE),control.idx),],
                     normFactor = rep(1, ncol(data)))
    data.norm[["Raw"]] <- data.raw
  }

  # RUV normalization
  if (ruv.norm) {
    ruv.ls <- list()
    for (i in names(data.norm)) {
      data.curr <- data.norm[[i]]$dataNorm
      for (k in 1:ruv.k) {
        if (!is.null(control.idx)) {
          ruv.ls[[paste0(i,"_RUVg_k",k)]] <- normRUV(data.curr,
                                                     control.idx = control.idx,
                                                     method = "RUVg",
                                                     k = k, drop = ruv.drop)
          ruv.ls[[paste0(i,"_RUVg_k",k)]]$normFactor <- data.norm[[i]]$normFactor
          
          if (!is.null(sc.idx)) {
            ruv.ls[[paste0(i,"_RUVs_k",k)]] <- normRUV(data.curr,
                                                       control.idx = control.idx,
                                                       sc.idx = sc.idx,
                                                       method = "RUVs",
                                                       k = k, drop = ruv.drop)
            ruv.ls[[paste0(i,"_RUVs_k",k)]]$normFactor <- data.norm[[i]]$normFactor
          }
          
          if (!is.null(enrich.idx)) {
            ruv.ls[[paste0(i,"_RUVse_k",k)]] <- normRUV(data.curr,
                                                        control.idx = control.idx,
                                                        sc.idx = enrich.idx,
                                                        method = "RUVse",
                                                        k = k, drop = ruv.drop)
            ruv.ls[[paste0(i,"_RUVse_k",k)]]$normFactor <- data.norm[[i]]$normFactor
            
          }
        }
      }
    }
    data.norm <- c(data.norm, ruv.ls)

    # return only non-spike-in counts
    data.norm <- lapply(data.norm, function(x) {
      x$dataNorm <- x$dataNorm[!rownames(x$dataNorm) %in% control.idx,]
      return(x)
      })
  }
  return(data.norm)
}

#' Perform total count normalization
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is 
#' the number of samples and p is the number of features. 
#'
#' @return List containing normalized counts and normalized factors for library size. 
#' @export
#'
normTC <- function(data) {
  
  normFactor <- rep(1,ncol(data))
  sizeFactor <- normFactor*colSums(data)/1e6
  dataNorm <- t(t(data)/sizeFactor)
  
  return(list(
    dataNorm = dataNorm,
    normFactor = normFactor
  ))
}

#' Perform upper-quartile normalization
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is 
#' the number of samples and p is the number of features. 
#'
#' @return List containing normalized counts and normalized factors for library size. 
#' @export
#'
#' @importFrom edgeR calcNormFactors
normUQ <- function(data) {
  
  normFactor <- edgeR::calcNormFactors(data, method = "upperquartile")
  sizeFactor <- normFactor*colSums(data)/1e6
  dataNorm <- t(t(data)/sizeFactor)
  
  return(list(
    dataNorm = dataNorm,
    normFactor = normFactor
  ))
}

#' Perform TMM normalization
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is 
#' the number of samples and p is the number of features.
#'
#' @return List containing normalized counts and normalized factors for library size. 
#' @export
#'
#' @importFrom edgeR calcNormFactors
normTMM <- function(data) {
  
  normFactor <- edgeR::calcNormFactors(data, method = "TMM")
  sizeFactor <- normFactor*colSums(data)/1e6
  dataNorm <- t(t(data)/sizeFactor)
  
  return(list(
    dataNorm = dataNorm,
    normFactor = normFactor
  ))
}

#' Perform DESeq normalization
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is 
#' the number of samples and p is the number of features.
#'
#' @return List containing normalized counts and normalized factors for library size. 
#' @export
#'
#' @importFrom edgeR calcNormFactors
normDESeq <- function(data) {
  
  normFactor <- edgeR::calcNormFactors(data, method = "RLE")
  sizeFactor <- normFactor*colSums(data)/1e6
  dataNorm <- t(t(data)/sizeFactor)
  
  return(list(
    dataNorm = dataNorm,
    normFactor = normFactor
  ))
}

#' Perform PossionSeq normalization
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is 
#' the number of samples and p is the number of features.
#' @param ... Additional parameters can be passed to \code{PS.Est.Depth}. 
#'
#' @return List containing normalized counts and normalized factors for library size. 
#' @export
#'
normPossionSeq <- function(data, ...) {
  
  sizeFactor <- PS.Est.Depth(data, ...)
  normFactor <- sizeFactor*1e6/colSums(data)
  dataNorm <- t(t(data)/sizeFactor)
  
  return(list(
    dataNorm = dataNorm,
    normFactor = normFactor
  ))
}

#' Perform RUV normalization 
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is 
#' the number of samples and p is the number of features.
#' @param control.idx Vector of control genes' id. 
#' @param sc.idx A numeric matrix specifying the replicate samples for which to 
#' compute the count differences used to estimate the factors of unwanted variation.
#' @param method Perform RUVg or RUVs normalization. 
#' @param k The number of factors of unwanted variation to be estimated from the data.
#' @param drop The number of singular values to drop in the estimation of 
#' unwanted variation, default drop the first singular value that represent the 
#' difference between enrichment and input. 
#' @param log Whether to perform log2-transformation with 1 offset on data matrix, 
#' default: TRUE. 
#'
#' @return List containing normalized counts and adjust factors for adjusting unwanted variation. 
#' @export
#'
normRUV <- function(data, 
                    control.idx = NULL, 
                    sc.idx = NULL, 
                    method = c("RUVg", "RUVs", "RUVse"),
                    k = 1,
                    drop = 0,
                    log = TRUE) {
  
  if (is.null(control.idx)) {
    control.idx <- rownames(data)
  }
 
  if (log) {
    dataNorm <- log2(data + 1)
  } else {
    dataNorm <- data
  }
  
  if (method == "RUVg") {
    ruv.set <- enRUVg(dataNorm, control.idx=control.idx, k=k, drop=drop, log=FALSE)
    dataNorm <- 2^(ruv.set$normalizedCounts)
  }
  
  if (method %in% c("RUVs","RUVse")) {
    ruv.set <- enRUVs(dataNorm, control.idx=control.idx, k=k, drop=drop, sc.idx=sc.idx, log=FALSE)
    dataNorm <- 2^(ruv.set$normalizedCounts)
  }
  return(list(
    dataNorm = dataNorm,
    adjustFactor = ruv.set$W,
    alpha = ruv.set$alpha
  ))
}
