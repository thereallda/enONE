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
#' @param synthetic.id Character or vector of string specifying the name of synthetic RNAs.  
#' 
#' @return List of objects containing normalized data and associated normalization factors. 
#' @export
#' 
#' @importFrom pbapply pblapply
ApplyNormalization <- function(data, 
                               scaling.method = c("TC", "UQ", "TMM", "DESeq", "PossionSeq"),
                               ruv.norm = TRUE, 
                               ruv.k = 1, 
                               ruv.drop = 0, 
                               control.idx = NULL,
                               sc.idx = NULL,
                               enrich.idx = NULL,
                               spike.in.prefix = NULL,
                               synthetic.id = NULL) {
  
  scaling.method <- match.arg(scaling.method,
                      choices = c("TC", "UQ", "TMM", "DESeq", "PossionSeq"),
                      several.ok = TRUE)
  
  # scaling
  if (is.null(spike.in.prefix)) {
    cat("- Scaling... \n")
    data.scale <- pbapply::pblapply(1:length(scaling.method), function(i) {
      norm.f <- get(paste0("norm", scaling.method[i])) 
      data.norm <- norm.f(data)
    })
    names(data.scale) <- scaling.method
    data.norm <- data.scale
    data.raw <- list(dataNorm = data, normFactor = rep(1, ncol(data)))
    data.norm[["Raw"]] <- data.raw
  } else {
    data.spike.in <- data[grep(spike.in.prefix, rownames(data)),]
    data.non.spike.in <- data[grep(paste(c(spike.in.prefix, synthetic.id), collapse = "|"), rownames(data), invert = TRUE),]
    
    cat("- Scaling... \n")
    data.scale <- pbapply::pblapply(1:length(scaling.method), function(i) {
      norm.f <- get(paste0("norm", scaling.method[i])) 
      data.norm.spike.in <- norm.f(data.spike.in)
      data.norm.non.spike.in <- norm.f(data.non.spike.in)
      
      # return non-spike-in and negative control counts
      dataNorm <- rbind(data.norm.non.spike.in$dataNorm, data.norm.spike.in$dataNorm[control.idx,])
      
      return(list(
        dataNorm = dataNorm,
        normFactor = data.norm.non.spike.in$normFactor
      ))
    })
    names(data.scale) <- scaling.method
    data.norm <- data.scale
    data.raw <- list(dataNorm = data[c(grep(spike.in.prefix, rownames(data), invert=TRUE,value=TRUE),control.idx),],
                     normFactor = rep(1, ncol(data)))
    data.norm[["Raw"]] <- data.raw
  }

  # RUV normalization
  if (ruv.norm & !is.null(control.idx)) {
    # stop when provided ruv.k larger than sample size
    if (ruv.k > ncol(data)) {
      stop("Number of `ruv.k` must not exceed the number of samples.")
    }
    
    cat("- Regression-based normalization... \n")
    # generate all integrated methods
    # only perform RUVs, when enrichment is the only covariate of interest
    if (identical(sc.idx, enrich.idx)) {
      integrated.methods <- paste(rep(paste(rep(c("Raw",scaling.method),each=2), c("RUVg","RUVs"), sep="_"), each=ruv.k), paste0("k", 1:ruv.k), sep="_")
    } else {
      integrated.methods <- paste(rep(paste(rep(c("Raw",scaling.method),each=3), c("RUVg","RUVs","RUVse"), sep="_"), each=ruv.k), paste0("k", 1:ruv.k), sep="_")
    }
    
    ruv.ls <- pbapply::pblapply(1:length(integrated.methods), function(i) {
      # method.curr[1]: scaling; method.curr[2]: RUV; method.curr[3]: number of k
      method.curr <- unlist(strsplit(integrated.methods[i], split = "_"))
      
      # get current scaled data and scaling factors
      data.curr <- data.norm[[method.curr[1]]]$dataNorm
      normFactor.curr <- data.norm[[method.curr[1]]]$normFactor
      
      # apply all RUV
      if (method.curr[2] %in% c("RUVg", "RUVs", "RUVse")) {
        # switch sc.idx
        sc.idx <- switch(method.curr[2],
                         "RUVg" = NULL,
                         "RUVs" = sc.idx,
                         "RUVse" = enrich.idx)
        ruv.curr <- normRUV(data.curr,
                            control.idx = control.idx,
                            sc.idx = sc.idx,
                            method = method.curr[2],
                            k = as.numeric(gsub("k", "", method.curr[3])))
        ruv.curr$normFactor <- normFactor.curr
      }
      return(ruv.curr)
    })
    names(ruv.ls) <- integrated.methods
      
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
#' @param log Whether to perform log2-transformation with 1 offset on data 
#' matrix (default: TRUE), while normalized counts are returned in non-log format.  
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
  
  # index for all zero
  zero.idx <- data == 0
  
  if (log) {
    dataNorm <- log2(data + 1)
  } else {
    dataNorm <- data
  }
  
  if (method == "RUVg") {
    ruv.set <- enRUVg(dataNorm, control.idx=control.idx, k=k, drop=drop, log=FALSE)
  }
  
  if (method %in% c("RUVs","RUVse")) {
    ruv.set <- enRUVs(dataNorm, control.idx=control.idx, k=k, drop=drop, sc.idx=sc.idx, log=FALSE)
  }
  
  # return natural values
  dataNorm <- 2^(ruv.set$normalizedCounts)
  # restore zero
  dataNorm[zero.idx] <- 0
  
  return(list(
    dataNorm = dataNorm,
    adjustFactor = ruv.set$W,
    alpha = ruv.set$alpha
  ))
}
