#' Assessments of normalization performance
#'
#' @param data.ls List containing normalized counts and adjust factors for 
#'   adjusting unwanted variation. Output of \code{ApplyNormalization}. 
#' @param bio.group Vector of index indicating the column index of samples of 
#'   each biological groups in the raw/normalized count data matrix. 
#' @param enrich.group Vector of index indicating the column index of 
#'   enrichment and input samples in the raw/normalized count data matrix. 
#' @param batch.group Vector of index indicating the column index of 
#'   each batch groups in the raw/normalized count data matrix. 
#' @param eval.pam.k Integer or vector of integers indicates the number of 
#'   clusters for PAM clustering in performance evaluation, default: 2:6. 
#' @param eval.pc.n Integer indicates the evaluation metrics will be calculated 
#'   in the first nth PCs, default: 3.
#' @param log Whether to perform log2-transformation with 1 offset on data matrix, 
#'   default: TRUE. 
#' @param pos.eval.set Vector of genes id.
#' @param neg.eval.set Vector of genes id.
#'
#' @return List containing the metrics matrix and the ranking matrix, 
#'   both sorted by the score of methods from top to bottom. 
#' @export
#'
#' @importFrom stats dist cor lm na.omit var
#' @importFrom cluster silhouette 
#' @importFrom MatrixGenerics rowMedians colMedians colIQRs
#' @importFrom fpc pamk
#' @importFrom pbapply pblapply
AssessNormalization <- function(data.ls, 
                                bio.group = NULL, 
                                enrich.group = NULL, 
                                batch.group = NULL,
                                eval.pam.k = 2:6, 
                                eval.pc.n = 3, 
                                log = TRUE, 
                                pos.eval.set = NULL, 
                                neg.eval.set = NULL) {
  
  metrics.ls <- pbapply::pblapply(1:length(data.ls), function(i) {
    data <- as.matrix(data.ls[[i]]$dataNorm)
    # Clustering properties
    if (log) {
      data.log <- log2(data + 1)
    } else {
      data.log <- data
    }
    # remove constant genes
    data.log <- data.log[apply(data.log, 1, var, na.rm=TRUE) !=0, ]
    # PCA on expression matrix
    pca.expr <- prcomp(scale(t(data.log)))
    # compute right singular value by svd
    expr_sv <- svd(scale(t(data.log), center = TRUE, scale = TRUE),
                   nu = eval.pc.n, nv = 0)$u
    # calculate euclidean distance in the space of first k PCs (default: 3)
    dist.pca.expr <- dist(scale(pca.expr$x[, 1:eval.pc.n]), method = "euclidean")
    # dist.pca.expr <- dist(expr_sv, method = "euclidean")
    # silhouette width
    if (length(bio.group) == ncol(data)) {
      bio_sil <- mean(cluster::silhouette(bio.group, dist.pca.expr)[,"sil_width"])
    } else {
      bio_sil <- 0
    }
    if (length(enrich.group) == ncol(data)) {
      en_sil <- mean(cluster::silhouette(enrich.group, dist.pca.expr)[,"sil_width"])
    } else {
      en_sil <- 0
    }
    if (length(batch.group) == ncol(data)) {
      batch_sil <- mean(cluster::silhouette(batch.group, dist.pca.expr)[,"sil_width"])
    } else {
      batch_sil <- 0
    }
    
    prk <- fpc::pamk(pca.expr$x[,1:eval.pc.n], krange=eval.pam.k) # PAM clustering with user specified k
    pam_sil <- prk$pamobject$silinfo$avg.width
    
    # Global distribution properties
    data.log.rle <- data.log - rowMedians(data.log)
    # Mean squared Median RLE
    rle_med <- mean(colMedians(data.log.rle)^2)
    # Variance of IQR of RLE
    rle_iqr <- var(colIQRs(data.log.rle))
    
    # Association with control genes
    # wanted factors from positive set 
    if (!is.null(pos.eval.set)) {
      wv_factors <- svd(scale(t(data.log[rownames(data.log) %in% pos.eval.set,]), center = TRUE, scale = TRUE),
                        nu = eval.pc.n, nv = 0)$u
      # weighted coefficient of determination
      wv_cor <- 1 - sum(unlist(apply(expr_sv, 2, function(y) {
        lm(y ~ wv_factors)$residual
      })) ^ 2) / sum(scale(expr_sv, scale = FALSE) ^ 2)
    } else {
      wv_cor <- 0
    }
    
    # unwanted factors from negative set
    if (!is.null(neg.eval.set)) {
      uv_factors <- svd(scale(t(data.log[rownames(data.log) %in% neg.eval.set,]), center = TRUE, scale = TRUE),
                        nu = eval.pc.n, nv = 0)$u
      # weighted coefficient of determination
      uv_cor <- 1 - sum(unlist(apply(expr_sv, 2, function(y) {
        lm(y ~ uv_factors)$residual
      })) ^ 2) / sum(scale(expr_sv, scale = FALSE) ^ 2)
    } else {
      uv_cor <- 0
    }
    
    metrics <- c(
      BIO_SIM = bio_sil,
      EN_SIM = en_sil,
      BATCH_SIM = batch_sil,
      PAM_SIM = pam_sil,
      RLE_MED = rle_med,
      RLE_IQR = rle_iqr,
      WV_COR = wv_cor,
      UV_COR = uv_cor
    )
  })
  
  # reduce list of metrics into table
  # with methods in row and measures in column
  metrics <- data.frame(do.call(rbind, metrics.ls))
  rownames(metrics) <- names(data.ls)
  
  # multiplying by +/- 1 so that large values correspond to good performance
  score <- t(t(metrics) * c(1,1,-1,1,-1,-1,1,-1))  # BIO_SIM,EN_SIM,BATCH_SIM,PAM_SIM,RLE_MED,RLE_IQR,WV_COR,UV_COR
  # rank score
  ranked_score <- apply(na.omit(score), 2, rank, ties.method = "min")
  # mean score rank
  if (is.null(dim(ranked_score))) {
    mean_score_rank <- ranked_score
  } else {
    # if score all 1, remove it before scoring
    metrics.keep <- colSums(ranked_score==1) != nrow(ranked_score)
    mean_score_rank <- rowMeans(ranked_score[,metrics.keep])
  }
  
  ranked_score <- as.data.frame(ranked_score)
  ranked_score$SCORE <- mean_score_rank
  ranked_score <- ranked_score[order(mean_score_rank, decreasing = TRUE),]
  
  metrics <- metrics[order(mean_score_rank, decreasing = TRUE), ]
    
  return(list(
    metrics = metrics,
    score = ranked_score
  ))
}


