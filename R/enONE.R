#' Normalization and assessment in one function
#'  
#' @param object Enone object.
#' @param auto Whether to automatically select negative control, positive 
#'   evaluation and negative evaluation genes, default: TRUE. 
#' @param n.neg.control Number of negative control genes for RUV normalization, default: 1000. 
#' @param n.pos.eval Number of positive evaluation genes for wanted variation assessment, default: 500.
#' @param n.neg.eval Number of negative evaluation genes for unwanted variation assessment, default: 500.
#' @param neg.control Vector of negative control genes' id for RUV normalization, default: NULL. 
#' @param pos.eval Vector of positive evaluation genes' id for wanted variation assessment, default: NULL.
#' @param neg.eval Vector of negative evaluation genes' id for unwanted variation assessment, default: NULL.
#' @param scaling.method Vector of normalization methods that are applied to the data.
#'   Available methods are: \code{c("TC", "UQ", "TMM", "DESeq", "PossionSeq")}. 
#'   Select one or multiple methods. By default all normalization methods will be applied.
#' @param ruv.norm Whether to perform RUV normalization. 
#' @param ruv.k The number of factors of unwanted variation to be estimated from the data.
#' @param ruv.drop The number of singular values to drop in the estimation of 
#'   unwanted variation, default: 0.  
#' @param pam.krange Integer or vector of integers indicates the number of 
#'   clusters for PAM clustering, default: 2:6. 
#' @param pc.k Integer indicates the metrics will be calculated in the first kth PCs, default: 3.
#'
#' @return Enone object.
#' @export
#'
#' @importFrom utils head
#' @importFrom stringr str_extract
#' @importFrom stats as.formula model.matrix setNames
enONE <- function(object,
                  auto = TRUE, 
                  n.neg.control = 1000, n.pos.eval = 500, n.neg.eval = 500,
                  neg.control = NULL, pos.eval = NULL, neg.eval = NULL,
                  scaling.method = c("TC", "UQ", "TMM", "DESeq", "PossionSeq"),
                  ruv.norm = TRUE, ruv.k = 1, ruv.drop = 0,
                  pam.krange = 2:6, pc.k = 3) {
  
  # retrieve parameters from Enone object
  bio.group <- object$condition
  enrich.group <- object$enrich
  
  if (!any(is.na(object$batch))) {
    batch.group <- object$batch
  } else {
    batch.group <- NULL
  }
  
  spike.in.prefix <- object@parameter$spike.in.prefix
  input.id <- object@parameter$input.id
  enrich.id <- object@parameter$enrich.id
  synthetic.id <- object@parameter$synthetic.id
  
  # create group matrix
  sc_mat <-  CreateGroupMatrix(bio.group)
  enrich_mat <- CreateGroupMatrix(enrich.group)
  
  # get counts
  data <- SummarizedExperiment::assay(object)
  counts_nsp <- data[grep(spike.in.prefix, rownames(data), invert = TRUE),]
  counts_sp <- data[grep(spike.in.prefix, rownames(data)),]
  
  ## gene selection 
  if (auto) {
    ### 1. negative control genes for RUV
    cat(paste("The number of negative control genes for RUV:",n.neg.control,"\n"))
    designMat <- model.matrix(~0+enrich.group)
    deg.en <- edgeRDE(counts_sp,
                      group = enrich.group,
                      design.formula = as.formula("~0+condition"),
                      contrast.df = data.frame(Group1=enrich.id, Group2=input.id)
    )
    # top 1000 (default) non-sig de 
    res_tab <- deg.en$res.ls[[paste(enrich.id, input.id, sep="_")]]
    # res_tab <- subset(res_tab, FDR > 0.05)
    neg.control.set <- head(res_tab[order(res_tab$FDR, decreasing = TRUE),]$GeneID, n=n.neg.control)
    
    ### 2. positive evaluation genes (default 500)
    # if provided, preclude synthetic RNA from evaluation set 
    cat(paste("The number of positive evaluation genes:",n.pos.eval,"\n"))
    deg.en <- edgeRDE(counts_nsp[!rownames(counts_nsp) %in% synthetic.id,],
                      group = enrich.group,
                      design.formula = as.formula("~0+condition"),
                      contrast.df = data.frame(Group1=enrich.id, Group2=input.id)
    )
    res_tab <- deg.en$res.ls[[paste(enrich.id, input.id, sep="_")]]
    pos.eval.set <- head(res_tab[order(res_tab$FDR),]$GeneID, n=n.pos.eval)
    
    ### 3. negative evaluation genes (default 500)
    # if provided, preclude synthetic RNA from evaluation set 
    cat(paste("The number of negative evaluation genes:",n.neg.eval,"\n"))
    de.all <- edgeRDE(counts_nsp[!rownames(counts_nsp) %in% synthetic.id,],
                      group = bio.group,
                      design.formula = as.formula("~condition"),
                      coef = 2:length(unique(bio.group))
    )
    res_tab <- de.all$res.ls[[1]]
    neg.eval.set <- head(res_tab[order(res_tab$FDR, decreasing = TRUE),]$GeneID, n = n.neg.eval)
    
  } else {
    
    if (!all(neg.control %in% rownames(data))) {
      stop("`neg.control` are not presented in the rownames of count matrix.")
    } else {
      neg.control.set <- neg.control  
    }
    
    if (!is.null(pos.eval) & !all(pos.eval %in% rownames(data))) {
      stop("`pos.eval` are not presented in the rownames of count matrix.")
    } else {
      pos.eval.set <- pos.eval  
    }
    
    if (!is.null(neg.eval) & !all(neg.eval %in% rownames(data))) {
      stop("`neg.eval` are not presented in the rownames of count matrix.")
    } else {
      neg.eval.set <- neg.eval  
    }
    
  }
  
  # save gene set to object
  SummarizedExperiment::rowData(object)$NegControl <- SummarizedExperiment::rowData(object)$GeneID %in% neg.control.set
  SummarizedExperiment::rowData(object)$NegEvaluation <- SummarizedExperiment::rowData(object)$GeneID %in% neg.eval.set
  SummarizedExperiment::rowData(object)$PosEvaluation <- SummarizedExperiment::rowData(object)$GeneID %in% pos.eval.set
  
  ## apply normalization 
  cat("Apply normalization...\n")
  norm.ls <- ApplyNormalization(data,
                                scaling.method = scaling.method, 
                                ruv.norm = ruv.norm, ruv.k = ruv.k, ruv.drop = ruv.drop,
                                spike.in.prefix = spike.in.prefix,
                                # below parameters are generated inside function
                                control.idx = neg.control.set, 
                                sc.idx = sc_mat, 
                                enrich.idx = enrich_mat)
  
  ## assessment 
  bio_group_index <- as.numeric(factor(bio.group, levels=unique(bio.group)))
  assay_group_index <- as.numeric(factor(enrich.group, levels=unique(enrich.group)))
  if (!is.null(batch.group)) {
    batch_group_index <- as.numeric(factor(batch.group, levels=unique(batch.group)))
  } else {
    batch_group_index <- NULL
  }
  cat("Perform assessment...\n")
  norm.eval <- AssessNormalization(norm.ls,
                                   pam.krange = pam.krange,
                                   pc.k = pc.k,
                                   batch.group = batch_group_index,
                                   # below parameters are created inside function
                                   bio.group = bio_group_index, 
                                   assay.group = assay_group_index, 
                                   pos.eval.set = pos.eval.set,
                                   neg.eval.set = neg.eval.set)
  
  ## save metrics to Enone object
  object@enone_metrics <- norm.eval$metrics
  ## save score to Enone object
  object@enone_score <- norm.eval$score
  ## add run parameter in Enone object
  parameter.run <- list(
    n.neg.control = n.neg.control,
    n.pos.eval = n.pos.eval,
    n.neg.eval = n.neg.eval,
    scaling.method = scaling.method,
    ruv.norm = ruv.norm,
    ruv.k = ruv.k,
    ruv.drop = ruv.drop,
    pam.krange = pam.krange,
    pc.k = pc.k
  )
  object@parameter <- c(object@parameter, parameter.run)
  
  # only store normalization method names in object for reducing memory cost
  norm.methods <- names(norm.ls)
  object@counts$sample <- stats::setNames(vector("list", length(norm.methods)), nm = norm.methods)
  object@enone_factor$sample <- stats::setNames(vector("list", length(norm.methods)), nm = norm.methods)
  # except 'Raw' matrix
  Counts(object, slot = "sample", method = "Raw") <- counts_nsp
  Counts(object, slot = "spike_in", method = "Raw") <- counts_sp
  
  return(object)
}