#' @rdname createEnone
#' @export
#'
setClass(
  Class = "Enone",
  contains = "SummarizedExperiment",
  slots = list(
    counts = "list",
    enone_factor = "list",
    enone_metrics = "data.frame",
    enone_score = "data.frame",
    enrichment = "list",
    enrichment_filtered = "list",
    parameter = "list"
  )
)

#' Enone object and constructor
#' 
#' @description \code{Enone} object extends the \code{SummarizedExperiment} class. 
#' The \code{createEnone} is a easy constructor of \code{Enone} object
#' 
#' @param data A un-normalized count data matrix of shape n x p, where n is the 
#'   number of samples and p is the number of features. 
#' @param bio.group Vector of samples group, e.g., c("Young.Input","Young.Enrich","Old.Input","Old.Enrich").
#' @param enrich.group Vector of enrichment group, e.g., c("Input","Enrich","Input","Enrich").
#' @param batch.group Vector of samples batch, e.g., c("A","A","B","B"), default: NULL. 
#' @param spike.in.prefix A character specify the prefix of spike-in id, e.g., "^FB" stands for fly spike-in id, default: NULL. 
#' @param input.id Input library id, must be consistent with \code{enrich.group}, e.g., "Input". 
#' @param enrich.id Enrich library id, must be consistent with \code{enrich.group}, e.g., "Enrich".  
#' @param synthetic.id Vector of synthetic RNA id, e.g. c("Syn1","Syn2"), default: NULL. 
#' 
#' @return Enone object
#' 
#' @details Description of each slot:
#'   \code{assay} \code{SummarizedExperiment::Assays} object, contains all counts. 
#'   \code{counts} list for holding raw/normalized counts of sample and spike_in. 
#'   \code{enone_factor} list for holding normalization factors of sample and spike_in. 
#'   \code{enone_metrics} data frame with normalization methods in row and metrics in columns. 
#'   \code{enone_score} data frame with normalization methods in row and scores in columns.  
#'   \code{enrichment} list for holding all differential analysis results. 
#'   \code{enrichment_filtered} lists for holding filtered differential analysis results. 
#'   \code{parameter} list of parameters.
#'   
#' @export
#' 
#' @importFrom methods validObject
#'
createEnone <- function(data, 
                        bio.group, 
                        enrich.group, 
                        batch.group = NULL,
                        spike.in.prefix = NULL, 
                        input.id = "Input", 
                        enrich.id = "Enrich",
                        synthetic.id = NULL
                        ) {
  
  # parameters 
  params <- list(
    spike.in.prefix = spike.in.prefix,
    synthetic.id = synthetic.id,
    input.id = input.id,
    enrich.id = enrich.id
  )
  
  # create Assay object
  ## count_assay is a list for holding raw/normalized counts of sample and spike_in
  count_assay <- list(sample = list(), spike_in = list())
  ## enrichment_assay and enrichment_filtered_assay are lists for holding all/filtered differential analysis results
  enrichment_assay <- enrichment_filtered_assay <- list(sample=list(), spike_in=list())
  
  # factor slot
  enone_factor <- list(sample=list(), spike_in=list())
  
  # SummarizedExperiment object
  ## rowData for mapping gene id 
  rowDf <- S4Vectors::DataFrame(GeneID=rownames(data),
                                SpikeIn=(rownames(data) %in% grep(spike.in.prefix,rownames(data),value=TRUE)) 
  )
  # if synthetic RNA id provided
  if (!is.null(synthetic.id)) {
    rowDf$Synthetic <- rowDf$GeneID %in% grep(paste(synthetic.id,collapse = "|"),rowDf$GeneID,value=TRUE)
  } else {
    rowDf$Synthetic <- rep(FALSE, nrow(rowDf))
  }
  ## colData for mapping samples
  colDf <- S4Vectors::DataFrame(
    id = colnames(data),
    condition = bio.group,
    enrich = enrich.group,
    replicate = countReplicate(bio.group),
    batch = NA_character_
  )
  
  if (!is.null(batch.group)) colDf$batch <- batch.group
  # create SummarizedExperiment object
  se <- SummarizedExperiment::SummarizedExperiment(S4Vectors::SimpleList(as.matrix(data)),
                                                   colData = colDf,
                                                   rowData = rowDf)
  # create Enone object
  Enone <- methods::new("Enone", 
                       se,
                       counts = count_assay,
                       enone_factor = enone_factor,
                       enone_metrics = data.frame(),
                       enone_score = data.frame(),
                       enrichment = enrichment_assay,
                       enrichment_filtered = enrichment_filtered_assay,
                       parameter = params
                       )
  validObject(Enone)
  return(Enone)
}

setValidity("Enone", function(object) {
  
  # check input.id and enrich.id match the enrich.group
  if (!all(c(object@parameter$input.id, object@parameter$enrich.id) %in% object$enrich)) {
    return("The `input.id` and/or `enrich.id` must be the same as `enrich.group`.")
  }
  
  # check spike.in.prefix match the rownames
  if (length(grep(object@parameter$spike.in.prefix, rownames(object))) < 1) {
    return("The `spike.in.prefix` does not match the rownames in the count matrix.")
  }
  
  # check synthetic.id match the rownames
  if (!is.null(object@parameter$synthetic.id) & !all(object@parameter$synthetic.id %in% rownames(object))) {
    return("The `synthetic.id` are not presented in the rownames of count matrix.")
  }
  
  # check length of bio.group match the sample size
  if (length(object$condition) != ncol(object)) {
    return("The number of elements in `bio.group` does not match the sample size.")
  }
  
  # check length of enrich.group match the sample size
  if (length(object$enrich) != ncol(object)) {
    return("The number of elements in enrich.group` does not match the sample size.")
  }
  
  # check length of batch.group match the sample size
  if (!all(is.na(object$batch)) & length(object$batch) != ncol(object)) {
    return("The number of elements in `batch.group` does not match the sample size.")
  }
  
})
