#' Accessor for the "counts" slot of Enone object. 
#'
#' @param object Enone. 
#' @param slot Which slot to get, one of \code{sample} or \code{spike_in}.  
#' @param method Which counts matrix to get, must be one of the raw or normalized counts matrix presented in the selected slot. 
#' @param value Raw or normalized counts matrix. 
#' @name Counts
#' @aliases Counts Counts,Enone,character,character-method 
#' Counts<-,Enone,character,character,matrix-method
#' 
#' @return matrix. 
#' @export
#'
setMethod("Counts", signature = signature(object="Enone", slot="character", method="character"), 
          function(object, slot=c("sample","spike_in"), method) {
            
            slot <- match.arg(slot, choices = c("sample","spike_in"))
            
            if (is.null(names(object@counts[[slot]]))) {
              stop("Normalizations for ", slot, " not found. At least one normalization should be performed.")
            }
            method <- match.arg(method, choices = names(object@counts[[slot]]))
            
            object@counts[[slot]][[method]]
            })

#' @rdname Counts
#' @name Counts
#' @export "Counts<-"
setReplaceMethod("Counts", signature = signature(object="Enone", slot="character", method="character", value="matrix"),
                 function(object, slot=c("sample","spike_in"), method, value) {
                   slot <- match.arg(slot, choices = c("sample","spike_in"))
                   object@counts[[slot]][[method]] <- value
                   methods::validObject(object)
                   return(object)
                 })

#' Accessor of normalization methods
#'
#' List all normalization 
#'
#' @param object Enone. 
#' @name listNormalization
#' @aliases listNormalization listNormalization,Enone-method
#' 
#' @return Vector of normalization methods
#' @export
#' 
setMethod("listNormalization", signature = signature(object="Enone"),
          function(object) {
            names(object@counts$sample)
          })

#' Accessor of Enone normalization factors
#'
#' @param object Enone. 
#' @param slot Which slot to get, one of \code{sample} or \code{spike_in}.  
#' @param method Which normalization methods to get, must be one of the methods presented in the selected slot. 
#' @name getFactor
#' @aliases getFactor getFactor,Enone,character,character-method
#'
#' @return vector or list of factors. 
#' @export
#'
setMethod("getFactor", signature = signature(object="Enone", slot="character", method="character"),
          function(object, slot=c("sample","spike_in"), method) {
            slot <- match.arg(slot, choices = c("sample","spike_in"))
            object@enone_factor[[slot]][[method]]
          })

#' Accessor of Enone metrics
#'
#' @param object Enone. 
#' @name getMetrics
#' @aliases getMetrics getMetrics,Enone-method
#'
#' @return data.frame
#' @export
#'
setMethod("getMetrics", signature = signature(object="Enone"),
          function(object) {
            object@enone_metrics
          })

#' Accessor of Enone score
#'
#' @param object Enone. 
#' @name getScore
#' @aliases getScore getScore,Enone-method
#' 
#' @return data.frame
#' @export
#'
setMethod("getScore", signature = signature(object="Enone"),
          function(object) {
            object@enone_score
          })

#' Accessor of Enone parameter
#'
#' List all parameters. 
#'
#' @param object Enone. 
#' @name getParameter
#' @aliases listParameter listParameter,Enone-method
#' 
#' @return Vector of parameter names
#' @export
setMethod("listParameter", signature = signature(object="Enone"),
          function(object) {
            names(object@parameter)
          })


#' Accessor of Enone parameter
#'
#' Get specific parameter.
#'  
#' @param object Enone. 
#' @param name Name of the parameter. 
#' @name getParameter
#' @aliases getParameter getParameter,Enone,character-method
#' 
#' @return Vector of parameter 
#' @export
#'
setMethod("getParameter", signature = signature(object="Enone", name="character"),
          function(object, name) {
            object@parameter[[name]]
          })

#' Accessor of gene set
#'
#' Get gene set. 
#'  
#' @param object Enone. 
#' @param name Name of the gene set, one of \code{NegControl}, 
#' \code{NegEvaluation}, or \code{PosEvaluation}
#' @name getGeneSet
#' @aliases getGeneSet getGeneSet,Enone,character-method
#' 
#' @return Vector of gene ID
#' @export
#'
setMethod("getGeneSet", signature = signature(object="Enone", name="character"),
          function(object, name) {
            
            if (!name %in% c("NegControl","NegEvaluation","PosEvaluation")) {
              stop(name, "not presented in data.")
            }
            
            idx <- SummarizedExperiment::rowData(object)[, name] == TRUE
            rownames(SummarizedExperiment::rowData(object)[idx,])
          })

#' Accessor of enrichment results
#'
#' Get all or filtered enrichment
#'  
#' @param object Enone. 
#' @param slot Which slot to get, one of \code{sample} or \code{spike_in}.  
#' @param filter Whether to get the filtered enrichment, default FALSE. 
#' @name getEnrichment
#' @aliases getEnrichment getEnrichment,Enone,character-method
#' 
#' @return list of enrichment table 
#' @export
#'
setMethod("getEnrichment", signature = signature(object="Enone", slot="character"),
          function(object, slot=c("sample","spike_in"), filter=FALSE) {
            if (filter) {
              object@enrichment_filtered[[slot]]
            } else {
              object@enrichment[[slot]]
            }
          })
