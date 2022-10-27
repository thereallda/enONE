#' @rdname Counts
#' @export
setGeneric("Counts", function(object, slot=c("sample","spike_in"), method) standardGeneric("Counts"))

#' @rdname Counts
#' @export
setGeneric("Counts<-", function(object, slot=c("sample","spike_in"), method, value) standardGeneric("Counts<-"))

#' @rdname listNormalization
#' @export
setGeneric("listNormalization", function(object) standardGeneric("listNormalization"))

#' @rdname getFactor
#' @export
setGeneric("getFactor", function(object, slot=c("sample","spike_in"), method) standardGeneric("getFactor"))

#' @rdname getMetrics
#' @export
setGeneric("getMetrics", function(object) standardGeneric("getMetrics"))

#' @rdname getScore
#' @export
setGeneric("getScore", function(object) standardGeneric("getScore"))

#' @rdname getParameter
#' @export
setGeneric("listParameter", function(object) standardGeneric("listParameter"))

#' @rdname getParameter
#' @export
setGeneric("getParameter", function(object, name) standardGeneric("getParameter"))

#' @rdname getGeneSet
#' @export
setGeneric("getGeneSet", function(object, name) standardGeneric("getGeneSet"))

#' @rdname getEnrichment
#' @export
setGeneric("getEnrichment", function(object, slot=c("sample","spike_in"), filter=FALSE) standardGeneric("getEnrichment"))
