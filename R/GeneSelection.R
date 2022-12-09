#' Gene set selection 
#'
#' @description Select negative control genes for RUV normalization, positive 
#' and negative evaluation genes for assessment. 
#' 
#' @param object Enone object
#' @param n.neg.control Number of negative control genes for RUV normalization, default: 1000. 
#' @param n.pos.eval Number of positive evaluation genes for wanted variation assessment, default: 500.
#' @param n.neg.eval Number of negative evaluation genes for unwanted variation assessment, default: 500.
#'
#' @return list of genes
#' @export
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom stats model.matrix as.formula
GeneSelection <- function(object,
                          n.neg.control,
                          n.pos.eval,
                          n.neg.eval) {
  
  # retrieve parameters from Enone object
  bio.group <- object$condition
  enrich.group <- object$enrich
  spike.in.prefix <- object@parameter$spike.in.prefix
  input.id <- object@parameter$input.id
  enrich.id <- object@parameter$enrich.id
  synthetic.id <- object@parameter$synthetic.id
  
  # get counts
  data <- SummarizedExperiment::assay(object)
  counts_nsp <- data[grep(spike.in.prefix, rownames(data), invert = TRUE),]
  counts_sp <- data[grep(spike.in.prefix, rownames(data)),]
  
  cat("Gene set selection for normalization and assessment...\n")
  ### 1. negative control genes for RUV
  cat(paste("- The number of negative control genes for normalization:",n.neg.control,"\n"))
  designMat <- stats::model.matrix(~0+enrich.group)
  deg.en <- edgeRDE(counts_sp,
                    group = enrich.group,
                    design.formula = stats::as.formula("~0+condition"),
                    contrast.df = data.frame(Group1=enrich.id, Group2=input.id)
  )
  # top 1000 (default) non-sig de 
  res_tab <- deg.en$res.ls[[paste(enrich.id, input.id, sep="_")]]
  # res_tab <- subset(res_tab, FDR > 0.05)
  neg.control.set <- head(res_tab[order(res_tab$FDR, decreasing = TRUE),]$GeneID, n=n.neg.control)
  
  ### 2. positive evaluation genes (default 500)
  # if provided, preclude synthetic RNA from evaluation set 
  cat(paste("- The number of positive evaluation genes:",n.pos.eval,"\n"))
  deg.en <- edgeRDE(counts_nsp[!rownames(counts_nsp) %in% synthetic.id,],
                    group = enrich.group,
                    design.formula = stats::as.formula("~0+condition"),
                    contrast.df = data.frame(Group1=enrich.id, Group2=input.id)
  )
  res_tab <- deg.en$res.ls[[paste(enrich.id, input.id, sep="_")]]
  pos.eval.set <- head(res_tab[order(res_tab$FDR),]$GeneID, n=n.pos.eval)
  
  ### 3. negative evaluation genes (default 500)
  # if provided, preclude synthetic RNA from evaluation set 
  cat(paste("- The number of negative evaluation genes:",n.neg.eval,"\n"))
  de.all <- edgeRDE(counts_nsp[!rownames(counts_nsp) %in% synthetic.id,],
                    group = bio.group,
                    design.formula = stats::as.formula("~condition"),
                    coef = 2:length(unique(bio.group))
  )
  res_tab <- de.all$res.ls[[1]]
  neg.eval.set <- head(res_tab[order(res_tab$FDR, decreasing = TRUE),]$GeneID, n = n.neg.eval)
  
  return(list("NegControl" = neg.control.set, 
              "PosEvaluation" = pos.eval.set, 
              "NegEvaluation" = neg.eval.set))
}