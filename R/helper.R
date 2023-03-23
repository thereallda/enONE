#' Find enriched genes between enrich and input samples
#'
#' @param object An object.
#' @param slot Which slot, one of \code{sample} or \code{spike_in}.  
#' @param method Which normalization methods, must be one of the methods presented in the selected slot.  
#' @param logfc.cutoff Filter genes by at least X-fold difference (log2-scale) 
#' between the two groups of samples, default: 1. 
#' @param p.cutoff Filter genes by no more than \code{p.cutoff} adjusted p-value, default: 0.05. 
#' @param ... Additional parameters can be passed to \code{edgeRDE}. 
#' 
#' @name FindEnrichment 
#' @aliases FindEnrichment  FindEnrichment,Enone,character,character,numeric,numeric-method
#' 
#' @return updated Enone 
#' @export
#'
FindEnrichment <-  function(object, slot=c("sample","spike_in"), method, 
                            logfc.cutoff=1, p.cutoff=0.05, ...) {
  
  contrast_df <- data.frame(Group1 = unique(grep(object@parameter$enrich.id, object$condition, value = TRUE)),
                            Group2 = unique(grep(object@parameter$input.id, object$condition, value = TRUE)))
  # extract sample or spike-in counts
  slot <- match.arg(slot, choices = c("sample","spike_in"))
  
  # test if chosen method in object 
  if (is.null(names(object@enone_factor[[slot]]))) {
    stop("Normalizations for ", slot, " not found. At least one normalization should be performed.")
  }
  method <- match.arg(method, choices = names(object@enone_factor[[slot]]))
  
  if (slot == "spike_in") {
    counts_df <- SummarizedExperiment::assay(object)[SummarizedExperiment::rowData(object)$SpikeIn,]
  } 
  else if (slot == "sample" & any(SummarizedExperiment::rowData(object)$Synthetic)) {
    counts_df <- SummarizedExperiment::assay(object)[!SummarizedExperiment::rowData(object)$SpikeIn & !SummarizedExperiment::rowData(object)$Synthetic,]
  } 
  else {
    counts_df <- SummarizedExperiment::assay(object)[!SummarizedExperiment::rowData(object)$SpikeIn,]
  }
  
  # get list of factors 
  factor.ls <- getFactor(object, slot=slot, method=method)
  if ( !is.null(factor.ls[["normFactor"]]) ) {
    if ( !is.null(factor.ls[["adjustFactor"]]) ) {
      # if norm factors and adjust factors were both provided
      de <- edgeRDE(counts = counts_df,
                    group = object$condition,
                    contrast.df = contrast_df,
                    norm.factors = factor.ls$normFactor,
                    adjust.factors = factor.ls$adjustFactor,
                    logfc.cutoff = logfc.cutoff, p.cutoff = p.cutoff, ...)
    } else {
      # if only norm factors were provided
      de <- edgeRDE(counts = counts_df,
                    group = object$condition,
                    contrast.df = contrast_df,
                    norm.factors = factor.ls$normFactor,
                    design.formula = as.formula("~0+condition"),
                    logfc.cutoff = logfc.cutoff, p.cutoff = p.cutoff, ...)
    }
  } else {
    stop("One or both of \"normFactor\" and \"adjustFacotr\" should be provided.")
  }
  # save enrichment in Enone object
  object@enrichment[[slot]] <- de$res.ls
  object@enrichment_filtered[[slot]] <- de$res.sig.ls
  
  return(object)
}

#' Apply specific normalization method
#'
#' @param object Enone object. 
#' @param slot Which slot, one of \code{sample} or \code{spike_in}.  
#' @param method Which normalization methods to perform. 
#'
#' @return updated Enone
#' @export
#' 
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom stringr str_extract
UseNormalization <- function(object, slot=c("sample","spike_in"), method) {
  
  # sample or spike-in 
  slot <- match.arg(slot, choices = c("sample","spike_in"))
  
  # get raw counts
  if (slot == "spike_in") {
    # only spike-in counts
    counts_df <- SummarizedExperiment::assay(object)[SummarizedExperiment::rowData(object)$SpikeIn,]
  } 
  else {
    # All counts
    counts_df <- SummarizedExperiment::assay(object)
  }
  
  # method.curr[1]: scaling; method.curr[2]: RUV; method.curr[3]: number of k
  method.curr <- unlist(strsplit(method, split = "_"))
  
  # check whether selected method can be provided
  if (!method.curr[1] %in% c("Raw","TC","UQ","DESeq","TMM","PossionSeq")) {
    stop("Scaling method: ", method.curr[1], " is not provided. It should be one of ", c("Raw","TC","UQ","DESeq","TMM","PossionSeq"))
  }
  if (length(method.curr) > 1 & !method.curr[2] %in% c("RUVg","RUVs","RUVse")) {
    stop("RUV method: ", method.curr[2], " is not provided. It should be one of ", c("RUVg","RUVs","RUVse"))
  }
  if (length(method.curr) > 2 & is.na(stringr::str_extract(method.curr[3],"k"))) {
    stop("Use x number of factors to estimate unwanted variation as \"kx\", but not ", method.curr[3])
  }  
  if (length(method.curr) > 2 & as.numeric(stringr::str_extract(method.curr[3],"\\d")) > ncol(object)) {
    stop("Number of required factors exceed, try least.")
  } 
  
  # scaling
  neg.control <- getGeneSet(object, "NegControl")
  if (method.curr[1] == "Raw") {
    counts_scale <- list(dataNorm=counts_df, normFactor=rep(1, ncol(counts_df)))
  } 
  else {
    normScaling <- get(paste0("norm", method.curr[1]))
    
    if (slot == "spike_in") {
      counts_scale <- normScaling(counts_df)
    } 
    else {
      counts_sp_scale <- normScaling(counts_df[SummarizedExperiment::rowData(object)$SpikeIn,])
      counts_nsp_scale <- normScaling(counts_df[!SummarizedExperiment::rowData(object)$SpikeIn,])
      dataNorm <- rbind(counts_nsp_scale$dataNorm, counts_sp_scale$dataNorm[neg.control,])
      # counts_scale contain both normalized data and normalization factors
      counts_scale <- list(dataNorm = dataNorm, 
                           normFactor = counts_nsp_scale$normFactor)
    }
  }
  
  # RUV
  sc_idx <- CreateGroupMatrix(object$condition)
  enrich_idx <- CreateGroupMatrix(object$enrich)
  if (method.curr[2] %in% c("RUVg", "RUVs", "RUVse")) {
    
    sc.idx <- switch(method.curr[2],
                     "RUVg" = NULL,
                     "RUVs" = sc_idx,
                     "RUVse" = enrich_idx)
    
    if (slot == "spike_in") {
      counts_norm <- normRUV(counts_scale$dataNorm,
                             control.idx = neg.control,
                             sc.idx = sc.idx,
                             method = method.curr[2],
                             k = as.numeric(gsub("k", "", method.curr[3])))
    } 
    else {
      counts_norm <- normRUV(counts_scale$dataNorm,
                             control.idx = neg.control,
                             sc.idx = sc.idx,
                             method = method.curr[2],
                             k = as.numeric(gsub("k", "", method.curr[3])))
      # return only non-spike-in counts
      counts_norm$dataNorm <- counts_norm$dataNorm[!rownames(counts_norm$dataNorm) %in% neg.control,]
    }
    
  } 
  else {
    counts_norm <- counts_scale
  }
  
  Counts(object, slot=slot, method=method) <- counts_norm$dataNorm
  
  object@enone_factor[[slot]][[method]] <- list(normFactor=counts_scale$normFactor,
                                                adjustFactor=counts_norm$adjustFactor)
  return(object)
}

#' Enrichment level of synthetic RNA 
#'
#' @param object Enone object. 
#' @param method Which normalization methods to perform. 
#' @param log Whether to return \code{log2} values, default: TRUE. 
#'
#' @return enrichment level of synthetic RNA
#' @export
#'
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom stringr str_extract
synEnrichment <- function(object, method="TC", log=TRUE) {
  
  if (!any(SummarizedExperiment::rowData(object)$Synthetic)) {
    stop("Synthetic RNA id not provided in the object")
  }
  
  # get raw counts
  counts_df <- SummarizedExperiment::assay(object)
  
  # method.curr[1]: scaling; method.curr[2]: RUV; method.curr[3]: number of k
  method.curr <- unlist(strsplit(method, split = "_"))
  
  # check whether selected method can be provided
  if (!method.curr[1] %in% c("Raw","TC","UQ","DESeq","TMM","PossionSeq")) {
    stop("Scaling method: ", method.curr[1], " is not provided. It should be one of ", c("Raw","TC","UQ","DESeq","TMM","PossionSeq"))
  }
  if (length(method.curr) > 1 & !method.curr[2] %in% c("RUVg","RUVs","RUVse")) {
    stop("RUV method: ", method.curr[2], " is not provided. It should be one of ", c("RUVg","RUVs","RUVse"))
  }
  if (length(method.curr) > 2 & is.na(stringr::str_extract(method.curr[3],"k"))) {
    stop("Use x number of factors to estimate unwanted variation as \"kx\", but not ", method.curr[3])
  }  
  if (length(method.curr) > 2 & as.numeric(stringr::str_extract(method.curr[3],"\\d")) > ncol(object)) {
    stop("Number of required factors exceed, try least.")
  } 
  
  # scaling
  if (method.curr[1] == "Raw") {
    counts_scale <- list(dataNorm=counts_df, normFactor=rep(1,ncol(counts_df)))
  } else {
    normScaling <- get(paste0("norm", method.curr[1]))
    counts_scale <- normScaling(counts_df)
  }
  
  # RUV
  sc_idx <- CreateGroupMatrix(object$condition)
  enrich_idx <- CreateGroupMatrix(object$enrich)
  if (method.curr[2] %in% c("RUVg", "RUVs", "RUVse")) {
    
    sc.idx <- switch(method.curr[2],
                     "RUVg" = NULL,
                     "RUVs" = sc_idx,
                     "RUVse" = enrich_idx)
    counts_norm <- normRUV(counts_scale$dataNorm,
                           control.idx = getGeneSet(object, "NegControl"),
                           sc.idx = sc.idx,
                           method = method.curr[2],
                           k = as.numeric(gsub("k", "", method.curr[3])))
    
  } else {
    counts_norm <- counts_scale
  }
  
  # calculate enrichment of synthetic RNA
  counts_norm <- counts_norm$dataNorm
  syn_id <- SummarizedExperiment::rowData(object)$Synthetic
  # add 1 offset to avoid zero division (in log-scale)
  syn_en <- log2(counts_norm[syn_id, enrich_idx[1,]]+1) - log2(counts_norm[syn_id, enrich_idx[2,]]+1)
  
  if (log) {
    syn_en <- syn_en
  } else {
    syn_en <- 2^syn_en
  }
  
  return(syn_en)
}

#' Wrapper of edgeR procedure
#'
#' @param counts A un-normalized counts data matrix. 
#' @param group Vector of length p mapping the columns of \code{counts} to 
#'   corresponding samples group. 
#' @param norm.factors Vector of normalization factors with p length. 
#' @param adjust.factors Matrix with each column indicates the adjusting factors 
#'   that estimated from RUV. 
#' @param design.formula Formula
#' @param contrast.df Data frame of contrast, where extracting results as 
#'   first column vs. second column. 
#' @param coef Integer or character vector indicating which coefficients of the 
#'   linear model are to be tested equal to zero. Values must be columns or column names of design. 
#' @param logfc.cutoff Filter genes by at least X-fold difference (log2-scale) 
#'   between the two groups of samples, default: 1. 
#' @param p.cutoff Filter genes by no more than Y adjusted p-value, default: 0.05. 
#' @param only.pos Only return positive genes in filtered results \code{res.sig.ls}, default: TRUE.
#' 
#' @return List containing differential analysis object, result table and filtered result table.  
#' @export
#'
#' @import edgeR
#' @import dplyr
#' @importFrom stats model.matrix
#' @importFrom limma makeContrasts
#' @importFrom tibble rownames_to_column
#' @importFrom pbapply pblapply
edgeRDE <- function(counts, 
                    group,  
                    norm.factors = NULL, 
                    adjust.factors = NULL, 
                    design.formula = NULL, 
                    contrast.df = NULL,
                    coef = NULL,
                    logfc.cutoff = 1, p.cutoff = 0.05, only.pos = TRUE) {
  
  degs <- edgeR::DGEList(counts, group = group)
  
  if (is.null(norm.factors)) {
    degs <- edgeR::calcNormFactors(degs, method = "RLE") # Default: perform RLE normalization
  } else {
    degs$samples$norm.factors <- norm.factors
  }
  
  if (is.null(adjust.factors)) {
    design.df <- data.frame(condition = group, row.names = colnames(counts))
    design.mat <- stats::model.matrix(design.formula, data = design.df)
  } else {
    design.df <- data.frame(condition = group, adjust.factors, row.names = colnames(counts))
    design.mat <- stats::model.matrix(as.formula(paste("~0+", paste(colnames(design.df), collapse = "+"))), data = design.df)
  }
  
  # estimate dispersion and fit GLM with progress bar
  cat("- Estimate dispersion & Fit GLM... \n")
  
  de.ls1 <- pbapply::pblapply(1, function(i) {
    degs <- edgeR::estimateDisp(degs, design = design.mat)
    fit.glm <- edgeR::glmFit(degs, design = design.mat)
    return(list(degs=degs, fit.glm=fit.glm))
  })
  # unlist
  degs <- de.ls1[[1]]$degs
  fit.glm <- de.ls1[[1]]$fit.glm
  
  # Testing DEGs 
  cat("- Testing differential genes... \n")
  if (is.null(coef) & !is.null(contrast.df)) {
    # construct contrast
    contrast.vec <- apply(contrast.df, 1, function(x) { paste(paste0("condition",x), collapse="-") })
    # progress bar
    de.ls2 <- pbapply::pblapply(1:length(contrast.vec), function(i) {
      contrast.mat <- limma::makeContrasts(contrasts = contrast.vec[i], levels=design.mat)
      # LRT 
      lrt.glm <- edgeR::glmLRT(fit.glm, contrast = contrast.mat)
      
      # extract DE results
      res1 <- edgeR::topTags(lrt.glm, n = Inf, adjust.method = "BH")
      res.tab <- res1$table %>% tibble::rownames_to_column(var = "GeneID")
      
      # filter significant DEGs
      if (only.pos) {
        res.sig.tab <- res.tab[res.tab$logFC >= logfc.cutoff & res.tab$FDR < p.cutoff,]
      } else {
        res.sig.tab <- res.tab[abs(res.tab$logFC) >= logfc.cutoff & res.tab$FDR < p.cutoff,]
      }
      return(list(res.tab = res.tab, res.sig.tab = res.sig.tab))
    })
    # unlist 
    res.ls <- lapply(de.ls2, function(x) x$res.tab)
    res.sig.ls <- lapply(de.ls2, function(x) x$res.sig.tab)
    names(res.ls) <- names(res.sig.ls) <- gsub("-","_",gsub("condition","",contrast.vec))
  } 
  if (!is.null(coef) & is.null(contrast.df)) { 
    de.ls2 <- pbapply::pblapply(1, function(i) {
      lrt.glm <- edgeR::glmLRT(fit.glm, coef = coef)
      res1 <- edgeR::topTags(lrt.glm, n = Inf, adjust.method = "BH")
      res.tab <- res1$table %>% tibble::rownames_to_column(var = "GeneID")
      # filter significant DEGs
      if (only.pos) {
        res.sig.tab <- res.tab[res.tab$logFC >= logfc.cutoff & res.tab$FDR < p.cutoff,]
      } else {
        res.sig.tab <- res.tab[abs(res.tab$logFC) >= logfc.cutoff & res.tab$FDR < p.cutoff,]
      }
      return(list(res.tab = res.tab, res.sig.tab = res.sig.tab))
    })
    res.ls <- lapply(de.ls2, function(x) x$res.tab)
    res.sig.ls <- lapply(de.ls2, function(x) x$res.sig.tab)
  }
  
  return(list(de.obj = degs, res.ls = res.ls, res.sig.ls = res.sig.ls))
}

#' Create a matrix for RUVSeq
#'
#' @param group.vec A vector indicating membership in a group.
#'
#' @return A matrix. 
#' @export
#'
#' @examples 
#' CreateGroupMatrix(c('a','b','b','c','c','c','a','d','d'))
CreateGroupMatrix <- function(group.vec) {
  group.vec <- factor(group.vec)
  group.mat <- matrix(-1, nrow = length(levels(group.vec)), ncol = max(table(group.vec)))
  for (i in 1:length(levels(group.vec))) {
    idxs <- which(group.vec == levels(group.vec)[i])
    group.mat[i, 1:length(idxs)] <- idxs
  }
  group.mat
}

#' Count numbers of each members
#'
#' @param group.vec Vector of members. 
#'
#' @return vector
#' @export
#'
#' @examples 
#' countReplicate(c('a','b','b','c','c','c','a','d','d'))
countReplicate <- function(group.vec) {
  group.vec <- factor(group.vec, levels = unique(group.vec)) # keep factor levels input order
  rep.vec <- vector("double")
  for (i in levels(group.vec)) {
    curr.group <- grep(i, group.vec, value = TRUE)
    curr.idx <- grep(i, group.vec)
    curr.reps <- which(curr.group == i)
    names(curr.reps) <- curr.idx
    rep.vec <- c(rep.vec, curr.reps)
  }
  rep.vec <- rep.vec[as.character(sort(as.integer(names(rep.vec))))] # sort as input order
  names(rep.vec) <- NULL
  return(rep.vec)
}


#' Combine list of DE results
#'
#' @param res.ls Named list of differential analysis results tables. 
#' Each elements in the list correspond to a table of differential analysis 
#' results between two groups of samples. 
#' @param logfc.col Column name of the log fold-change. 
#' @param levels Factor levels of the groups, default order by the element order of \code{res.ls}. 
#'
#' @return data.frame
#' @export
#'
#' @import dplyr
reduceRes <- function(res.ls, logfc.col, levels=names(res.ls)) {
  df <- data.frame()
  for (id in names(res.ls)) {
    curr <- res.ls[[grep(id, names(res.ls), value=TRUE)]] 
    df1 <- curr %>% 
      dplyr::mutate(Group = factor(rep(id, nrow(curr)), levels = levels)) %>% 
      dplyr::select(GeneID, !!sym(logfc.col), Group)
    df <- rbind(df, df1)
  }
  return(df)
}

##--Visualization--##

#' PCA plot from counts matrix
#' 
#' @details Perform PCA based on matrix using \code{prcomp}, and visualized with scatter plot.
#' 
#' @param object A count matrix.
#' @param use.pc Which two PCs to be used, default PC1 in x-axis and PC2 in y-axis.
#' @param color Vector indicates the color mapping of samples, default NULL.
#' @param label Vector of sample names or labels, default NULL.
#' @param shape Vector indicates the shape mapping of samples, default NULL.
#' @param title Plot title, default NULL. 
#' @param vst.norm Whether to perform \code{vst} transformation, default FALSE.
#' @param palette The color palette for different groups.
#' @param repel Whether to use \code{ggrepel} to avoid overlapping text labels or not, default TRUE.
#'
#' @return ggplot2 object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom DESeq2 vst
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats prcomp
#' @importFrom paintingr paint_palette
PCAplot <- function(object, use.pc=c(1,2),
                  color=NULL, label=NULL, shape=NULL, title=NULL,
                  vst.norm=FALSE, palette=NULL, repel=TRUE) {
  if (vst.norm) {
    counts_norm <- DESeq2::vst(as.matrix(object))
  } else {
    counts_norm <- object
  }
  
  # perform PCA
  pca <- prcomp(t(counts_norm))
  pc.var <- round(summary(pca)$importance[2,], 3)
  pca_dat <- as.data.frame(pca$x)
  
  # check if use.pc exceed the range of pcs
  use.pc <- paste0("PC", use.pc)
  if (!all(use.pc %in% colnames(pca_dat))) {
    stop(use.pc, "exceed the range of PCs.")
  }
  # mapping data
  var.ls <- list(color = color,
                 shape = shape
  )
  var.length <- unlist(lapply(var.ls, length))
  var.ls <- var.ls[var.length == max(var.length)]
  map_df <- as.data.frame(Reduce(cbind, var.ls))
  colnames(map_df) <- names(var.ls)
  # combine with pca_dat if not empty
  if (!any(dim(map_df) == 0)) {
    # re-assign class
    for (i in 1:ncol(map_df)) {
      class(map_df[,i]) <- class(var.ls[[i]])
    }
    pca_dat <- cbind(pca_dat, map_df)
  }
  
  # generate color palette
  if (is.null(palette)) {
    palette <- paintingr::paint_palette("Spring", length(unique(pca_dat$color)), "continuous")
  }
  
  # create aes mapping
  map_ls <- list(x = use.pc[1],
                 y = use.pc[2],
                 color = "color",
                 shape = "shape")
  mapping <- do.call(ggplot2::aes_string, map_ls)
  
  p <- ggplot(pca_dat, mapping) +
    geom_point(size=3) +
    geom_vline(xintercept=0, color="grey80", lty=2) +
    geom_hline(yintercept=0, color="grey80", lty=2) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "top",
          axis.text = element_text(color="black")) +
    labs(x=paste0(use.pc[1], ": ", pc.var[ use.pc[1] ]*100, "%"),
         y=paste0(use.pc[2], ": ", pc.var[ use.pc[2] ]*100, "%"),
         title=title)
  
  # tune color
  if (is.character(color)) {
    p <- p + scale_color_manual(values = palette)
  } else if (is.numeric(color)) {
    p <- p + scale_color_gradientn(colors = palette)
  }

  # add text label
  if (!is.null(label)) {
    if (repel) {
      p <- p + ggrepel::geom_text_repel(label=label, max.overlaps = 20, color="black")
    } else {
      p <- p + geom_text(label=label, color="black")
    }
  }
  return(p)
}


#' Biplot of individuals and variables
#'
#' @param object PCA object returned by \code{prcomp}. 
#' @param score Vector of performance scores from \code{Enone} evaluation score. 
#' @param pt.label Whether to plot the point labels, default: TRUE. 
#' @param interactive Whether to demonstrate the plot interactively, default: FALSE. 
#'
#' @return plot
#' @export
#'
#' @import ggplot2
#' @importFrom paintingr paint_palette
#' @importFrom ggrepel geom_text_repel
#' @importFrom plotly ggplotly
PCA_Biplot <- function(object, score, pt.label=TRUE, interactive=FALSE) {
  
  # get data matrix 
  X <- object$x %*% solve(object$rotation)
  
  pc.var <- round(summary(object)$importance[2,], 3)
  pc.score <- as.data.frame(object$x)
  pc.score$method.id <- rownames(pc.score)
  pc.score$Performance <- score
  
  p <- ggplot(pc.score, aes(PC1, PC2, color=Performance, text=method.id)) +
    geom_point(size=3) +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.text = element_text(color="black")) +
    geom_hline(yintercept=0, lty="dashed") +
    geom_vline(xintercept=0, lty="dashed") +
    scale_color_gradientn(colors = paint_palette("Vesuvius", 100, "continuous")) +
    labs(x=paste0("Dim1: ", pc.var[1]*100, "%"),
         y=paste0("Dim2: ", pc.var[2]*100, "%"))
  
  # add arrow, inspired by: 
  # https://stats.stackexchange.com/questions/276645/arrows-of-underlying-variables-in-pca-biplot-in-r
  for (i in 1:ncol(X)) {
    x.cord <- cor(X[,i], object$x[,1]) * sqrt(nrow(X)-1) * 0.4
    y.cord <- cor(X[,i], object$x[,2]) * sqrt(nrow(X)-1) * 0.4
    p <- p +   
      annotate("segment", x=0, y=0,
               xend=x.cord*0.9, 
               yend=y.cord*0.9, 
               arrow=arrow(), color="#4F99B4") +
      annotate("text", x=x.cord, y=y.cord, 
               label=colnames(X)[i], color="#4F99B4") 
  }
  
  if (pt.label) {
    p <- p + ggrepel::geom_text_repel(aes(label=method.id), max.overlaps = 20, color="black")
  }
  
  if (interactive) {
    p <- plotly::ggplotly(p)
  }
  
  return(p)
}

#' Box-violin plot comparing values between groups
#'
#' @param data A data frame (or a tibble).
#' @param x The grouping variable from the \code{data}.
#' @param y The value variable from the \code{data}.
#' @param color The color variable from the \code{data}.
#' @param palette The color palette for different groups.
#' @param test Perform wilcoxon rank sum test or t-test or no test, must be one of c("wilcox.test", "t.test", "none").
#' @param add.p Label p-value or adjusted p-value, must be one of c("p", "p.adj").
#' @param step.increase Numeric vector with the increase in fraction of total height for every additional comparison to minimize overlap.
#' @param comparisons	A list of length-2 vectors specifying the groups of interest to be compared. For example to compare groups "A" vs "B" and "B" vs "C", the argument is as follow: comparisons = list(c("A", "B"), c("B", "C"))
#' @param title Plot title, default NULL. 
#' @return ggplot2 object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom rstatix wilcox_test t_test adjust_pvalue p_format add_xy_position
#' @importFrom paintingr paint_palette
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom stats as.formula
BetweenStatPlot <- function(data, x, y, color, palette = NULL,
                            test = c("wilcox.test", "t.test", "none"),
                            add.p = c("p", "p.adj"),
                            comparisons = NULL,
                            step.increase = 0.3,
                            title = NULL) {
  stat.formula <- stats::as.formula(paste(y, "~", x))

  test <- match.arg(test, choices = c("wilcox.test", "t.test", "none"))
  if (test != "none") {
    if (test == "wilcox.test") {
      stat_dat <- rstatix::wilcox_test(data, stat.formula, comparisons = comparisons)
    }
    if (test == "t.test") {
      stat_dat <- rstatix::t_test(data, stat.formula, comparisons = comparisons)
    }
    add.p <- match.arg(add.p, choices = c("p", "p.adj"))
    stat_dat <- stat_dat %>%
      rstatix::adjust_pvalue() %>%
      rstatix::p_format(!!dplyr::sym(add.p), digits = 2, leading.zero = FALSE,
               trailing.zero = TRUE, add.p = TRUE, accuracy = 2e-16) %>%
      rstatix::add_xy_position(x = x, dodge=0.8, step.increase=step.increase)
  }
  data[,x] <- as.factor(data[,x])
  x.labs <- paste0(levels(data[,x]), "\n(n=", tabulate(data[,x]),")")
  x.num <- length(unique(data[,color])) # number of x types
  if (is.null(palette)) palette <- paintingr::paint_palette("Spring", x.num, "continuous")

  p <- data %>%
    ggplot(aes_string(x, y, color = color)) +
    geom_violin(width = 0.8) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text = element_text(color="black")) +
    scale_color_manual(values = palette) +
    scale_x_discrete(labels = x.labs) +
    labs(x="", title=title)

  if (exists("stat_dat")) {
    p <- p + ggpubr::stat_pvalue_manual(data = stat_dat, label = add.p, tip.length = 0.01, size = 3)
  }

  return(p)
}

#' Dot-plot with mean_sd bar
#'
#' @param data A data.frame (or a tibble).
#' @param x The grouping variable from the \code{data}.
#' @param y The value variable from the \code{data}.
#' @param fill The fill variable from the \code{data}.  
#' @param palette The fill palette for different groups.
#'
#' @return ggplot2 object
#' @export
#' 
#' @import ggplot2
#' @importFrom paintingr paint_palette
DotPlot <- function(data, x, y, fill = NULL, palette = NULL) {
  
  if (!is.null(fill) & is.null(palette)) {
    palette <- paintingr::paint_palette("Splash",length(unique(data[[fill]])),"continuous")
  }
  
  if (!is.null(fill)) {
    dp <- ggplot(data, aes_string(x, y, fill = fill)) +
      scale_fill_manual(values = palette)
  } else {
    dp <- ggplot(data, aes_string(x, y))
  }
  
  dp <- dp +
    geom_dotplot(binaxis = "y", stackdir = "center", color = NA, 
                 dotsize = 0.6, position = "dodge") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.3, size = 0.3,
                 position = position_dodge(width=0.9), show.legend = FALSE) +
    stat_summary(fun.data = .mean_sd, width = 0.15, geom = "errorbar", size = 0.8,
                 position = position_dodge(width=0.9), show.legend = FALSE) +
    theme_classic() +
    theme(strip.background = element_blank()) +
    theme(axis.text = element_text(color = 'black')) +
    labs(x="")
  return(dp)
}

#' Statistics summary (mean and +/- sd)
#'
#' @param x Value
#'
#' @return Vector of mean and mean +/- sd 
#'
#' @importFrom stats sd
.mean_sd <- function(x) {
  m <- mean(x)
  ymin <- m - stats::sd(x)
  ymax <- m + stats::sd(x)
  return(c(y=m, ymin=ymin, ymax=ymax))
}

#' Filter Low Expressed Genes
#'
#' @param x Enone object
#' @param group Vector or factor giving group membership for a oneway layout, 
#' if appropriate, default: NULL. 
#' @param min.count Minimum count required for at least some samples, default: 10. 
#'
#' @return updated Enone object
#' @export
#'
FilterLowExprGene <- function(x, group=NULL, min.count=10) {
  # get counts
  data <- x@assays@data@listData[[1]]
  
  # size of minimum samples
  if (is.null(group)) {
    min_sample <- ncol(data)
  } else {
    group <- as.factor(group)
    n_sample <- tabulate(group)
    min_sample <- min(n_sample[n_sample > 0])
  }
  
  # filter low genes
  keep <- rowSums(data >= min.count) >= min_sample
  x <- x[keep, ]
  return(x)
  
}

#' Outlier Test
#' 
#' @description Rosner’s outlier test on principal component 1 to assess and 
#' remove potential outliers.
#'
#' @param x Enone object
#' @param return Whether to return object or simply perform test, default FALSE
#' @param remove Whether to remove outliers. If TRUE, must be paired with 
#' "return=TRUE", default: FALSE
#'
#' @return updated Enone object
#' @export
#'
#' @importFrom DESeq2 vst
#' @importFrom EnvStats rosnerTest
OutlierTest <- function(x, return=FALSE, remove=FALSE) {
  # get counts
  data <- as.matrix(x@assays@data@listData[[1]])
  
  # variance stabilizing transformation
  vt <- DESeq2::vst(data)
  
  # pricipal component analysis
  pc <- prcomp(t(vt))
  
  # Rosner’s outlier test on principal component 1
  test <- EnvStats::rosnerTest(pc$x[,1])$all.stats
  
  if (sum(test$Outlier) > 0 & return & remove) {
    outlier.idx <- test[test$Outlier == TRUE, ]$Obs.Num
    x <- x[, -outlier.idx]
  }
  
  # print test message
  cat("Rosner's outlier test\n")
  print(test)
  
  if (return) {
    return(x)
  }
}

# For adjusting no visible binding
## reduceRes
utils::globalVariables(c("GeneID", "Group"))
## ggPCA
utils::globalVariables(c("PC1", "PC2", "group"))
## ggPCA_Biplot
utils::globalVariables(c("Performance", "method.id"))
## edgeRDE
utils::globalVariables(c("rowname"))