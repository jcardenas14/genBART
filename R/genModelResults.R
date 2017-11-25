#' Generate formatted results file
#' 
#' Generate formatted results file from completed differential analysis in 
#' limma, DESeq2, or edgeR
#' @param y data frame the model was run on.
#' @param data.type type of data being analyzed ("microarray", "rnaseq", 
#'   "flow", "metab"). Default is data.type="microarray".
#' @param object result object generated from \code{limma}, \code{DESeq2}, or 
#'   \code{edgeR}.
#' @param lm.Fit linear model fit object generated from \code{limma}, 
#'   \code{DESeq2}, or \code{edgeR}
#' @param method Which method is being used? ("limma","deseq2","edgeR")
#' @param comp.names comp.names=NULL is the default. A vector of comparison 
#'   (contrast) names
#' @param var.names var.names=NULL is the default. Optional vector of variable 
#'   names. Otherwise, the rownames of the expression data frame or matrix are 
#'   used.
#' @param var.symbols var.symbols=NULL is the default. Optional vector of 
#'   annotations for the variables. Otherwise, the rownames of the 
#'   expression dataframe are used.
#' @param gene.sets A list of gene sets.
#' @param annotations A data frame of additional annotations for the gene sets. 
#' @details This function formats the results obtained from running differential
#'   analysis in either one of \code{limma}, \code{DESeq2}, or \code{edgeR}. The
#'   parameter \code{object} accepts as input a results object from functions 
#'   \code{\link[limma]{eBayes}} (limma), \code{\link[DESeq2]{results}} 
#'   (DESeq2), \code{\link[edgeR]{glmLRT}} or \code{\link[edgeR]{glmQLFTest}} 
#'   (edgeR). The parameter \code{lm.Fit} accepts a fitted model object from 
#'   functions \code{\link[limma]{lmFit}} (limma), \code{\link[DESeq2]{DESeq}} 
#'   (DESeq2), \code{\link[edgeR]{glmFit}} or \code{\link[edgeR]{glmQLFit}} 
#'   (edgeR).
#' @return \code{data.type} string denoting the type of data that was analyzed
#' @return \code{results} the formatted results returned as a data frame
#' @return \code{resids} data frame of residuals. Returned only if 
#'   data.type="microarray" or "rnaseq" and method="limma". Used to estimate VIFs 
#'   when running the QUSAGE algorithim in \code{\link{qBart}}.
#' @return \code{gene.sets} list of gene sets provided by the user. NULL if no 
#'   list provided.
#' @return \code{annotations} data frame of gene set annotations provided by the
#'   user. Null if no annotations are provided.
#' @examples
#' # Example data
#' data(tb.expr)
#' data(tb.design)
#' 
#' # Only use first 100 genes to demonstrate
#' dat <- tb.expr[1:100,]
#' 
#' # Generate lmFit and eBayes (limma) objects needed for genModelResults
#' tb.design$Group <- paste(tb.design$clinical_status,tb.design$timepoint,sep = "")
#' grp <- factor(tb.design$Group)
#' design2 <- model.matrix(~0+grp)
#' colnames(design2) <- levels(grp)
#' 
#' dupcor <- limma::duplicateCorrelation(dat, design2, block = tb.design$monkey_id)
#' fit <- limma::lmFit(dat, design2, block = tb.design$monkey_id, 
#'              correlation = dupcor$consensus.correlation)
#' contrasts <- limma::makeContrasts(A_20vsPre = Active20-Active0, A_42vsPre = Active42-Active0, 
#'                                   levels=design2)
#' fit2 <- limma::contrasts.fit(fit, contrasts)
#' fit2 <- limma::eBayes(fit2, trend = FALSE)
#' 
#' # Format results
#' model.results <- genModelResults(y = dat, data.type = "microarray", object = fit2,
#'                                  lm.Fit = fit, method = "limma")
#' @import limma
#' @import qusage
#' @export
genModelResults <- function (y, data.type, object, lm.Fit, method = "limma", 
                             comp.names = NULL, var.names = NULL, 
                             var.symbols = NULL, gene.sets = NULL, 
                             annotations = NULL) 
{
  if (method != "limma") {
    estimates <- zstats <- pvals <- fdr <- list()
    if (method == "deseq2") {
      for (i in 1:length(object)) {
        object[[i]] <- data.frame(object[[i]])
        estimates[[i]] <- object[[i]]$log2FoldChange
        zstats[[i]] <- object[[i]]$stat
        pvals[[i]] <- object[[i]]$pvalue
        fdr[[i]] <- object[[i]]$padj
      }
    }
    if (method == "edgeR") {
      for (i in 1:length(object)) {
        object[[i]] <- data.frame(object[[i]]$table)
        estimates[[i]] <- object[[i]]$logFC
        zstats[[i]] <- object[[i]][,3]
        pvals[[i]] <- object[[i]]$PValue
        fdr[[i]] <- p.adjust(object[[i]]$PValue, method = "fdr")
      }
    }
    if(length(object) == 1){
      estimates <- data.frame(estimates)
      zstats <- data.frame(zstats)
      pvals <- data.frame(pvals)
      fdr <- data.frame(fdr)
    } else {
      estimates <- data.frame(do.call("cbind", estimates))
      zstats <- data.frame(do.call("cbind", zstats))
      pvals <- data.frame(do.call("cbind", pvals))
      fdr <- data.frame(do.call("cbind", fdr))
    }
    if (is.null(comp.names)) {
      colnames(estimates) <- paste0("Estimate.", colnames(estimates))
      colnames(zstats) <- paste0("Test.statistic.", colnames(tstats))
      colnames(pvals) <- paste0("P.Value.", colnames(pvals))
      colnames(fdr) <- paste0("FDR.P.Value.", colnames(fdr))
    }
    if (!is.null(comp.names)) {
      colnames(estimates) <- paste0("Estimate.", comp.names)
      colnames(zstats) <- paste0("Test.statistic.", comp.names)
      colnames(pvals) <- paste0("P.Value.", comp.names)
      colnames(fdr) <- paste0("FDR.P.Value.", comp.names)
    }
    if (is.null(var.names)) {
      var.names <- as.character(rownames(as.data.frame(object[[1]])))
    }
    if (is.null(var.symbols)) {
      var.symbols <- as.character(var.names)
    }
    results <- cbind(var.names, var.symbols, estimates, zstats, pvals, fdr)
    colnames(results)[1:2] <- c("Transcript.ID", "Gene.Symbol")
    rownames(results) <- NULL
    results <- data.frame(results)
    for (i in 3:ncol(results)) {
      results[, i] <- as.numeric(as.character(results[, i]))
    }
    results <- list(data.type = data.type, results = results, gene.sets = 
                      gene.sets, annotations = annotations)
  } else {
    comp.count <- ncol(object)
    df <- data.frame(matrix(rep(object$df.total, comp.count), 
                            ncol = comp.count))
    estimates <- object$coefficients
    tstats <- object$t
    pvals <- object$p.value
    fdr <- apply(pvals, 2, p.adjust, method = "fdr")
    if (is.null(comp.names)) {
      colnames(df) <- paste0("DF.", colnames(object$coefficients))
      colnames(estimates) <- paste0("Estimate.", colnames(estimates))
      colnames(tstats) <- paste0("Test.statistic.", colnames(tstats))
      colnames(pvals) <- paste0("P.Value.", colnames(pvals))
      colnames(fdr) <- paste0("FDR.P.Value", colnames(fdr))
    }
    if (!is.null(comp.names)) {
      colnames(df) <- paste0("DF.", comp.names)
      colnames(estimates) <- paste0("Estimate.", comp.names)
      colnames(tstats) <- paste0("Test.statistic.", comp.names)
      colnames(pvals) <- paste0("P.Value.", comp.names)
      colnames(fdr) <- paste0("FDR.P.Value.", comp.names)
    }
    if (is.null(var.names)) {
      var.names <- as.character(rownames(object))
    }
    if (is.null(var.symbols)) {
      var.symbols <- as.character(var.names)
    }
    results <- cbind(var.names, var.symbols, df, estimates, tstats, pvals, fdr)
    colnames(results)[1:2] <- c("Transcript.ID", "Gene.Symbol")
    rownames(results) <- NULL
    results <- data.frame(results)
    for (i in 3:ncol(results)) {
      results[, i] <- as.numeric(as.character(results[, i]))
    }
    resids <- residuals.MArrayLM(lm.Fit, y)
    resids <- data.frame(resids)
    if (data.type %in% c("flow", "metab")) {
      results <- results[, -1]
      if (data.type == "flow") {
        colnames(results)[1] <- "Flow.variable"
      }
      else {
        colnames(results)[1] <- "Metab.variable"
      }
      results <- list(data.type = data.type, results = results, y = y)
    } else {
      results <- list(data.type = data.type, results = results, resids = resids, 
                      gene.sets = gene.sets, annotations = annotations)
    }
  }
  return(results)
}
