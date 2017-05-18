#' Generate formatted results file
#' 
#' Generate formatted results file from completed differential analysis in 
#' limma, DESeq2, or edgeR
#' @param design_info list generated from \code{\link{desInfo}}
#' @param object result object generated from \code{limma}, \code{DESeq2}, or 
#'   \code{edgeR}
#' @param lm_Fit linear model fit object generated from \code{limma}, 
#'   \code{DESeq2}, or \code{edgeR}
#' @param method Which method is being used? ("limma","deseq2","edgeR")
#' @param comp_names comp_names=NULL is the default. A vector of comparison 
#'   (contrast) names
#' @param var_names var_names=NULL is the default. Optional vector of variable 
#'   names. Otherwise, the rownames of the expression dataframe are used
#' @param var_symbols var_symbols=NULL is the default. Optional vector of 
#'   secondary annotations for the variables. Otherwise, the rownames of the 
#'   expression dataframe are used
#' @details This function formats the results obtained from running differential
#'   analysis in either one of \code{limma}, \code{DESeq2}, or \code{edgeR}. The
#'   parameter \code{object} accepts as input a results object from functions 
#'   \code{\link[limma]{eBayes}} (limma), \code{\link[DESeq2]{results}} 
#'   (DESeq2), \code{\link[edgeR]{glmLRT}} or \code{\link[edgeR]{glmQLFTest}} 
#'   (edgeR). The parameter \code{lm_Fit} accepts a fitted model object from 
#'   functions \code{\link[limma]{lmFit}} (limma), \code{\link[DESeq2]{DESeq}} 
#'   (DESeq2), \code{\link[edgeR]{glmFit}} or \code{\link[edgeR]{glmQLFit}} 
#'   (edgeR).
#' @return \code{data_type} String denoting the type of data that was analyzed
#' @return \code{results} The formatted results file
#' @return \code{resids} Matrix of residuals. Returned only if data_type="micro"
#'   or "rna". Used to estimate VIFs when running the QUSAGE algorithim in 
#'   \code{\link{qBart}}
#' @examples
#' # Example data
#' data(tb.expr)
#' data(tb.design)
#' 
#' # Only use first 100 probes to demonstrate
#' dat <- tb.expr[1:100,]
#' 
#' # Create desInfo object
#' des.info <- desInfo(y = dat, design = tb.design, data_type = "micro", 
#'                     columnname = "columnname", long = TRUE, patient_id = "monkey_id",
#'                     baseline_var = "timepoint", baseline_val = 0, time_var = "timepoint", 
#'                     responder_var = "clinical_status", sample_id = "sample_id", 
#'                     project_name = "TB")
#' 
#' # Generate lmFit and eBayes (limma) objects needed for genModelResults
#' tb.design$Group <- paste(tb.design$clinical_status,tb.design$timepoint,
#'                          sep = "")
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
#' model.results <- genModelResults(design_info = des.info, object = fit2, lm_Fit = fit, 
#'                                method = "limma")
#' @import limma
#' @import qusage
#' @importFrom SummarizedExperiment assays
#' @export
genModelResults <- function(design_info, object, lm_Fit, method = "limma",
                            comp_names = NULL, var_names = NULL,
                            var_symbols = NULL) {
  y <- design_info$y
  data_type <- design_info$data_type
  if (method != "limma") {
    if (length(object) > 1) {
      estimates <- list()
      zstats <- list()
      pvals <- list()
      if (method == "deseq2") {
        for (i in 1:length(object)) {
          object[[i]] <- as.data.frame(object[[i]])
          estimates[[i]] <- object[[i]]$log2FoldChange
          zstats[[i]] <- object[[i]]$stat
          pvals[[i]] <- object[[i]]$pvalue
        }
      }
      if (method == "edgeR") {
        for (i in 1:length(object)) {
          object[[i]] <- as.data.frame(object[[i]])
          estimates[[i]] <- object[[i]]$logFC
          zstats[[i]] <- object[[i]]$F
          pvals[[i]] <- object[[i]]$PValue
        }
      }
      estimates <- do.call("cbind", estimates)
      zstats <- do.call("cbind", zstats)
      pvals <- do.call("cbind", pvals)
      if (is.null(comp_names)) {
        colnames(estimates) <- paste("Estimate of ", colnames(estimates),
                                     sep = "")
        colnames(zstats) <- paste("Test.statistic for ", colnames(tstats),
                                  sep = "")
        colnames(pvals) <- paste("P.Value for ", colnames(pvals),
                                 sep = "")
      }
      if (!is.null(comp_names)) {
        colnames(estimates) <- paste("Estimate of ", comp_names,
                                     sep = "")
        colnames(zstats) <- paste("Test.statistic for ", comp_names,
                                  sep = "")
        colnames(pvals) <- paste("P.Value for ", comp_names, sep = "")
      }
      if (is.null(var_names)) {
        var_names <- as.character(rownames(as.data.frame(object[[1]])))
      }
      if (is.null(var_symbols)) {
        var_symbols <- as.character(var_names)
      }
      mat <- do.call("cbind", list(estimates, zstats, pvals))
      results <- cbind(var_names, var_symbols, mat)
      colnames(results)[1:2] <- c("PROBE_ID", "SYMBOL")
      rownames(results) <- NULL
      results <- as.data.frame(results)
      for (i in 3:ncol(results)) {
        results[, i] <- as.numeric(as.character(results[, i]))
      }
    }
    if (length(object) == 1) {
      if (method == "deseq2") {
        object <- as.data.frame(object)
        estimates <- object$log2FoldChange
        zstats <- object$stat
        pvals <- object$pvalue
      }
      if (method == "edgeR") {
        object <- as.data.frame(object)
        estimates <- object$logFC
        zstats <- object$F
        pvals <- object$PValue
      }
      mat <- do.call("cbind", list(estimates, zstats, pvals))
      if (is.null(comp_names)) {
        colnames(mat)[1] <- paste("Estimate of ", colnames(mat)[1],
                                  sep = "")
        colnames(mat)[2] <- paste("Test.statistic for ", colnames(mat)[2],
                                  sep = "")
        colnames(mat)[3] <- paste("P.Value for ", colnames(mat)[3],
                                  sep = "")
      }
      if (!is.null(comp_names)) {
        colnames(mat)[1] <- paste("Estimate of ", comp_names, sep = "")
        colnames(mat)[2] <- paste("Test.statistic for ", comp_names,
                                  sep = "")
        colnames(mat)[3] <- paste("P.Value for ", comp_names, sep = "")
      }
      if (is.null(var_names)) {
        var_names <- as.character(rownames(as.data.frame(object[[1]])))
      }
      if (is.null(var_symbols)) {
        var_symbols <- as.character(var_names)
      }
      results <- cbind(var_names, var_symbols, mat)
      colnames(results)[1:2] <- c("PROBE_ID", "SYMBOL")
      rownames(results) <- NULL
      results <- as.data.frame(results)
      for (i in 3:ncol(results)) {
        results[, i] <- as.numeric(as.character(results[, i]))
      }
    }
    if (method == "deseq2") {
      fitted <- t(t(SummarizedExperiment::assays(lm_Fit)[["mu"]])/
                    DESeq2::sizeFactors(lm_Fit))
      resids <- DESeq2::counts(lm_Fit, normalized = TRUE) - fitted
      resids <- as.data.frame(resids)
      results <- list(data_type = data_type, results = results, resids = resids)
    }
    if (method == "edgeR") {
      fitted <- lm_Fit$fitted.values
      resids <- lm_Fit$counts - fitted
      resids <- as.data.frame(resids)
      results <- list(data_type = data_type,results = results, resids = resids)
    }
  } else {
    comp_count <- ncol(object)
    df <- as.data.frame(matrix(rep(object$df.total, comp_count),
                               ncol = comp_count))
    estimates <- object$coefficients
    tstats <- object$t
    pvals <- object$p.value
    if (is.null(comp_names)) {
      colnames(df) <- paste("DF.", colnames(object$coefficients),
                            sep = "")
      colnames(estimates) <- paste("Estimate of ", colnames(estimates),
                                   sep = "")
      colnames(tstats) <- paste("Test.statistic for ", colnames(tstats),
                                sep = "")
      colnames(pvals) <- paste("P.Value for ", colnames(pvals), sep = "")
    }
    if (!is.null(comp_names)) {
      colnames(df) <- paste("DF.", comp_names)
      colnames(estimates) <- paste("Estimate of ", comp_names, sep = "")
      colnames(tstats) <- paste("Test.statistic for ", comp_names,
                                sep = "")
      colnames(pvals) <- paste("P.Value for ", comp_names, sep = "")
    }
    if (is.null(var_names)) {
      var_names <- as.character(rownames(object))
    }
    if (is.null(var_symbols)) {
      var_symbols <- as.character(var_names)
    }
    results <- cbind(var_names, var_symbols, df, estimates, tstats,
                     pvals)
    colnames(results)[1:2] <- c("PROBE_ID", "SYMBOL")
    rownames(results) <- NULL
    results <- as.data.frame(results)
    for (i in 3:ncol(results)) {
      results[, i] <- as.numeric(as.character(results[, i]))
    }
    resids <- residuals.MArrayLM(lm_Fit, y)
    resids <- as.data.frame(resids)
    if (data_type %in% c("flow", "metab")) {
      results <- results[, -1]
      if (data_type == "flow") {
        colnames(results)[1] <- "Flow.variable"
      } else {
        colnames(results)[1] <- "Metab.variable"
      }
      results <- list(data_type = data_type, results = results, y = y)
    } else {
      results <- list(data_type = data_type, results = results, resids = resids)
    }
  }
  return(results)
}
