qusageGen <- function(resids, estimates, dof, std.errors, gene.sets) {
  labels <- rep("Resids", ncol(resids))
  if (nrow(resids) != length(estimates)) {
    return("Error: Number of rows in residual matrix do not equal length of
           estimate vectors")
  }
  if (nrow(resids) != length(dof)) {
    return("Error: Number of rows in residual matrix do not equal length of dof
           vectors")
  }
  if (nrow(resids) != length(std.errors)) {
    return("Error: Number of rows in residual matrix do not equal length of
           std.errors vectors")
  }
  names(estimates) <- rownames(resids)
  names(dof) <- rownames(resids)
  names(std.errors) <- rownames(resids)
  qlist <- list(mean = estimates, SD = std.errors, dof = dof, labels = labels)
  results <- newQSarray(qlist)
  cat("Aggregating gene data for gene sets.")
  results <- aggregateGeneSet(results, gene.sets, n.points = 2 ^ 14)
  cat("Done. \nCalculating VIF's on residual matrix.")
  results <- calcVIF(resids, results, useCAMERA = FALSE)
  cat("\nQ-Gen analysis complete.")
  return(results)
}

#' Run Qusage algorithm using gene level statistics
#' 
#' @param model.results object returned by \code{genModelResults}.
#' @param gene.sets list of gene sets. See \code{\link{genModelResults}} for more
#'   formatting details.
#' @param annotations A data frame of additional annotations for the gene sets. 
#'   See \code{\link{genModelResults}} for more formatting details.
#' @details This function takes the gene level comparison estimates and test 
#'   statistics contained in the object returned by 
#'   \code{\link{genModelResults}} and runs the Qusage algorithm across all of 
#'   the comparisons. The VIFs are estimated using the raw residuals, which are 
#'   also contained in the output of \code{\link{genModelResults}}.
#' @return \code{qusage.results} Tall formatted matrix of results
#' @return \code{lower.ci} Matrix of gene level lower 95\% confidence intervals
#' @return \code{upper.ci} Matrix gene level upper 95\% confidence intervals
#' @return \code{gene.sets} List of gene sets provided to \code{gene.sets}
#' @return \code{annotations} data frame of gene set annotations. Default is 
#'   NULL.
#' @examples
#' # Example data
#' data(tb.expr)
#' data(tb.design)
#' 
#' # Use first 100 probes to demonstrate
#' dat <- tb.expr[1:100,]
#' 
#' # Create desInfo object
#' meta.data <- metaData(y = dat, design = tb.design, data.type = "microarray", 
#'                     columnname = "columnname", long = TRUE, subject.id = "monkey_id",
#'                     baseline.var = "timepoint", baseline.val = 0, time.var = "timepoint", 
#'                     sample.id = "sample_id")
#' 
#' # Generate lmFit and eBayes (limma) objects needed for genModelResults
#' tb.design$Group <- paste(tb.design$clinical_status,tb.design$timepoint, sep = "")
#' grp <- factor(tb.design$Group)
#' design2 <- model.matrix(~0+grp)
#' colnames(design2) <- levels(grp)
#' dupcor <- limma::duplicateCorrelation(dat, design2, block = tb.design$monkey_id)
#' fit <- limma::lmFit(dat, design2, block = tb.design$monkey_id, 
#'                     correlation = dupcor$consensus.correlation)
#' contrasts <- limma::makeContrasts(A_20vsPre = Active20-Active0, A_42vsPre = Active42-Active0, 
#'                                   levels=design2)
#' fit2 <- limma::contrasts.fit(fit, contrasts)
#' fit2 <- limma::eBayes(fit2, trend = FALSE)
#' 
#' # Create model results object for qBart
#' model.results <- genModelResults(y = dat, data.type = "microarray", object = fit2, lm.Fit = fit, 
#'                                  method = "limma")
#'                                
#' # Run qusage on baylor modules                             
#' data(modules)
#' qus.results <- qBart(model.results, modules)
#' @export
qBart <- function(model.results, gene.sets, annotations = NULL) {
  results <- model.results$results
  rownames(results) <- results$Transcript.ID
  mod.overlap <- sapply(gene.sets, function(x) {
    sum(x %in% results$Transcript.ID)
  })
  module.list <- gene.sets[mod.overlap > 3]
  comparisons <- results[, grep("^Estimate", colnames(results))]
  colnames(comparisons) <- gsub("Estimate.of.", "", colnames(comparisons))
  tstats <- results[, grep("^Test.statistic", colnames(results))]
  stde <- abs(comparisons / tstats)
  for (i in 1:ncol(stde)) {
    if (length(which(is.na(stde[, i]))) > 0) {
      stde[, i][which(is.na(stde[, i]))] <- min(stde[, i][stde[, i] > 10 ^ -6], 
                                                na.rm = TRUE)
    }
    stde[, i][stde[, i] < (10) ^ -6] <- min(stde[, i][stde[, i] > (10 ^ -6)])
  }
  colnames(stde) <- paste("Std.error.", colnames(stde), sep = "")
  df <- results[, grep("DF.", colnames(results))]
  lim95 <- apply(df, 2, function(x) qt(p = 0.975, df = x)) * stde
  lower.ci <- comparisons - lim95
  upper.ci <- comparisons + lim95
  residuals <- as.matrix(model.results$resids)
  comparisons <- as.matrix(comparisons)
  std.error <- as.matrix(stde)
  q.results <- list()
  for (i in 1:ncol(comparisons)) {
    q.results[[i]] <- qusageGen(resids = residuals,
                                estimates = comparisons[, i], dof = df[[i]],
                                std.errors = std.error[, i],
                                gene.sets = module.list)
  }
  names(q.results) <- colnames(comparisons)
  for (i in 0:(length(q.results) - 1)) {
    mytable <- qsTable(q.results[[i + 1]], number = length(module.list))
    mytable <- mytable[order(as.numeric(rownames(mytable))), ]
    myCI <- calcBayesCI(q.results[[i + 1]], low = 0.025, up = 0.975,
                        addVIF = !is.null(q.results[[i + 1]]$vif))
    final.table <- cbind(mytable, t(myCI))
    names(final.table) <- paste(names(q.results[i + 1]), names(final.table),
                                sep = "_")
    if (i == 0) {
      master <- cbind(final.table[order(final.table[, 1]), ],
                      rep(names(q.results[i + 1]), dim(final.table)[1]))
      names(master) <- c(paste(names(mytable), sep = ""), "low", "up",
                         "Comparison")
    }
    if (i > 0) {
      dummy <- cbind(final.table[order(final.table[, 1]), ],
                     rep(names(q.results[i + 1]), dim(final.table)[1]))
      names(dummy) <- c(paste(names(mytable), sep = ""), "low", "up",
                        "Comparison")
      master <- rbind(master, dummy)
    }
  }
  qusage.results <- master
  z <- list(qusage.results = qusage.results, lower.ci = lower.ci, 
            upper.ci = upper.ci, gene.sets = gene.sets, 
            annotations = annotations)
  return(z)
}

#' Run ROAST method and format results for BART
#'
#' @param meta list generated from \code{metaData}
#' @param gene.sets list of gene sets
#' @param design design matrix
#' @param contrast numeric matrix with rows corresponding to coefficients and
#'   columns containing contrasts. May be a vector if there is only one
#'   contrast.
#' @param block vector or factor specifying a blocking variable on the samples
#' @param correlation the inter-duplicate or inter-technical replicate
#'   correlation
#' @param gene.weights numeric vector of directional (positive or negative)
#'   probewise weights. Must have length equal to the number of rows in the
#'   expression dataset.
#' @param var.prior prior value for residual variances. If not provided, this is
#'   estimated from all the data using squeezeVar.
#' @param df.prior prior degrees of freedom for residual variances. If not
#'   provided, this is estimated using squeezeVar.
#' @param nrot number of rotations used to compute the p-values.
#' @details This function runs the ROAST gene set tests across all summary set
#'   statistics ("mean", "floormean", "mean50", and "msq") and returns formatted
#'   results for BART.
#'
#' @export
rBart <- function(meta, gene.sets, design, contrast, block = NULL,
                  correlation = NULL, gene.weights = NULL, var.prior = NULL,
                  df.prior = NULL, nrot = 999) {
  y <- meta$y
  rst.mean <- rst.floormean <- rst.mean50 <- rst.msq <- list()
  contrast.t <- t(contrast)
  n <- ncol(contrast)
  for (i in 1:n) {

    rst.mean[[i]] <- mroast(y = y, index = gene.sets, design = design,
                            contrast = unname(contrast.t[i, ]), block = block,
                            correlation = correlation, var.prior = var.prior,
                            gene.weights = gene.weights, df.prior = df.prior,
                            nrot = nrot)
    rst.floormean[[i]] <- mroast(y = y, index = gene.sets, design = design,
                                 contrast = unname(contrast.t[i, ]),
                                 correlation = correlation, df.prior = df.prior,
                                 gene.weights = gene.weights, block = block,
                                 var.prior = var.prior, nrot = nrot,
                                 set.statistic = "floormean")
    rst.mean50[[i]] <- mroast(y = y, index = gene.sets, design = design,
                              contrast = unname(contrast.t[i, ]),
                              correlation = correlation, df.prior = df.prior,
                              gene.weights = gene.weights, block = block,
                              var.prior = var.prior, nrot = nrot,
                              set.statistic = "mean50")
    rst.msq[[i]] <- mroast(y = y, index = gene.sets, design = design,
                           contrast = unname(contrast.t[i, ]), block = block,
                           correlation = correlation, var.prior = var.prior,
                           gene.weights = gene.weights, df.prior = df.prior,
                           nrot = nrot, set.statistic = "msq")
  }
  names(rst.mean) <- names(rst.floormean) <- names(rst.mean50) <-
    names(rst.msq) <- colnames(contrast)
  roast.results <- list(mean = rst.mean, floormean = rst.floormean,
                        mean50 = rst.mean50, msq = rst.msq)
  return(roast.results)
}
