qusageGen <- function(resids, labels, estimates, dof, std_errors, gene_sets,
                      var.equal = TRUE) {
  if (var.equal) {
    labels <- rep("Resids", ncol(resids))
  }
  if (nrow(resids) != length(estimates)) {
    return("Error: Number of rows in residual matrix do not equal length of
           estimate vectors")
  }
  if (nrow(resids) != length(dof)) {
    return("Error: Number of rows in residual matrix do not equal length of dof
           vectors")
  }
  if (nrow(resids) != length(std_errors)) {
    return("Error: Number of rows in residual matrix do not equal length of
           std_errors vectors")
  }
  names(estimates) <- rownames(resids)
  names(dof) <- rownames(resids)
  names(std_errors) <- rownames(resids)
  qlist <- list(mean = estimates, SD = std_errors, dof = dof, labels = labels)
  results <- newQSarray(qlist)
  cat("Aggregating gene data for gene sets.")
  results <- aggregateGeneSet(results, gene_sets, n.points = 2 ^ 14)
  cat("Done. \nCalculating VIF's on residual matrix.")
  results <- calcVIF(resids, results, useCAMERA = FALSE)
  cat("\nQ-Gen analysis complete.")
  results
  }

#' Run Qusage algorithm using gene level statistics
#'
#' @param object object generated from \code{genModelResults}
#' @param gene_sets list of gene_sets
#' @details This function takes the gene level comparison estimates and test
#'   statistics contained in the object produced from
#'   \code{\link{genModelResults}} and runs the Qusage algorithm across all of
#'   the comparisons. The VIFs are estimated using the raw residuals, which are
#'   also contained in the output of \code{\link{genModelResults}}.
#' @return \code{qusage_results} tall formatted matrix of results
#' @return \code{lowerCI} Matrix of gene level lower 95\% confidence intervals
#' @return \code{upperCI} Matrix gene level upper 95\% confidence intervals
#' @examples
#' # Example data
#' data(tb.expr)
#' data(tb.design)
#' 
#' # Use first 100 probes to demonstrate
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
#' model.results <- genModelResults(design_info = des.info, object = fit2, lm_Fit = fit, 
#'                                method = "limma")
#'                                
#' # Run qusage on baylor modules                             
#' data(modules)
#' qus <- qBart(model.results, modules)
#' @export
qBart <- function(object, gene_sets) {
  results <- object$results
  mod_overlap <- sapply(gene_sets, function(x) {
    sum(x %in% results$PROBE_ID)
  })
  module_list <- gene_sets[mod_overlap > 3]
  comparisons <- results[, grep("^Estimate", colnames(results))]
  tstats <- results[, grep("^Test.statistic", colnames(results))]
  stde <- abs(comparisons / tstats)
  for (i in 1:ncol(stde)) {
    stde[, i][stde[, i] < (10) ^ -6] <- min(stde[, i][stde[, i] > (10 ^ -6)])
  }
  colnames(stde) <- paste("Std_error.", colnames(stde), sep = "")
  df <- as.matrix(results[, grep("DF.", colnames(results))])
  my_df <- list()
  for (i in 1:dim(df)[2]) {
    my_df[[i]] <- df[, i]
  }
  lowerCI <- upperCI <- data.frame(matrix(nrow = nrow(comparisons),
                                          ncol = ncol(comparisons)))
  for (i in 1:ncol(comparisons)) {
    lowerCI[, i] <- comparisons[, i] - qt(0.975, my_df[[i]][1]) * stde[, i]
    upperCI[, i] <- comparisons[, i] + qt(0.975, my_df[[i]][1]) * stde[, i]
  }
  colnames(lowerCI) <- colnames(upperCI) <- sub("Estimate of ", "",
                                                colnames(comparisons))
  rownames(lowerCI) <- rownames(upperCI) <- rownames(comparisons)
  final_residual <- as.matrix(object$resids)
  comparisons <- as.matrix(comparisons)
  std_error <- as.matrix(stde)
  q_results <- list()
  for (i in 1:ncol(comparisons)) {
    q_results[[i]] <- qusageGen(resids = final_residual,
                                estimates = comparisons[, i], dof = my_df[[i]],
                                std_errors = std_error[, i],
                                gene_sets = module_list, var.equal = TRUE)
  }
  names(q_results) <- sub("Estimate of ", "", colnames(comparisons))
  for (i in 0:(length(q_results) - 1)) {
    mytable <- qsTable(q_results[[i + 1]], number = length(module_list))
    mytable <- mytable[order(as.numeric(rownames(mytable))), ]
    myCI <- calcBayesCI(q_results[[i + 1]], low = 0.025, up = 0.975,
                        addVIF = !is.null(q_results[[i + 1]]$vif))
    final_table <- cbind(mytable, t(myCI))
    names(final_table) <- paste(names(q_results[i + 1]), names(final_table),
                                sep = "_")
    if (i == 0) {
      master <- cbind(final_table[order(final_table[, 1]), ],
                      rep(names(q_results[i + 1]), dim(final_table)[1]))
      names(master) <- c(paste(names(mytable), sep = ""), "low", "up",
                         "Comparison")
    }
    if (i > 0) {
      dummy <- cbind(final_table[order(final_table[, 1]), ],
                     rep(names(q_results[i + 1]), dim(final_table)[1]))
      names(dummy) <- c(paste(names(mytable), sep = ""), "low", "up",
                        "Comparison")
      master <- rbind(master, dummy)
    }
  }
  qusage_results <- master
  z <- list(qusage_results = qusage_results, lowerCI = lowerCI,
            upperCI = upperCI)
  return(z)
}

#' Run ROAST method and format results for BART
#'
#' @param design_info design_info list generated from \code{desInfo}
#' @param gene_sets list of gene sets
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
rBart <- function(design_info, gene_sets, design, contrast, block = NULL,
                  correlation = NULL, gene.weights = NULL, var.prior = NULL,
                  df.prior = NULL, nrot = 999) {
  y <- design_info$y
  rst_mean <- rst_floormean <- rst_mean50 <- rst_msq <- list()
  contrast.t <- t(contrast)
  n <- ncol(contrast)
  for (i in 1:n) {

    rst_mean[[i]] <- mroast(y = y, index = gene_sets, design = design,
                            contrast = unname(contrast.t[i, ]), block = block,
                            correlation = correlation, var.prior = var.prior,
                            gene.weights = gene.weights, df.prior = df.prior,
                            nrot = nrot)
    rst_floormean[[i]] <- mroast(y = y, index = gene_sets, design = design,
                                 contrast = unname(contrast.t[i, ]),
                                 correlation = correlation, df.prior = df.prior,
                                 gene.weights = gene.weights, block = block,
                                 var.prior = var.prior, nrot = nrot,
                                 set.statistic = "floormean")
    rst_mean50[[i]] <- mroast(y = y, index = gene_sets, design = design,
                              contrast = unname(contrast.t[i, ]),
                              correlation = correlation, df.prior = df.prior,
                              gene.weights = gene.weights, block = block,
                              var.prior = var.prior, nrot = nrot,
                              set.statistic = "mean50")
    rst_msq[[i]] <- mroast(y = y, index = gene_sets, design = design,
                           contrast = unname(contrast.t[i, ]), block = block,
                           correlation = correlation, var.prior = var.prior,
                           gene.weights = gene.weights, df.prior = df.prior,
                           nrot = nrot, set.statistic = "msq")
  }
  names(rst_mean) <- names(rst_floormean) <- names(rst_mean50) <-
    names(rst_msq) <- colnames(contrast)
  roast_results <- list(mean = rst_mean, floormean = rst_floormean,
                        mean50 = rst_mean50, msq = rst_msq)
  return(roast_results)
}
