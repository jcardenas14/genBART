#' Cross Correlations
#' 
#' Perform pairwise correlations formatted for BART
#' @param x A numeric matrix or data frame.
#' @param y A second numeric matrix or data frame with the same number of rows as \code{x}.
#' @param by Default is NULL. A grouping vector with length equal to \code{nrow(x)}. Allows
#'   correlation by group (e.g. time).
#' @param by_name Default is NULL. String denoting the name of the grouping factor.
#' @param description String giving description of what's being correlated. Default is "X vs Y".
#' @param x_var String denoting type of variables in \code{x}. Default is "X".
#' @param y_var String denoting type of variables in \code{y}. Default is "Y".
#' @param method Default is set to "pearson". The alternatives are "spearman" and "kendall".
#' @param order_by Order by p-value ("p") or correlation value ("r"). Default is "p".
#' @param decreasing logical. Should the sort be increasing or decreasing? Default is TRUE.
#' @details This function uses the \code{\link[psych]{corr.test}} function in the \code{pysch}
#'   package to find the correlations and p-values. It then formats the results in a tall format for
#'   BART.
#' @return \code{correlations} Tall dataframe of correlations and p-values.
#' @return \code{correlation_files} A dataframe obtained by \code{cbind(by,x,y)}.
#' @return \code{correlation_names} A string denoting the description of the variables being
#'   correlated.
#' @return \code{x_var} String denoting type of variables in \code{x}.
#' @return \code{y_var} String denoting type of variables in \code{y}.
#' @return \code{correlation_method} One of "pearson", "spearman", or "kendall".
#' @examples
#' # Example data
#' data(tb.flow)
#' data(module.as)
#' 
#' # Define time variable
#' time <- module.as$time
#' module.as$time <- NULL
#' 
#' # Format flow data to run correlations and match flow samples with module.as 
#' flow <- data.frame(t(tb.flow))
#' flow <- flow[match(rownames(module.as), rownames(flow), nomatch = 0), ]
#' 
#' # Run correlations and format for BART
#' corrs <- crossCorr(x = module.as, y = flow, by = time, by_name = "days", 
#'                    description = "Mod.Act.Score vs Flow", x_var = "Mod.Act.Score", 
#'                    y_var = "Flow", method = "spearman")
#' @export
crossCorr <- function(x, y, by = NULL, by_name = NULL, description = "X vs Y",
                      x_var = "X", y_var = "Y", method = "pearson",
                      order_by = "p", decreasing = TRUE) {
  if (nrow(x) != nrow(y)) {
    stop("x and y do not have the same number of rows")
  }
  x_sub <- y_sub <- correlations <- r <- p <- num <- list()
  if (is.null(by)) {
    n <- 1
    time <- 0
    corr_design <- cbind(time = time, x, y)
  } else {
    n <- length(unique(by))
    corr_design <- cbind(by, x, y)
    colnames(corr_design)[1] <- by_name
  }
  for (i in 1:n) {
    if (is.null(by)) {
      rindex <- nrow(x)
    } else {
      rindex <- which(by %in% unique(by)[i])
    }
    x_sub[[i]] <- x[rindex, ]
    y_sub[[i]] <- y[rindex, ]
    correlations[[i]] <- psych::corr.test(x_sub[[i]], y_sub[[i]],
                                          method = method, adjust = "none",
                                          ci = FALSE)
    r[[i]] <- reshape2::melt(correlations[[i]]$r)
    if (n == 1) {
      r[[i]] <- cbind(r[[i]][, c(1, 2)], time = 0, r[[i]][, 3])
    } else {
      r[[i]] <- cbind(r[[i]][, c(1, 2)], unique(by)[i], r[[i]][, 3])
    }
    p[[i]] <- reshape2::melt(correlations[[i]]$p)
    p[[i]][, 3] <- -log10(p[[i]][, 3])
    num[[i]] <- reshape2::melt(correlations[[i]]$n)
    if (ncol(num[[i]]) == 1) {
      correlations[[i]] <- cbind(r[[i]], num[[i]][, 1], p[[i]][, 3])
    } else {
      correlations[[i]] <- cbind(r[[i]], num[[i]][, 3], p[[i]][, 3])
    }
    if (is.null(by_name)) {
      colnames(correlations[[i]]) <- c("Variable", "With", "time",
                                       paste(tools::toTitleCase(method),
                                             "_Correlation", sep = ""), "NObs",
                                       "NegLog10p")
    } else {
      colnames(correlations[[i]]) <- c("Variable", "With", by_name,
                                       paste(tools::toTitleCase(method),
                                             "_Correlation", sep = ""), "NObs",
                                       "NegLog10p")
    }
  }
  corrs <- do.call("rbind", correlations)
  corrs$Variable <- as.character(corrs$Variable)
  corrs$With <- as.character(corrs$With)
  if (length(which(corrs$Variable == corrs$With)) > 0) {
    corrs <- corrs[-which(corrs$Variable == corrs$With), ]
  }
  if (length(which(is.na(corrs$NegLog10p))) > 0) {
    corrs <- corrs[-which(is.na(corrs$NegLog10p)), ]
  }
  if (length(which(corrs$NegLog10p == "Inf")) > 0) {
    corrs$NegLog10p[which(corrs$NegLog10p == "Inf")] <- 312
  }
  if (order_by == "p") {
    corrs <- corrs[order(corrs[, 6], decreasing = decreasing), ]
    rownames(corrs) <- NULL
  }
  if (order_by == "r") {
    corrs <- corrs[order(corrs[, 4], decreasing = decreasing), ]
    rownames(corrs) <- NULL
  }
  corrs$Base_subtracted <- FALSE
  corrs$NegLog10p <- 1 / (10 ^ corrs$NegLog10p)
  colnames(corrs)[colnames(corrs) == "NegLog10p"] <- "Raw.P.Value"
  z <- list(correlations = corrs, correlation_names = description,
            x_var = x_var, y_var = y_var, correlation_method = method,
            correlation_files = corr_design)
  return(z)
}
