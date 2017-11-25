corr.test <- function (x, y = NULL, use = "pairwise", method = "pearson", 
                        adjust = "holm", alpha = 0.05, ci = TRUE) 
{
  cl <- match.call()
  if (is.null(y)) {
    r <- cor(x, use = use, method = method)
    sym <- TRUE
    n <- t(!is.na(x)) %*% (!is.na(x))
  }
  else {
    r <- cor(x, y, use = use, method = method)
    sym = FALSE
    n <- t(!is.na(x)) %*% (!is.na(y))
  }
  if ((use == "complete") | (min(n) == max(n))) 
    n <- min(n)
  t <- (r * sqrt(n - 2))/sqrt(1 - r^2)
  p <- 2 * (1 - pt(abs(t), (n - 2)))
  se <- sqrt((1 - r * r)/(n - 2))
  nvar <- ncol(r)
  p[p > 1] <- 1
  if (adjust != "none") {
    if (is.null(y)) {
      lp <- upper.tri(p)
      pa <- p[lp]
      pa <- p.adjust(pa, adjust)
      p[upper.tri(p, diag = FALSE)] <- pa
    }
    else {
      p[] <- p.adjust(p, adjust)
    }
  }
  z <- psych::fisherz(r[lower.tri(r)])
  if (ci) {
    if (min(n) < 4) {
      warning("Number of subjects must be greater than 3 to find confidence intervals.")
    }
    alpha <- 1 - alpha/2
    dif <- qnorm(alpha)
    if (sym) {
      if (is.matrix(n)) {
        se <- 1/sqrt(n[lower.tri(n)] - 3)
      }
      else {
        se <- 1/sqrt(n - 3)
      }
      lower <- psych::fisherz2r(z - dif * se)
      upper <- psych::fisherz2r(z + dif * se)
      ci <- data.frame(lower = lower, r = r[lower.tri(r)], 
                       upper = upper, p = p[lower.tri(p)])
      cnR <- abbreviate(colnames(r), minlength = 5)
      k <- 1
      for (i in 1:(nvar - 1)) {
        for (j in (i + 1):nvar) {
          rownames(ci)[k] <- paste(cnR[i], cnR[j], sep = "-")
          k <- k + 1
        }
      }
    }
    else {
      z <- psych::fisherz(r)
      se <- 1/sqrt(n - 3)
      lower <- as.vector(psych::fisherz2r(z - dif * se))
      upper <- as.vector(psych::fisherz2r(z + dif * se))
      ci <- data.frame(lower = lower, r = as.vector(r), 
                       upper = upper, p = as.vector(p))
      cnR <- abbreviate(rownames(r), minlength = 5)
      cnC <- abbreviate(colnames(r), minlength = 5)
      k <- 1
      for (i in 1:ncol(y)) {
        for (j in 1:ncol(x)) {
          rownames(ci)[k] <- paste(cnR[j], cnC[i], sep = "-")
          k <- k + 1
        }
      }
    }
  }
  else {
    ci <- NULL
  }
  result <- list(r = r, n = n, t = t, p = p, se = se, adjust = adjust, 
                 sym = sym, ci = ci, Call = cl)
  return(result)
}

#' Cross Correlations
#' 
#' Perform pairwise correlations formatted for BART
#' @param x A numeric matrix or data frame.
#' @param y A second numeric matrix or data frame with the same number of rows
#'   as \code{x}.
#' @param by Default is NULL. A grouping vector with length equal to
#'   \code{nrow(x)}. Allows correlation by group (e.g. time).
#' @param by.name Default is NULL. String denoting the name of the grouping
#'   factor.
#' @param description String giving description of what's being correlated.
#'   Default is "X vs Y".
#' @param x.var String denoting type of variables in \code{x}. Default is "X".
#' @param y.var String denoting type of variables in \code{y}. Default is "Y".
#' @param method Default is set to "pearson". The alternatives are "spearman"
#'   and "kendall".
#' @param order.by Order by p-value ("p") or correlation value ("r"). Default is
#'   "p".
#' @param decreasing logical. Should the sort be increasing or decreasing?
#'   Default is TRUE.
#' @details This function uses the \code{\link[psych]{corr.test}} function in
#'   the \code{pysch} package to find the correlations and p-values. It then
#'   formats the results in a tall format for BART.
#' @return \code{corrs} Tall dataframe of correlations and p-values.
#' @return \code{corr.files} A dataframe obtained by
#'   \code{cbind(by,x,y)}.
#' @return \code{corr.names} A string denoting the description of the
#'   variables being correlated.
#' @return \code{x.var} String denoting type of variables in \code{x}.
#' @return \code{y.var} String denoting type of variables in \code{y}.
#' @return \code{corr.method} One of "pearson", "spearman", or "kendall".
#' @importFrom stats cor p.adjust pt qnorm
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
#' corrs <- crossCorr(x = module.as, y = flow, by = time, by.name = "days", 
#'                    description = "Mod.Act.Score vs Flow", x.var = "Mod.Act.Score", 
#'                    y.var = "Flow", method = "spearman")
#' @export
crossCorr <- function(x, y, by = NULL, by.name = NULL, description = "X vs Y",
                      x.var = "X", y.var = "Y", method = "pearson",
                      order.by = "p", decreasing = TRUE) {
  if (nrow(x) != nrow(y)) {
    stop("x and y do not have the same number of rows")
  }
  x_sub <- y_sub <- corrs <- r <- p <- num <- list()
  if (is.null(by)) {
    n <- 1
    time <- 0
    corr_design <- cbind(time = time, x, y)
  } else {
    n <- length(unique(by))
    corr_design <- cbind(by, x, y)
    colnames(corr_design)[1] <- by.name
  }
  for (i in 1:n) {
    if (is.null(by)) {
      rindex <- nrow(x)
    } else {
      rindex <- which(by %in% unique(by)[i])
    }
    x_sub[[i]] <- x[rindex, ]
    y_sub[[i]] <- y[rindex, ]
    corrs[[i]] <- corr.test(x_sub[[i]], y_sub[[i]], method = method, 
                            adjust = "none", ci = FALSE)
    r[[i]] <- reshape2::melt(corrs[[i]]$r)
    if (n == 1) {
      r[[i]] <- cbind(r[[i]][, c(1, 2)], time = 0, r[[i]][, 3])
    } else {
      r[[i]] <- cbind(r[[i]][, c(1, 2)], unique(by)[i], r[[i]][, 3])
    }
    p[[i]] <- reshape2::melt(corrs[[i]]$p)
    p[[i]][, 3] <- -log10(p[[i]][, 3])
    num[[i]] <- reshape2::melt(corrs[[i]]$n)
    if (ncol(num[[i]]) == 1) {
      corrs[[i]] <- cbind(r[[i]], num[[i]][, 1], p[[i]][, 3])
    } else {
      corrs[[i]] <- cbind(r[[i]], num[[i]][, 3], p[[i]][, 3])
    }
    if (is.null(by.name)) {
      colnames(corrs[[i]]) <- c("Variable", "With", "time",
                                       paste(tools::toTitleCase(method),
                                             "_Correlation", sep = ""), "NObs",
                                       "NegLog10p")
    } else {
      colnames(corrs[[i]]) <- c("Variable", "With", by.name,
                                       paste(tools::toTitleCase(method),
                                             "_Correlation", sep = ""), "NObs",
                                       "NegLog10p")
    }
  }
  corrs <- do.call("rbind", corrs)
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
  if (order.by == "p") {
    corrs <- corrs[order(corrs[, 6], decreasing = decreasing), ]
    rownames(corrs) <- NULL
  }
  if (order.by == "r") {
    corrs <- corrs[order(corrs[, 4], decreasing = decreasing), ]
    rownames(corrs) <- NULL
  }
  corrs$Base_subtracted <- FALSE
  corrs$NegLog10p <- 1 / (10 ^ corrs$NegLog10p)
  colnames(corrs)[colnames(corrs) == "NegLog10p"] <- "Raw.P.Value"
  z <- list(corrs = corrs, corr.names = description,
            x.var = x.var, y.var = y.var, corr.method = method,
            corr.files = corr_design)
  return(z)
}
