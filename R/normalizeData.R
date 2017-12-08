manipulateData <- function(y, x, colname, norm.method = "mean", ref.var = NULL, 
                           ref.val = NULL, long = FALSE, subject.id = NULL, 
                           keep.ref = TRUE) {
  y <- y[, match(x[, colname], colnames(y), nomatch = 0)]
  x <- x[match(colnames(y), x[, colname], nomatch = 0), ]
  if (!is.null(ref.var) & !is.null(ref.val)) {
    if (!long) {
      y.norm <- y - apply(y[, x[, ref.var] == ref.val], 1, norm.method, 
                          na.rm = TRUE)
      x.norm <- x
      if (!keep.ref) {
        y.norm <- y.norm[, x[, ref.var] != ref.val]
        x.norm <- x[match(colnames(y.norm), x[, colname], nomatch = 0), ]
      }
    }
    if (long) {
      if (!is.null(subject.id)) {
        subjects.ref <- as.character(
          unique(x[, subject.id][x[, ref.var] == ref.val])
        )
        subjects.not.ref <- as.character(
          unique(x[, subject.id][x[, ref.var] != ref.val])
        )
        subjects <- intersect(subjects.ref, subjects.not.ref)
        y.norm <- list()
        for (i in 1:length(subjects)) {
          index <- x[, subject.id] %in% subjects[i]
          index.ref <- index & x[, ref.var] == ref.val
          y.norm[[i]] <- y[, index] - y[, index.ref]
        }
        y.norm <- data.frame(do.call("cbind", y.norm))
        x.norm <- x[match(colnames(y.norm), x[, colname], nomatch = 0), ]
        if (!keep.ref) {
          y.norm <- y.norm[, x.norm[, ref.var] != ref.val]
          x.norm <- x[match(colnames(y.norm), x[, colname], nomatch = 0), ]
        }
      }
      if (is.null(subject.id)) {
        return(
          warning("Must specify subject.id when long = TRUE")
        )
      }
    }
  } else {
    y.norm <- y - apply(y, 1, norm.method, na.rm = TRUE)
    x.norm <- x
  }
  normData <- list(exprs.norm = y.norm, design.norm = x.norm)
  return(normData)
}

#' Data Normalization
#' 
#' Perform various normalizations of expression data
#' @param meta list generated from \code{metaData}
#' @param norm.method string denoting whether to normalize to the mean or median 
#'   of all samples or a control group specified in \code{meta}. Default is
#'   norm.method = "mean".  
#' @details This function performs various normalizations of the expression 
#'   data, depending on the study design and the parameters defined in 
#'   \code{\link{metaData}}. For all study designs, the data is normalized to 
#'   the mean (or median) of all the samples. For cross-sectional studies with 
#'   controls, an additional normalization to the mean (or median) of the 
#'   controls is performed. For longitudinal designs, baseline normalization (
#'   subtract out each subject's baseline) and normalization to the mean (or 
#'   median) of controls (if present) is performed. In addition, separate 
#'   normalizations on baseline samples is performed.
#' @return \code{y1b} data frame of baseline samples normalized according to 
#'   \code{norm.method}. NULL if baseline samples are not specified in 
#'   \code{meta}.
#' @return \code{y2b} data frame of baseline samples normalized to controls 
#'   according to \code{norm.method}. NULL if control samples are not specified 
#'   in \code{meta}.
#' @return \code{y1} data frame of all samples normalized according to 
#'   \code{norm.method}. 
#' @return \code{y2} data frame of all samples normalized to controls according 
#'   to \code{norm.method}. NULL if control samples are not specified in 
#'   \code{meta}.
#' @return \code{y3} data frame of all samples normalized to their baseline. 
#'   NULL if study is not longitudinal or if baseline samples are not specified 
#'   in \code{meta}.
#' @return \code{norm.method} string describing normalization method used.
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
#'                     columnname = "columnname", long = TRUE, sample.id = "sample_id",
#'                     subject.id = "monkey_id", time.var = "timepoint",
#'                     baseline.var = "timepoint", baseline.val = 0)
#' 
#' # Normalize and cluster data
#' data.norm <- normalizeData(meta = meta.data)
#' @importFrom stats as.dendrogram dist qt sd
#' @export
normalizeData <- function(meta, norm.method = "mean") {
  long <- meta$long
  subject.id <- meta$subject.id
  baseline.var <- meta$baseline.var
  baseline.val <- meta$baseline.val
  control.var <- meta$control.var
  control.val <- meta$control.val
  design <- meta$design
  exprs <- meta$y
  y1 <- y2 <- y1b <- y2b <- y3 <- NULL
  y1 <- manipulateData(y = exprs, x = design, colname = "columnname", 
                       norm.method = norm.method)$exprs.norm
  if (!is.null(control.var)) {
    y2 <- manipulateData(y = exprs, x = design, colname = "columnname", 
                         ref.var = control.var, ref.val = control.val, 
                         long = FALSE, keep.ref = TRUE)$exprs.norm
    desCase <- design[design[, control.var] != control.val, ]
    expCase <- exprs[, match(desCase$columnname, colnames(exprs), nomatch = 0)]
    if (long) {
      if (!is.null(baseline.var) & !is.null(baseline.val)) {
        desBase <- desCase[desCase[, baseline.var] == baseline.val, ]
        expBase <- exprs[, match(desBase$columnname, colnames(exprs), 
                                 nomatch = 0)]
        y1b <- manipulateData(y = expBase, x = desBase, 
                              colname = "columnname")$exprs.norm
        desBaseCtrl <- design[design[, baseline.var] == baseline.val | 
                                design[, control.var] == control.val, ]
        expBaseCtrl <- exprs[, match(desBaseCtrl$columnname, colnames(exprs), 
                                     nomatch = 0)]
        y2b <- manipulateData(y = expBaseCtrl, x = desBaseCtrl, 
                              colname = "columnname", ref.var = control.var, 
                              ref.val = control.val, long = FALSE, 
                              keep.ref = TRUE)$exprs.norm
        if (is.null(subject.id)) {
          message("subject.id is not defined. Cannot produce baseline normalized 
data.")
        } else {
          y3 <- manipulateData(y = expCase, x = desCase, colname = "columnname", 
                               ref.var = baseline.var, ref.val = baseline.val, 
                               long = TRUE, subject.id = subject.id, 
                               keep.ref = FALSE)$exprs.norm
        }
      } else {
        message("baseline.var and/or baseline.val is unspecified. Cannot produce 
baseline ", norm.method, " normalized, baseline healthy normalized, or all 
samples baseline normalized data.")
      }
    }
  }
  if (is.null(control.var)) {
    if (long) {
      if (!is.null(baseline.var) & !is.null(baseline.val)) {
        desBase <- design[design[, baseline.var] == baseline.val, ]
        expBase <- exprs[, match(desBase$columnname, colnames(exprs), 
                                 nomatch = 0)]
        y1b <- manipulateData(y = expBase, x = desBase, 
                              colname = "columnname")$exprs.norm
        if (is.null(subject.id)) {
          message("subject.id is not defined. Cannot produce baseline normalized 
data.")
        } else {
          y3 <- manipulateData(y = exprs, x = design, colname = "columnname", 
                               ref.var = baseline.var, ref.val = baseline.val, 
                               long = TRUE, subject.id = subject.id, 
                               keep.ref = FALSE)$exprs.norm
        }
      } else {
        message("baseline.var and/or baseline.val is unspecified. Cannot produce 
baseline ",norm.method, " normalized or all samples baseline normalized data.")
      }
    } 
  }
  y <- list(y1b = y1b, y2b = y2b, y1 = y1, y2 = y2, y3 = y3, norm.method = 
              norm.method)
  return(y)
}
