#' Declare meta data information for downstream analysis
#' 
#' Match design and expression data frames. Declare and store design parameters.
#' @param y An expression data frame.
#' @param design A sample annotation data frame containing sample information
#'   (e.g. age, condition, timepoint, etc.).
#' @param data.type Type of data that's going to be analyzed ("microarray", 
#'   "rnaseq","flow", "metab"). Default is data.type="microarray".
#' @param columnname Name of column in design that contains the column names of 
#'   y.
#' @param long logical; Is the study longitudinal?.
#' @param time.var For longitudinal studies. Name of column in design that 
#'   contains the study time points.
#' @param subject.id Name of column in design that contains ids for individual
#'   subjects.
#' @param baseline.var Name of column in design that contains values referring 
#'   to baseline observations.
#' @param baseline.val String or numeric value denoting baseline observations.
#' @param control.var Name of column in design that contains values referring to
#'   controls.
#' @param control.val String or numeric value denoting controls.
#' @param sample.id Name of column in design that contains unique sample 
#'   identification.
#' @return A list containing the matched design and expression data and all of 
#'   the design parameters specified.
#' @examples
#' # Using example data
#' data(tb.expr)
#' data(tb.design)
#' meta.data <- metaData(y = tb.expr, design = tb.design, data.type = "microarray", 
#'                     columnname = "columnname", long = TRUE, sample.id = "sample_id",
#'                     subject.id = "monkey_id", time.var = "timepoint",
#'                     baseline.var = "timepoint", baseline.val = 0)
#' @export
metaData <- function(y, design, data.type = "microarray", columnname = NULL, 
                    long = FALSE, time.var = NULL, sample.id = NULL, 
                    subject.id = NULL, baseline.var = NULL, baseline.val = NULL, 
                    control.var = NULL, control.val = NULL) 
{
  if (!data.type %in% c("microarray", "rnaseq", "flow", "metab")) {
    return(
      warning("data.type must be one of 'microarray', 'rnaseq', 'metab', or 
              'flow'.")
    )
  }
  if (is.null(columnname)) {
    return(warning("Must enter columnname to match design and expression"))
  } else {
    if (is.null(design[[columnname]])) {
      return(
        warning("columnname parameter given is not a column name in design. 
Please check spelling.")
      )
    }
    if (!any(design[[columnname]] %in% colnames(y))) {
      return(
        warning("columnname values in design do not match any columnnames in y")
      )
    }
    if (!is.null(sample.id)) {
      if (is.null(design[[sample.id]])) {
        return(
          warning("sample.id parameter given is not a column name in design. 
Please check spelling.")
        )
      }
      if (!is.null(subject.id)) {
        if (is.null(design[[subject.id]])) {
          return(
            warning("subject.id parameter given is not a column name in design. 
Please check spelling.")
          )
        }
        if (columnname == sample.id) {
          design$columnname <- as.character(design[[sample.id]])
        }
        if (columnname == subject.id) {
          design$columnname <- as.character(design[[subject.id]])
        }
        if (columnname != sample.id & columnname != subject.id) {
          colnames(design)[colnames(design) %in% columnname] <- "columnname"
        }
      } else {
        if (columnname == sample.id) {
          design$columnname <- as.character(design[[sample.id]])
        } else {
          colnames(design)[colnames(design) %in% columnname] <- "columnname"
        }
      }
    } else {
      if (!is.null(subject.id)) {
        if (is.null(design[[subject.id]])) {
          return(
            warning("subject.id parameter given is not a column name in design. 
Please check spelling.")
          )
        }
        if (columnname == subject.id) {
          design$columnname <- as.character(design[[subject.id]])
        } else {
          colnames(design)[colnames(design) %in% columnname] <- "columnname"
        }
      } else {
        colnames(design)[colnames(design) %in% columnname] <- "columnname"
      }
    }
  }
  if (nrow(design) > ncol(y)) {
    if (all(colnames(y) %in% design$columnname)) {
      message("More samples in design than y. Throwing out excess samples.")
      design <- design[match(colnames(y), design$columnname, nomatch = 0), ]
    } else {
      message("More samples in design than y, but some samples in y are not in 
design. Throwing out excess and unmatched samples from design and y.")
      design <- design[match(colnames(y), design$columnname, nomatch = 0), ]
      y <- y[, match(design$columnname, colnames(y), nomatch = 0)]
    }
  } else if (nrow(design) < ncol(y)) {
    if (all(design$columnname %in% colnames(y))) {
      message("More samples in y than design. Throwing out excess samples.")
      y <- y[, match(design$columnname, colnames(y), nomatch = 0)]
    } else {
      message("More samples in y than design, but some samples in design are not 
in y. Throwing out excess and unmatched samples from design and y.")
      y <- y[, match(design$columnname, colnames(y), nomatch = 0)]
      design <- design[match(colnames(y), design$columnname, nomatch = 0), ]
    }
  } else {
    if (all(design$columnname %in% colnames(y))) {
      y <- y[, match(design$columnname, colnames(y), nomatch = 0)]
    } else {
      message("There are an equal number of samples in design and y, but some 
sample names do not match. Throwing out unmatched samples.")
      y <- y[, match(design$columnname, colnames(y), nomatch = 0)]
      design <- design[match(colnames(y), design$columnname, nomatch = 0), ]
    }
  }
  if (is.null(control.var) || is.null(control.val)) {
    if (!is.null(control.var) & is.null(control.val)) {
      if (is.null(design[[control.var]])) {
        return(
          warning("control.var parameter given is not a column name in design. 
Please check spelling.")
        )
      }
      return(
        warning("control.var specified but not control.val. Please specify 
control.val")
      )
    }
    if (!is.null(control.val) & is.null(control.var)) {
      return(
        warning("control.val specified but not control.var. Please specify 
control.var")
      )
    }
  } 
  if (long) {
    if (is.null(time.var)) {
      warning("long is TRUE but time.var not specified")
    } else {
      if (is.null(design[[time.var]])) {
        return(
          warning("time.var parameter given is not a column name in design. 
Please check spelling.")
        )
      }
    }
    if (is.null(baseline.var)) {
      warning("long is TRUE but baseline.var not specified")
    } else {
      if (is.null(design[[baseline.var]])) {
        return(
          warning("baseline.var parameter given is not a column name in design. 
Please check spelling.")
        )
      }
    }
    if (is.null(baseline.val)) {
      warning("long is TRUE but baseline.val not specified")
    }
  } else {
    if (!is.null(time.var)) {
      if (is.null(design[[time.var]])) {
        return(
          warning("time.var parameter given is not a column name in design. 
Please check spelling.")
        )
      }
      message("time.var specified. Setting long to TRUE.")
      if (is.null(baseline.var) || is.null(baseline.val)) {
        if (!is.null(baseline.var) & is.null(baseline.val)) {
          if (is.null(design[[baseline.var]])) {
            return(
              warning("baselie.var parameter given is not a column name in 
design. Please check spelling.")
            )
          }
          return(
            warning("baseline.var specified but not baseline.val. Please specify 
baseline.val.")
          )
        }
        if (!is.null(baseline.val) & is.null(baseline.var)) {
          return(
            warning("baseline.val specified but not baseline.var. Please specify 
baseline.var.")
          )
        }
        if (is.null(baseline.var) & is.null(baseline.val)) {
          message("time.var specified but not baseline.var and baseline.val")
        }
      }
    }
  }
  z <- list(y = y, design = design, data.type = data.type, long = long,
            columnname = columnname, time.var = time.var, sample.id = sample.id, 
            subject.id = subject.id, baseline.var = baseline.var, 
            baseline.val = baseline.val, control.var = control.var, 
            control.val = control.val)
  return(z)
}
