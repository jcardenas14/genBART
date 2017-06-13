#' Declare design information for downstream analysis
#' 
#' Match design and expression matrices. Declare and store design parameters.
#' @param y An expression matrix or dataframe
#' @param design A design matrix or dataframe containing sample information 
#'   (e.g. age, condition, timepoint, etc.)
#' @param data_type What type of data is being analyzed? ("micro", "rna", 
#'   "flow", "metab"). Default is data_type="micro"
#' @param columnname Name of column in design that contains the column names of 
#'   y
#' @param long logical; Is the study longitudinal?
#' @param patient_id Name of column in design that contains patient ids
#' @param baseline_var Name of column in design that contains values referring 
#'   to baseline observations
#' @param baseline_val string or numeric value denoting baseline observations
#' @param control_var Name of column in design that contains values referring to
#'   controls
#' @param control_val string or numeric value denoting controls
#' @param time_var For longitudinal studies. Name of column in design that 
#'   contains the study timepoints
#' @param responder_var Name of column in design that contains subject responder
#'   status or condition
#' @param summary_var Column name in design denoting a numeric variable over 
#'   which sample design summary tables are initialized (e.g. Age)
#' @param sample_id Name of column in design that contains unique sample 
#'   identification
#' @param project_name String denoting the name of the project or study
#' @return A list containing matched design and expression and all of the design
#'   parameters specified
#' @examples
#' # Using example data
#' data(tb.expr)
#' data(tb.design)
#' des.info <- desInfo(y = tb.expr, design = tb.design, data_type = "micro", 
#'                     columnname = "columnname", long = TRUE, sample_id = "sample_id",
#'                     patient_id = "monkey_id", time_var = "timepoint",
#'                     baseline_var = "timepoint", baseline_val = 0,
#'                     responder_var = "clinical_status", project_name = "TB")
#' @export
desInfo <- function(y, design, data_type = "micro", long = FALSE, 
                    columnname = NULL, sample_id = NULL, patient_id = NULL, 
                    time_var = NULL, baseline_var = NULL, baseline_val = NULL, 
                    control_var = NULL, control_val = NULL, summary_var = NULL, 
                    responder_var = NULL, project_name = NULL) {
  if (is.null(columnname)) {
    return(warning("Must enter columnname to match design and expression"))
  } else {
    if (!any(design[[columnname]] %in% colnames(y))) {
      return(
        warning("columnname values in design do not match any columnnames in y")
      )
    }
    if (!is.null(sample_id)) {
      if (!is.null(patient_id)) {
        if (columnname == sample_id) {
          design$columnname <- as.character(design[[sample_id]])
        }
        if (columnname == patient_id) {
          design$columnname <- as.character(design[[patient_id]])
        }
        if (columnname != sample_id & columnname != patient_id) {
          colnames(design)[colnames(design) %in% columnname] <- "columnname"
        }
      } else {
        if (columnname == sample_id) {
          design$columnname <- as.character(design[[sample_id]])
        } else {
          colnames(design)[colnames(design) %in% columnname] <- "columnname"
        }
      }
    } else {
      if (!is.null(patient_id)) {
        if (columnname == patient_id) {
          design$columnname <- as.character(design[[patient_id]])
        } else {
          colnames(design)[colnames(design) %in% columnname] <- "columnname"
        }
      } else {
        colnames(design)[colnames(design) %in% columnname] <- "columnname"
      }
    }
  }
  if (nrow(design) > ncol(y)) {
    message("More samples in design than y. Throwing out excess samples.")
    design <- design[match(colnames(y), design$columnname, nomatch = 0), ]
  }
  if (nrow(design) < ncol(y)) {
    message("More samples in y than design. Throwing out excess samples.")
    y <- y[, match(design$columnname, colnames(y), nomatch = 0)]
  } else {
    y <- y[, match(design$columnname, colnames(y), nomatch = 0)]
  }
  if (is.null(control_var) || is.null(control_val)) {
    if (!is.null(control_var) & is.null(control_val)) {
      warning("control_var specified without control_val")
    }
    if (!is.null(control_val) & is.null(control_var)) {
      warning("control_val specified without control_var")
    }
    hc <- FALSE
  } else {
    hc <- TRUE
  }
  if (long) {
    if (is.null(time_var)) {
      warning("long is TRUE but time_var not specified")
    }
    if (is.null(baseline_var)) {
      warning("long is TRUE but baseline_var not specified")
    }
    if (is.null(baseline_val)) {
      warning("long is TRUE but baseline_val not specified")
    }
  } else {
    if (!is.null(time_var)) {
      message("time_var specified. Setting long to TRUE.")
      if (is.null(baseline_var) || is.null(baseline_val)) {
        if (!is.null(baseline_var) & is.null(baseline_val)) {
          warning("baseline_var specified without baseline_val")
        }
        if (!is.null(baseline_val) & is.null(baseline_var)) {
          warning("baseline_val specified without baseline_var")
        }
        if (is.null(baseline_var) & is.null(baseline_val)) {
          message("time_var specified without baseline_var and baseline_val")
        }
      }
    }
  }
  z <- list(y = y, data_type = data_type, design = design, long = long,
            columnname = columnname, sample_id = sample_id, patient_id = patient_id,
            time_var = time_var, baseline_var = baseline_var, 
            baseline_val = baseline_val, control_var = control_var,
            control_val = control_val, summary_var = summary_var, 
            responder_var = responder_var, hc = hc, project_name = project_name)
  return(z)
}
