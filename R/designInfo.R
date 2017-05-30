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
#'                     columnname = "columnname", long = TRUE, patient_id = "monkey_id",
#'                     baseline_var = "timepoint", baseline_val = 0, time_var = "timepoint", 
#'                     responder_var = "clinical_status", sample_id = "sample_id", 
#'                     project_name = "TB")
#' @export
desInfo <- function(y, design, data_type = "micro", long = FALSE,
                    columnname = NULL, patient_id = NULL, baseline_var = NULL,
                    baseline_val = NULL, control_var = NULL, control_val = NULL,
                    time_var = NULL, responder_var = NULL, summary_var = NULL,
                    sample_id = NULL, project_name = NULL) {
  if (is.null(control_var) || is.null(control_val)) {
    hc <- FALSE
  } else {
    hc <- TRUE
  }
  if (is.null(columnname)) {
    return("Must enter columnname to match design and expression")
  } else {
    if (!is.null(sample_id)) {
      if (columnname == sample_id) {
        design$columnname <- as.character(design[[sample_id]])
      }
    }
    if (!is.null(patient_id)) {
      if (columnname == patient_id) {
        design$columnname <- as.character(design[[patient_id]])
      }
    }
    if (columnname != sample_id & columnname != patient_id) {
      colnames(design)[which(colnames(design) %in% columnname)] <- "columnname"
    }
  }
  if (nrow(design) > ncol(y)) {
    print("More samples in design than expression. Throwing out excess samples
          from design and matching.")
    design <- design[match(colnames(y), design$columnname, nomatch = 0), ]
  }
  if (nrow(design) < ncol(y)) {
    print("More samples in expression than design. Throwing out excess samples
          from expression and matching.")
    y <- y[, match(design$columnname, colnames(y), nomatch = 0)]
  } else {
    y <- y[, match(design$columnname, colnames(y), nomatch = 0)]
  }
  z <- list(y = y, data_type = data_type, design = design,
            columnname = columnname, long = long, patient_id = patient_id,
            sample_id = sample_id, baseline_var = baseline_var,
            baseline_val = baseline_val, control_var = control_var,
            control_val = control_val, time_var = time_var,
            responder_var = responder_var, summary_var = summary_var, hc = hc,
            project_name = project_name)
  return(z)
}
