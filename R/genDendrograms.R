getNorm <- function(x, y, colname, id, mynames, index_sid, index_refvar,
                    ref_level, keep = TRUE) {
  dat <- x[which(x[, index_sid] == id), mynames]
  lev_name <- dat[dat[, 3] == ref_level, colname]
  if (keep == TRUE) {
    lev_norm <- y[, as.character(dat[, colname])] - y[, as.character(lev_name)]
  }
  if (keep == FALSE) {
    index <- setdiff(as.character(dat[, colname]), as.character(lev_name))
    lev_norm <- y[, index] - y[, as.character(lev_name)]
    if (length(index) == 1) {
      lev_norm <- data.frame(lev_norm)
      names(lev_norm) <- as.character(dat[, colname][dat[, colname] %in% index])
    }
  }
  return(lev_norm)
}

getNorm2 <- function(y, colnames_lev, colnames_lev_rm, keep = TRUE) {
  index_lev <- which(names(y) %in% colnames_lev)
  index_lev_rm <- which(names(y) %in% colnames_lev_rm)
  if (keep == FALSE) {
    normdat <- y[, index_lev_rm] - apply(y[, index_lev], 1, mean)
  }
  if (keep == TRUE) {
    normdat <- y[, c(index_lev, index_lev_rm)] - apply(y[, index_lev], 1, mean)
  }
  return(normdat)
}

dataManipulate <- function(y, x, colname, ref_var, ref_level, long = FALSE,
                           subjects, keep = TRUE, format = "Probes",
                           allsamples = TRUE) {
  index_samples <- which(names(y) %in% x[, colname])
  if (format == "Probes") {
    index_ps <- which(names(y) %in% c("PROBE_ID", "SYMBOL"))
  }
  if (format == "Modules") {
    index_ps <- which(names(y) %in% c("Module"))
  }
  y <- y[, c(index_ps, index_samples)]
  if (allsamples == FALSE) {
    if (long == TRUE) {
      index_refvar <- which(names(x) == ref_var)
      lev <- x[which(x[, index_refvar] == ref_level), ]
      lev_rm <- x[which(x[, index_refvar] != ref_level), ]
      index_subid <- which(names(x) == subjects)
      lev_subids <- unique(lev[, index_subid])
      lev_rm_subids <- unique(lev_rm[, index_subid])
      index_include <- lev_subids %in% lev_rm_subids
      complete_ids <- lev_subids[index_include]
      lev_norm <- vector("list", length(complete_ids))
      for (i in 1:(length(lev_norm))) {
        lev_norm[[i]] <- getNorm(y = y, x = x, colname = colname,
                                 id = complete_ids[i], index_sid = index_subid,
                                 mynames = c(subjects, "columnname", ref_var),
                                 index_refvar = index_refvar,
                                 ref_level = ref_level, keep = keep)
      }
      lev_norm <- do.call("cbind", lev_norm)
      design_norm <- x[(x[, index_subid] %in% complete_ids), ]
      if (keep == FALSE) {
        design_norm <- design_norm[(design_norm[, index_refvar] != ref_level), ]
      }
    }
    if (long == FALSE) {
      index_refvar <- which(names(x) == ref_var)
      lev <- x[which(x[, index_refvar] == ref_level), ]
      lev_rm <- x[which(x[, index_refvar] != ref_level), ]
      lev_names <- lev[, colname]
      lev_rm_names <- lev_rm[, colname]
      lev_norm <- getNorm2(y = y, lev_names, lev_rm_names, keep = keep)
      if (keep == FALSE) {
        design_norm <- lev_rm
      }
      if (keep == TRUE) {
        design_norm <- x
      }
    }
  } else {
    means_y <- as.vector(apply(y[, -(1:2)], 1, mean))
    lev_norm <- y[, -c(1:2)] - means_y
    design_norm <- x
  }
  final_norm_shuff <- as.matrix(lev_norm)
  colnames(final_norm_shuff) <- names(lev_norm)
  if (format == "Probes") {
    rownames(final_norm_shuff) <- y$SYMBOL
  }
  if (format == "Modules") {
    rownames(final_norm_shuff) <- y$Module
  }
  return(list(heatexp = final_norm_shuff, heatdes = design_norm))
}

#' Normalization and hierarchical clustering
#' 
#' Performs hierarchical clustering and various normalizations of expression 
#' matrix
#' @param design_info list generated from \code{desInfo}
#' @details This function performs various normalizations of the expression 
#'   data, depending on the study design and the parameters defined in 
#'   \code{\link{desInfo}}. For all study designs, the data is normalized to the
#'   mean of all the samples. For non-logitudinal studies with controls, an 
#'   additional normalization to the mean of the controls is performed. For 
#'   longitudinal studies, baseline normalization is performed (each subject's 
#'   baseline is subtracted out from later time points), along with separate 
#'   normalization of baseline samples only (proceeding as with non-longitudinal
#'   studies). For each normalized expression matrix, hierarchical clustering of
#'   the rows is performed.
#' @return Returns dendrogram objects obtained through hierarchical clustering 
#'   of the various normalized expression datasets.
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
#'                     columnname = "columnname", long = TRUE, sample_id = "sample_id",
#'                     patient_id = "monkey_id", time_var = "timepoint",
#'                     baseline_var = "timepoint", baseline_val = 0, 
#'                     responder_var = "clinical_status", project_name = "TB")
#' 
#' # Normalize and cluster data
#' dendros <- genDendrograms(des.info)
#' @importFrom stats as.dendrogram dist qt sd
#' @export
genDendrograms <- function(design_info) {
  hc <- design_info$hc
  long <- design_info$long
  patient_id <- design_info$patient_id
  baseline_var <- design_info$baseline_var
  baseline_val <- design_info$baseline_val
  control_var <- design_info$control_var
  control_val <- design_info$control_val
  design <- design_info$design
  final_expression <- design_info$y
  PROBE_ID <- SYMBOL <- rownames(final_expression)
  final_expression <- cbind(PROBE_ID = PROBE_ID, SYMBOL = SYMBOL,
                            final_expression)
  h1_coldendro <- h1_rowdendro <- NULL
  h2_coldendro <- h2_rowdendro <- NULL
  h1b_coldendro <- h1b_rowdendro <- NULL
  h2b_coldendro <- h2b_rowdendro <- NULL
  h3_coldendro <- h3_rowdendro <- NULL
  print("Normalizing Expression and Clustering....")
  y1 <- dataManipulate(y = final_expression, x = design, 
                       colname = "columnname", format = "Probes",
                       allsamples = TRUE)
  print("Clustering All Samples Median Normalized Heatmap")
  h1_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(y1$heatexp))))
  h1_rowdendro <- as.dendrogram(fastcluster::hclust(dist(y1$heatexp)))
  if (hc) {
    y2 <- dataManipulate(y = final_expression, x = design, 
                         colname = "columnname", ref_var = control_var,
                         ref_level = control_val, long = FALSE, keep = TRUE,
                         format = "Probes", allsamples = FALSE)
    print("Clustering All Samples Healthy Normalized Heatmap")
    h2_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(y2$heatexp))))
    h2_rowdendro <- as.dendrogram(fastcluster::hclust(dist(y2$heatexp)))
    if (long) {
      if (!is.null(baseline_var) & !is.null(baseline_val)) {
        base_sample_id <- design$columnname[which(design[, baseline_var] == 
                                                    baseline_val &
                                                    design[, control_var] != 
                                                    control_val)]
        ind1 <- which(colnames(final_expression) %in% c("PROBE_ID", "SYMBOL"))
        ind2 <- which(colnames(final_expression) %in% base_sample_id)
        exp_base <- final_expression[, c(ind1, ind2)]
        des_base <- design[which(design$columnname %in% colnames(exp_base)), ]
        y1b <- dataManipulate(y = exp_base, x = des_base, 
                              colname = "columnname", format = "Probes", 
                              allsamples = TRUE)
        print("Clustering Baseline Median Normalized Heatmap")
        h1b_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(y1b$heatexp)))
        )
        h1b_rowdendro <- as.dendrogram(fastcluster::hclust(dist(y1b$heatexp)))
        base_sample_id <- design$columnname[which(design[, baseline_var] == 
                                                    baseline_val |
                                                    design[, control_var] == 
                                                    control_val)]
        ind1 <- which(colnames(final_expression) %in% c("PROBE_ID", "SYMBOL"))
        ind2 <- which(colnames(final_expression) %in% base_sample_id)
        exp_base <- final_expression[, c(ind1, ind2)]
        des_base <- design[which(design$columnname %in% colnames(exp_base)), ]
        y2b <- dataManipulate(y = exp_base, x = des_base,
                              colname = "columnname", ref_var = control_var, 
                              ref_level = control_val, long = FALSE, 
                              keep = TRUE, format = "Probes", allsamples = FALSE
                              )
        print("Clustering Baseline Healthy Normalized Heatmap")
        h2b_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(y2b$heatexp)))
                                       )
        h2b_rowdendro <- as.dendrogram(fastcluster::hclust(dist(y2b$heatexp)))
        des_wo_controls <- design[-which(design[, control_var] == control_val), 
                                  ]
        h5index <- c(1, 2, which(colnames(final_expression) %in%
                                   des_wo_controls$columnname))
        if (is.null(patient_id)) {
          message("patient_id is not defined. Cannot produce baseline normalized 
data.")
        } else {
          y3 <- dataManipulate(y = final_expression[, h5index], 
                               x = des_wo_controls, colname = "columnname", 
                               ref_var = baseline_var, ref_level = baseline_val, 
                               long = TRUE, subjects = patient_id, keep = FALSE,
                               format = "Probes", allsamples = FALSE)
          print("Clustering All Samples Baseline Normalized Heatmap")
          h3_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(y3$heatexp)))
                                        )
          h3_rowdendro <- as.dendrogram(fastcluster::hclust(dist(y3$heatexp)))
        }
      } else {
        message("baseline_var and/or baseline_val is unspecified. Cannot produce
baseline median normalized, baseline healthy normalized, or all
samples baseline normalized data.")
      }
    } else {
      des_wo_controls <- design[-which(design[, control_var] == control_val), ]
      h5index <- c(1, 2, which(colnames(final_expression) %in% 
                                  des_wo_controls$columnname))
      y1b <- dataManipulate(y = final_expression[, h5index], 
                            x = des_wo_controls, colname = "columnname",
                            format = "Probes", allsamples = TRUE)
      print("Clustering Non-Healthy Samples Median Normalized Heatmap")
      h1b_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(y1b$heatexp))))
      h1b_rowdendro <- as.dendrogram(fastcluster::hclust(dist(y1b$heatexp)))
    }
  }
  if (!hc) {
    if (long) {
      if (!is.null(baseline_var) & !is.null(baseline_val)) {
        base_sample_id <- design$columnname[design[, baseline_var] == 
                                              baseline_val]
        ind1 <- which(colnames(final_expression) %in% c("PROBE_ID", "SYMBOL"))
        ind2 <- which(colnames(final_expression) %in% base_sample_id)
        exp_base <- final_expression[, c(ind1, ind2)]
        des_base <- design[which(design$columnname %in% colnames(exp_base)), ]
        y1b <- dataManipulate(y = exp_base, x = des_base, 
                              colname = "columnname", format = "Probes", 
                              allsamples = TRUE)
        print("Clustering Baseline Median Normalized Heatmap")
        h1b_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(y1b$heatexp)))
        )
        h1b_rowdendro <- as.dendrogram(fastcluster::hclust(dist(y1b$heatexp)))
        if (is.null(patient_id)) {
          message("patient_id is not defined. Cannot produce baseline normalized 
data.")
        } else {
          y3 <- dataManipulate(y = final_expression, x = design,
                               colname = "columnname", ref_var = baseline_var,
                               ref_level = baseline_val, long = TRUE,
                               subjects = patient_id, keep = FALSE,
                               format = "Probes", allsamples = FALSE)
          print("Clustering All Samples Baseline Normalized Heatmap")
          h3_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(y3$heatexp)))
          )
          h3_rowdendro <- as.dendrogram(fastcluster::hclust(dist(y3$heatexp)))
        }
      } else {
        message("baseline_var and/or baseline_val is unspecified. Cannot produce
baseline median normalized or all samples baseline normalized data.")
      }
    }
  }
  z <- list(h1b_coldendro = h1b_coldendro, h1b_rowdendro = h1b_rowdendro,
            h2b_coldendro = h2b_coldendro, h2b_rowdendro = h2b_rowdendro,
            h1_coldendro = h1_coldendro, h1_rowdendro = h1_rowdendro,
            h2_coldendro = h2_coldendro, h2_rowdendro = h2_rowdendro,
            h3_coldendro = h3_coldendro, h3_rowdendro = h3_rowdendro)
  return(z)
}
