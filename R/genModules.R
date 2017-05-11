#' Generate module (gene set) maps for plotting
#'
#' @param design_info list generated from \code{desInfo}
#' @param gene_sets list of gene sets
#' @details This function calculates modules scores whose values are constructed
#'   on individual samples at baseline and for all time points (for longitudinal
#'   studies). For non-longitudinal studies with controls, the control samples
#'   are used to determine an upper and lower threshold (mean HC +/- 2 sd). The
#'   module proportion for each sample is then calculated based on the
#'   percentage of probes within a module that are above or below this
#'   threshold. For example, if 40\% of the probes within a module are above the
#'   threshold and 15\% are below it, then the final module score would be 40\%
#'   - 15\% = 25\% up. For longitudinal studies, additional module scores are
#'   calculated using baseline samples to determine the threshold. For
#'   non-longitudinal studies without controls, \code{genModules} cannot be
#'   used, since there are no reference samples with which to calculate a
#'   threshold.
#' @return \code{base_mod} Percentage matrix of baseline samples with respect to
#'   controls
#' @return \code{long_mod} Percentage matrix of all time point samples with
#'   respect to their baseline
#' @return \code{long_mod2} Percentage matrix of all time point samples with
#'   respect to controls
#' @examples
#' # des.info object is obtained using the desInfo function
#' data(des.info)
#' data(modules)
#' mods <- genModules(des.info, modules)
#' @export
genModules <- function(design_info, gene_sets) {
  y <- design_info$y
  design <- design_info$design
  PROBE_ID <- SYMBOL <- rownames(y)
  final_expression <- cbind(PROBE_ID = PROBE_ID, SYMBOL = SYMBOL, y)
  baseline_var <- design_info$baseline_var
  baseline_val <- design_info$baseline_val
  control_val <- design_info$control_val
  control_var <- design_info$control_var
  sample_id <- design_info$sample_id
  subject_id <- design_info$patient_id
  time_var <- design_info$time_var
  long <- design_info$long
  genes <- unlist(gene_sets)
  gset_names <- rep(names(gene_sets), times = lapply(gene_sets, length))
  geneset2 <- data.frame(SYMBOL = genes, Module = gset_names)
  rownames(geneset2) <- NULL
  modexp <- merge(final_expression, geneset2, by = c("SYMBOL"))
  modordnum <- table(factor(geneset2$Module))
  modnames <- names(modordnum)
  if (is.null(control_var)) {
    base_mod <- NULL
    if (long) {
      sort_des <- design[order(design[, time_var], design[, subject_id]), ]
      data_ordered <- modexp[, match(sort_des[, "columnname"], names(modexp),
                                     nomatch = 0)]
      names(data_ordered) <- sort_des[, sample_id]
      baseline_data <- data_ordered[, which(sort_des[, baseline_var] ==
                                              baseline_val)]
      baseline_mean <- apply(baseline_data, 1, mean)
      baseline_sd <- apply(baseline_data, 1, sd)
      final_data <- as.matrix(data_ordered[, -which(names(data_ordered) %in%
                                                      names(baseline_data))])
      rownames(final_data) <- modexp$Module
      SignMatrix <- matrix(rep(0, dim(final_data)[1] * dim(final_data)[2]),
                           dim(final_data)[1], dim(final_data)[2])
      SignMatrix[final_data < (baseline_mean - 2 * baseline_sd)] <- -1
      SignMatrix[final_data > (baseline_mean + 2 * baseline_sd)] <- 1
      rownames(SignMatrix) <- rownames(final_data)
      colnames(SignMatrix) <- colnames(final_data)
      count_matrix <- c()
      n <- length(modordnum)
      for (i in 1:n) {
        subset <- SignMatrix[which(rownames(SignMatrix) == modnames[i]), ]
        if (is.vector(subset)) {
          count_matrix <- rbind(count_matrix, subset)
        }
        if (is.matrix(subset)) {
          if (dim(subset)[1] == 0) {
            count_matrix <- rbind(count_matrix, rep(0, dim(SignMatrix)[2]))
          }
          if (dim(subset)[1] > 1) {
            count_matrix <- rbind(count_matrix, apply(subset, 2, sum))
          }
        }
      }
      rownames(count_matrix) <- modnames
      PercentMatrix <- count_matrix / as.vector(modordnum)
      long_mod <- PercentMatrix + 1
      long_mod2 <- NULL
    }  else {
      long_mod <- NULL
      long_mod2 <- NULL
    }
  } else {
    base_index <- design[, baseline_var] == baseline_val
    base_sample_name <- design$columnname[base_index]
    base_sample_grp <- design[, control_var][base_index]
    exp_base_sam <- modexp[, colnames(modexp) %in% base_sample_name]
    grp <- base_donorid <- base_sample <- c()
    for (i in 1:ncol(exp_base_sam)) {
      grp[i] <- as.character(base_sample_grp[base_sample_name ==
                                               colnames(exp_base_sam)[i]])
      base_donorid[i] <- as.character(design[, subject_id][base_index][
        base_sample_name == colnames(exp_base_sam)[i]])
      base_sample[i] <- as.character(design[, sample_id][base_index][
        base_sample_name == colnames(exp_base_sam)[i]])
    }
    n_module <- length(modnames)
    n_case <- sum(grp != control_val)
    exp_base_case <- exp_base_sam[, grp != control_val]
    n_base_case <- ncol(exp_base_case)
    bhc_mean <- apply(exp_base_sam[, grp == control_val], 1, mean)
    bhc_sd <- apply(exp_base_sam[, grp == control_val], 1, sd)
    module_l <- c()
    bhc_lower_prop <- bhc_upper_prop <- bhc_sign <- matrix(NA, n_module,
                                                           n_base_case)
    for (i in 1:n_module) {
      index_i <- which(modexp$Module %in% modnames[i])
      module_l[i] <- modordnum[i]
      for (j in 1:n_base_case) {
        bhc_lower_prop[i, j] <- sum(exp_base_case[index_i, j] <
                                      bhc_mean[index_i] - 2 * bhc_sd[index_i]) /
          module_l[i]
        bhc_upper_prop[i, j] <- sum(exp_base_case[index_i, j] >
                                      bhc_mean[index_i] + 2 * bhc_sd[index_i]) /
          module_l[i]
        if (bhc_lower_prop[i, j] > 0.1 & bhc_upper_prop[i, j] < 0.1) {
          bhc_sign[i, j] <- 1 - bhc_lower_prop[i, j]
        } else if (bhc_lower_prop[i, j] < 0.1 & bhc_upper_prop[i, j] > 0.1) {
          bhc_sign[i, j] <- 1 + bhc_upper_prop[i, j]
        } else {
          bhc_sign[i, j] <- 1 + bhc_upper_prop[i, j] - bhc_lower_prop[i, j]
        }
      }
    }
    base_case_sample <- base_sample[grp != control_val]
    colnames(bhc_sign) <- base_case_sample
    rownames(bhc_sign) <- modnames
    base_mod <- bhc_sign
    if (long) {
      exp_sam <- modexp[, colnames(modexp) %in% design$columnname]
      column_name <- colnames(exp_sam)
      grp <- tp <- donorid <- sample_name <- c()
      for (i in 1:ncol(exp_sam)) {
        grp[i] <- as.character(design[, control_var][design$columnname ==
                                                       column_name[i]])
        tp[i] <- as.character(design[, baseline_var][design$columnname ==
                                                       column_name[i]])
        donorid[i] <- as.character(design[, subject_id][design$columnname ==
                                                          column_name[i]])
        sample_name[i] <- as.character(design[, sample_id][design$columnname ==
                                                             column_name[i]])
      }
      exp_case_sample <- exp_sam[, grp != control_val]
      n_case <- ncol(exp_case_sample)
      module_l <- c()
      bhc_lower_prop <- bhc_upper_prop <- bhc_sign <- matrix(NA, n_module,
                                                             n_case)
      for (i in 1:n_module) {
        index_i <- (modexp$Module %in% modnames[i])
        module_l[i] <- modordnum[i]
        for (j in 1:n_case) {
          bhc_lower_prop[i, j] <- sum(exp_case_sample[index_i, j] <
                                      bhc_mean[index_i] - 2 * bhc_sd[index_i]) /
            module_l[i]
          bhc_upper_prop[i, j] <- sum(exp_case_sample[index_i, j] >
                                      bhc_mean[index_i] + 2 * bhc_sd[index_i]) /
            module_l[i]
          if (bhc_lower_prop[i, j] > 0.1 & bhc_upper_prop[i, j] < 0.1) {
            bhc_sign[i, j] <- 1 - bhc_lower_prop[i, j]
          } else if (bhc_lower_prop[i, j] < 0.1 & bhc_upper_prop[i, j] > 0.1) {
            bhc_sign[i, j] <- 1 + bhc_upper_prop[i, j]
          } else {
            bhc_sign[i, j] <- 1 - bhc_lower_prop[i, j] + bhc_upper_prop[i, j]
          }
        }
      }
      case_grp <- grp[grp != control_val]
      case_donor <- donorid[grp != control_val]
      case_wk <- as.numeric(tp[grp != control_val])
      case_sample <- sample_name[grp != control_val]
      case_sample_s <- case_sample[order(case_grp, case_donor, case_wk)]
      bhc_sign_s <- bhc_sign[, order(case_grp, case_donor, case_wk)]
      colnames(bhc_sign_s) <- case_sample_s
      rownames(bhc_sign_s) <- modnames
      long_mod <- bhc_sign_s
      control_cols <- design$columnname[which(design[, control_var] ==
                                                control_val)]
      final_expression2 <- final_expression[, -which(colnames(final_expression)
                                                     %in% control_cols)]
      design2 <- design[-which(design$columnname %in% control_cols),
                        ]
      sort_des <- design2[order(design2[, time_var], design2[, subject_id]),
                          ]
      data_ordered <- modexp[, match(sort_des[, "columnname"], names(modexp),
                                     nomatch = 0)]
      names(data_ordered) <- sort_des[, sample_id]
      baseline_data <- data_ordered[, which(sort_des[, baseline_var] ==
                                              baseline_val)]
      baseline_mean <- apply(baseline_data, 1, mean)
      baseline_sd <- apply(baseline_data, 1, sd)
      final_data <- as.matrix(data_ordered[, -which(names(data_ordered) %in%
                                                      names(baseline_data))])
      rownames(final_data) <- modexp$Module
      SignMatrix <- final_data
      SignMatrix[final_data < (baseline_mean - 2 * baseline_sd)] <- -1
      SignMatrix[final_data > (baseline_mean + 2 * baseline_sd)] <- 1
      count_matrix <- c()
      for (i in 1:length(modordnum)) {
        subset <- SignMatrix[which(rownames(SignMatrix) == modnames[i]), ]
        if (is.vector(subset)) {
          count_matrix <- rbind(count_matrix, subset)
        }
        if (is.matrix(subset)) {
          if (dim(subset)[1] == 0) {
            count_matrix <- rbind(count_matrix, rep(0, dim(SignMatrix)[2]))
          }
          if (dim(subset)[1] > 1) {
            count_matrix <- rbind(count_matrix, apply(subset, 2, sum))
          }
        }
      }
      rownames(count_matrix) <- modnames
      PercentMatrix <- count_matrix / as.vector(modordnum)
      long_mod2 <- PercentMatrix + 1
    } else {
      long_mod <- NULL
      long_mod2 <- NULL
    }
  }
  z <- list(base_mod = base_mod, long_mod = long_mod, long_mod2 = long_mod2)
  return(z)
}
