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
#' @return \code{base_ctrl} Percentage matrix of baseline samples with respect
#'   to controls
#' @return \code{long_base} Percentage matrix of all time point samples with 
#'   respect to their baseline
#' @return \code{long_ctrl} Percentage matrix of all time point samples with 
#'   respect to controls
#' @examples
#' # Example data
#' data(tb.expr)
#' data(tb.design)
#' data(modules)
#' 
#' # Demonstrate on first 100 probes
#' dat <- tb.expr[1:100, ]
#' 
#' # Create desInfo object
#' des.info <- desInfo(y = dat, design = tb.design, data_type = "micro", 
#'                     columnname = "columnname", long = TRUE, patient_id = "monkey_id",
#'                     baseline_var = "timepoint", baseline_val = 0, time_var = "timepoint", 
#'                     responder_var = "clinical_status", sample_id = "sample_id", 
#'                     project_name = "TB")
#' 
#' # Generate module maps                                      
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
    base_ctrl <- NULL
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
      SignMatrix <- final_data
      SignMatrix[is.numeric(final_data)] <- 0
      SignMatrix[final_data < (baseline_mean - 2 * baseline_sd)] <- -1
      SignMatrix[final_data > (baseline_mean + 2 * baseline_sd)] <- 1
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
      long_base <- PercentMatrix + 1
      long_ctrl <- NULL
    }  else {
      long_base <- NULL
      long_ctrl <- NULL
    }
  } else {
    base_index <- design[, baseline_var] == baseline_val
    base_sample_name <- design$columnname[base_index]
    base_sample_grp <- design[, control_var][base_index]
    exp_base_sam <- modexp[, colnames(modexp) %in% base_sample_name]
    grp <- base_sample <- c()
    for (i in 1:ncol(exp_base_sam)) {
      grp[i] <- as.character(base_sample_grp[base_sample_name ==
                                               colnames(exp_base_sam)[i]])
      base_sample[i] <- as.character(design[, sample_id][base_index][
        base_sample_name == colnames(exp_base_sam)[i]])
    }
    exp_base_case <- as.matrix(exp_base_sam[, grp != control_val])
    rownames(exp_base_case) <- modexp$Module
    colnames(exp_base_case) <- base_sample[grp != control_val]
    bhc_mean <- apply(exp_base_sam[, grp == control_val], 1, mean)
    bhc_sd <- apply(exp_base_sam[, grp == control_val], 1, sd)
    SignMatrix <- exp_base_case
    SignMatrix[is.numeric(exp_base_case)] <- 0
    SignMatrix[exp_base_case < (bhc_mean - 2 * bhc_sd)] <- -1
    SignMatrix[exp_base_case > (bhc_mean + 2 * bhc_sd)] <- 1
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
    base_ctrl <- PercentMatrix + 1
    if (long) {
      sort_des <- design[order(design[, time_var], design[, subject_id]), ]
      data_ordered <- modexp[, match(sort_des[, "columnname"], names(modexp),
                                     nomatch = 0)]
      names(data_ordered) <- sort_des[, sample_id]
      exp_case_sample <- as.matrix(data_ordered[, -which(sort_des[[control_var]] 
                                                         == control_val)])
      rownames(exp_case_sample) <- modexp$Module
      SignMatrix <- exp_case_sample
      SignMatrix[is.numeric(exp_case_sample)] <- 0
      SignMatrix[exp_case_sample < (bhc_mean - 2 * bhc_sd)] <- -1
      SignMatrix[exp_case_sample > (bhc_mean + 2 * bhc_sd)] <- 1
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
      long_ctrl <- PercentMatrix + 1
      control_cols <- design$columnname[which(design[[control_var]] ==
                                                       control_val)]
      design2 <- design[-which(design$columnname %in% control_cols), ]
      sort_des <- design2[order(design2[, time_var], design2[, subject_id]), ]
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
      SignMatrix[is.numeric(final_data)] <- 0
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
      long_base <- PercentMatrix + 1
    } else {
      long_base <- NULL
      long_ctrl <- NULL
    }
  }
  z <- list(base_ctrl = base_ctrl, long_base = long_base, long_ctrl = long_ctrl)
  return(z)
}
