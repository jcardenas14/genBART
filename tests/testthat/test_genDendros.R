library(testthat)
library(genBART)

test_that("baseline data only produces one mean normalized dendrogram", {
  test.design <- tb.design
  test.expr <- tb.expr
  
})


genDendrograms2 <- function(design_info) {
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
  print("Normalizing Expression and Clustering....")
  y1 <- dataManipulate(y = final_expression, x = design, 
                       colname = "columnname", format = "Probes",
                       allsamples = TRUE)
  print("Clustering All Samples Median Normalized Heatmap")
  h1_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(y1$heatexp))))
  h1_rowdendro <- as.dendrogram(fastcluster::hclust(dist(y1$heatexp)))
  if (long) {
    if (hc) {
      y2 <- dataManipulate(y = final_expression, x = design,
                           colname = "columnname", ref_var = control_var,
                           ref_level = control_val, long = FALSE, keep = TRUE,
                           format = "Probes", allsamples = FALSE)
      print("Clustering All Samples Healthy Normalized Heatmap")
      h2_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(y2$heatexp))))
      h2_rowdendro <- as.dendrogram(fastcluster::hclust(dist(y2$heatexp)))
    } else {
      h2_coldendro <- NULL
      h2_rowdendro <- NULL
    }
    if (!is.null(baseline_var) & !is.null(baseline_val)) {
      base_sample_id <- design$columnname[design[, baseline_var] == baseline_val
                                          ]
      ind1 <- which(colnames(final_expression) %in% c("PROBE_ID", "SYMBOL"))
      ind2 <- which(colnames(final_expression) %in% base_sample_id)
      exp_base <- final_expression[, c(ind1, ind2)]
      des_base <- design[which(design$columnname %in% colnames(exp_base)), ]
      y1b <- dataManipulate(y = exp_base, x = des_base, colname = "columnname",
                            format = "Probes", allsamples = TRUE)
      print("Clustering Baseline Median Normalized Heatmap")
      h1b_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(y1b$heatexp))))
      h1b_rowdendro <- as.dendrogram(fastcluster::hclust(dist(y1b$heatexp)))
      if (hc) {
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
        y3 <- dataManipulate(y = final_expression[, h5index], 
                             x = des_wo_controls, colname = "columnname", 
                             ref_var = baseline_var, ref_level = baseline_val, 
                             long = TRUE, subjects = patient_id, keep = FALSE,
                             format = "Probes", allsamples = FALSE)
        print("Clustering All Samples Baseline Normalized Heatmap")
        h3_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(y3$heatexp))))
        h3_rowdendro <- as.dendrogram(fastcluster::hclust(dist(y3$heatexp)))
      } else {
        h2b_coldendro <- NULL
        h2b_rowdendro <- NULL
        message("baseline_var and/or baseline_val is unspecified. Cannot produce
baseline median normalized, baseline healthy normalized, or all
samples baseline normalized data.")
        y3 <- dataManipulate(y = final_expression, x = design,
                             colname = "columnname", ref_var = baseline_var,
                             ref_level = baseline_val, long = TRUE,
                             subjects = patient_id, keep = FALSE,
                             format = "Probes", allsamples = FALSE)
        print("Clustering All Samples Baseline Normalized Heatmap")
        h3_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(y3$heatexp))))
        h3_rowdendro <- as.dendrogram(fastcluster::hclust(dist(y3$heatexp)))
      }
    } else {
      h1b_coldendro <- NULL
      h1b_rowdendro <- NULL
      h3_coldendro <- NULL
      h3_rowdendro <- NULL
      if (!hc) {
        message("baseline_var and/or baseline_val is unspecified. Cannot produce
baseline median normalized or all samples baseline normalized data.")
      }
    }
  } else {
    h1b_coldendro <- NULL
    h1b_rowdendro <- NULL
    h2b_coldendro <- NULL
    h2b_rowdendro <- NULL
    h3_coldendro <- NULL
    h3_rowdendro <- NULL
    if (hc) {
      y2 <- dataManipulate(y = final_expression, x = design,
                           colname = "columnname", ref_var = control_var,
                           ref_level = control_val, long = FALSE, keep = TRUE,
                           format = "Probes", allsamples = FALSE)
      print("Clustering All Samples Healthy Normalized Heatmap")
      h2_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(y2$heatexp))))
      h2_rowdendro <- as.dendrogram(fastcluster::hclust(dist(y2$heatexp)))
    } else {
      h2_coldendro <- NULL
      h2_rowdendro <- NULL
    }
  }
  z <- list(h1b_coldendro = h1b_coldendro, h1b_rowdendro = h1b_rowdendro,
            h2b_coldendro = h2b_coldendro, h2b_rowdendro = h2b_rowdendro,
            h1_coldendro = h1_coldendro, h1_rowdendro = h1_rowdendro,
            h2_coldendro = h2_coldendro, h2_rowdendro = h2_rowdendro,
            h3_coldendro = h3_coldendro, h3_rowdendro = h3_rowdendro)
  return(z)
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  