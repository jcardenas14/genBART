library(testthat)
library(genBART)

design.base <- tb.design[which(tb.design$timepoint == 0),]
expr.base <- tb.expr[,which(tb.design$timepoint == 0)]

des.info <- desInfo(y = expr.base, design = design.base, columnname = "columnname")
genDendrograms2(des.info)



test_that("baseline data with no controls outputs correct dendrogram", {
  design.base <- tb.design[which(tb.design$timepoint == 0),]
  expr.base <- tb.expr[1:100,which(tb.design$timepoint == 0)]
  base.means <- apply(expr.base, 1, mean, na.rm = TRUE)
  expr.base.norm <- expr.base - base.means
  h1_test_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(expr.base.norm))))
  h1_test_rowdendro <- as.dendrogram(fastcluster::hclust(dist(expr.base.norm)))
  des.info <- desInfo(y = expr.base, design = design.base, columnname = "columnname")
  dends <- genDendrograms2(des.info)
  expect_equal(dends$h1b_coldendro, NULL)
  expect_equal(dends$h1b_rowdendro, NULL)
  expect_equal(dends$h2b_coldendro, NULL)
  expect_equal(dends$h2b_rowdendro, NULL)
  expect_equal(dends$h2_coldendro, NULL)
  expect_equal(dends$h2_rowdendro, NULL)
  expect_equal(dends$h3_coldendro, NULL)
  expect_equal(dends$h3_rowdendro, NULL)
  expect_equal(dends$h1_rowdendro, h1_test_rowdendro)
  expect_equal(dends$h1_coldendro, h1_test_coldendro)
})

test_that("baseline data with controls outputs correct dendrograms", {
  design.base <- tb.design[which(tb.design$timepoint == 0 & tb.design$clinical_status == "Active"),]
  expr.base <- tb.expr[1:100,which(tb.design$timepoint == 0 & tb.design$clinical_status == "Active")]
  design.hc <- tb.design[which(tb.design$timepoint == 0 & tb.design$clinical_status == "Latent"),]
  expr.hc <- tb.expr[1:100,which(tb.design$timepoint == 0 & tb.design$clinical_status == "Latent")]
  design.all <- rbind(design.hc, design.base)
  expr.all <- cbind(expr.hc, expr.base)
  base.means <- apply(expr.base, 1, mean, na.rm = TRUE)
  expr.base.norm <- expr.base - base.means
  h1b_test_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(expr.base.norm))))
  h1b_test_rowdendro <- as.dendrogram(fastcluster::hclust(dist(expr.base.norm)))
  all.means <- apply(expr.all, 1, mean, na.rm = TRUE)
  expr.all.norm <- expr.all - all.means
  h1_test_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(expr.all.norm))))
  h1_test_rowdendro <- as.dendrogram(fastcluster::hclust(dist(expr.all.norm)))
  hc.means <- apply(expr.hc, 1, mean, na.rm = TRUE)
  expr.all.hcnorm <- expr.all - hc.means
  h2_test_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(expr.all.hcnorm))))
  h2_test_rowdendro <- as.dendrogram(fastcluster::hclust(dist(expr.all.hcnorm)))
  des.info <- desInfo(y = expr.all, design = design.all, columnname = "columnname",
                      control_var = "clinical_status", control_val = "Latent")
  dends <- genDendrograms2(des.info)
  expect_equal(dends$h2b_coldendro, NULL)
  expect_equal(dends$h2b_rowdendro, NULL)
  expect_equal(dends$h3_coldendro, NULL)
  expect_equal(dends$h3_rowdendro, NULL)
  expect_equal(dends$h1b_rowdendro, h1b_test_rowdendro)
  expect_equal(dends$h1b_coldendro, h1b_test_coldendro)
  expect_equal(dends$h1_rowdendro, h1_test_rowdendro)
  expect_equal(dends$h1_coldendro, h1_test_coldendro)
  expect_equal(dends$h2_rowdendro, h2_test_rowdendro)
  expect_equal(dends$h2_coldendro, h2_test_coldendro)
})

test_that("longitudinal data with no controls outputs correct dendrograms", {
  design.all <- tb.design[-c(1,13,22),]
  expr.all <- tb.expr[1:100,-c(1,13,22)]
  means.all <- apply(expr.all, 1, mean, na.rm = TRUE)
  expr.all.norm <- expr.all - means.all
  h1_test_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(expr.all.norm))))
  h1_test_rowdendro <- as.dendrogram(fastcluster::hclust(dist(expr.all.norm)))
  
  design.sub <- tb.design[-c(1:3,13:15,22:24),]
  expr.sub <- tb.expr[1:100,-c(1:3,13:15,22:24)]
  design.base <- design.sub[which(design.sub$timepoint == 0),]
  expr.base <- expr.sub[,which(design.sub$timepoint == 0)]
  expr.sub.norm <- expr.sub
  for(i in 1:length(unique(design.sub$monkey_id))){
    index <- which(design.sub$monkey_id %in% unique(design.sub$monkey_id)[i])
    expr.sub.norm[,index] <- expr.sub[,index] - expr.sub[,index[1]]
  }
  expr.sub.norm <- expr.sub.norm[,-which(design.sub$timepoint == 0)]
  h3_test_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(expr.sub.norm))))
  h3_test_rowdendro <- as.dendrogram(fastcluster::hclust(dist(expr.sub.norm)))
  
  means.base <- apply(expr.base, 1, mean, na.rm = TRUE)
  expr.base.norm <- expr.base - means.base
  h1b_test_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(expr.base.norm))))
  h1b_test_rowdendro <- as.dendrogram(fastcluster::hclust(dist(expr.base.norm)))
  
  des.info <- desInfo(y = expr.all, design = design.all, 
                      columnname = "columnname",sample_id = "sample_id",patient_id = "monkey_id",
                      long = TRUE,time_var = "timepoint", baseline_var = "timepoint", baseline_val = 0)
  dends <- genDendrograms2(des.info)
  expect_equal(dends$h2b_coldendro, NULL)
  expect_equal(dends$h2b_rowdendro, NULL)
  expect_equal(dends$h2_coldendro, NULL)
  expect_equal(dends$h2_rowdendro, NULL)
  expect_equal(dends$h1b_rowdendro, h1b_test_rowdendro)
  expect_equal(dends$h1b_coldendro, h1b_test_coldendro)
  expect_equal(dends$h1_rowdendro, h1_test_rowdendro)
  expect_equal(dends$h1_coldendro, h1_test_coldendro)
  expect_equal(dends$h3_rowdendro, h3_test_rowdendro)
  expect_equal(dends$h3_coldendro, h3_test_coldendro)
})

test_that("longitudinal data with controls outputs correct dendrograms", {
  design.base <- tb.design[-c(4,19,25,which(tb.design$timepoint != 0 | tb.design$clinical_status != "Active")),]
  expr.base <- tb.expr[1:100,-c(4,19,25,which(tb.design$timepoint != 0 | tb.design$clinical_status != "Active"))]
  design.hc <- tb.design[-c(4,19,25,which(tb.design$timepoint != 0 | tb.design$clinical_status != "Latent")),]
  expr.hc <- tb.expr[1:100,-c(4,19,25,which(tb.design$timepoint != 0 | tb.design$clinical_status != "Latent"))]
  
  design.all <- tb.design[-c(4,19,25),]
  expr.all <- tb.expr[1:100,-c(4,19,25)]
  means.all <- apply(expr.all, 1, mean, na.rm = TRUE)
  expr.all.norm <- expr.all - means.all
  h1_test_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(expr.all.norm))))
  h1_test_rowdendro <- as.dendrogram(fastcluster::hclust(dist(expr.all.norm)))
  
  means.hc <- apply(expr.hc, 1, mean, na.rm = TRUE)
  expr.hc.all <- cbind(expr.hc, expr.all[,-which(colnames(expr.all) %in% design.hc$columnname)])
  expr.all.hcnorm <- expr.all - means.hc
  h2_test_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(expr.all.hcnorm))))
  h2_test_rowdendro <- as.dendrogram(fastcluster::hclust(dist(expr.all.hcnorm)))
  
  design.sub <- tb.design[-c(c(4:6,19:21,25:27), which(tb.design$clinical_status == "Latent")),]
  expr.sub <- tb.expr[1:100,-c(c(4:6,19:21,25:27), which(tb.design$clinical_status == "Latent"))]
  expr.sub.norm <- expr.sub
  for(i in 1:length(unique(design.sub$monkey_id))){
    index <- which(design.sub$monkey_id %in% unique(design.sub$monkey_id)[i])
    expr.sub.norm[,index] <- expr.sub[,index] - expr.sub[,index[1]]
  }
  expr.sub.norm <- expr.sub.norm[,-which(design.sub$timepoint == 0)]
  h3_test_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(expr.sub.norm))))
  h3_test_rowdendro <- as.dendrogram(fastcluster::hclust(dist(expr.sub.norm)))
  
  means.base <- apply(expr.base, 1, mean, na.rm = TRUE)
  expr.base.norm <- expr.base - means.base
  h1b_test_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(expr.base.norm))))
  h1b_test_rowdendro <- as.dendrogram(fastcluster::hclust(dist(expr.base.norm)))
  
  design.base.hc <- rbind(design.hc, design.base)
  expr.base.hc <- cbind(expr.hc, expr.base)
  expr.base.hcnorm <- expr.base.hc - means.hc
  h2b_test_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(expr.base.hcnorm))))
  h2b_test_rowdendro <- as.dendrogram(fastcluster::hclust(dist(expr.base.hcnorm)))
  
  des.info <- desInfo(y = expr.hc.all, design = design.all, 
                      columnname = "columnname",sample_id = "sample_id",patient_id = "monkey_id",
                      long = TRUE,time_var = "timepoint", baseline_var = "timepoint", baseline_val = 0,
                      control_var = "clinical_status", control_val = "Latent")
  dends <- genDendrograms2(des.info)
  expect_equal(dends$h2b_coldendro, h2b_test_coldendro)
  expect_equal(dends$h2b_rowdendro, h2b_test_rowdendro)
  expect_equal(dends$h2_coldendro, h2_test_coldendro)
  expect_equal(dends$h2_rowdendro, h2_test_rowdendro)
  expect_equal(dends$h1b_rowdendro, h1b_test_rowdendro)
  expect_equal(dends$h1b_coldendro, h1b_test_coldendro)
  expect_equal(dends$h1_rowdendro, h1_test_rowdendro)
  expect_equal(dends$h1_coldendro, h1_test_coldendro)
  expect_equal(dends$h3_rowdendro, h3_test_rowdendro)
  expect_equal(dends$h3_coldendro, h3_test_coldendro)
  
  
})

design_info <- des.info
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
      if (hc) {
        base_sample_id <- design$columnname[which(design[, baseline_var] == 
                                                    baseline_val &
                                                    design[, control_var] != 
                                                    control_val)]
      }
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
        base_sample_id <- design$columnname[which(design[, baseline_var] == 
                                                    baseline_val &
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
      } else {
        message("baseline_var and/or baseline_val is unspecified. Cannot produce
baseline median normalized, baseline healthy normalized, or all
samples baseline normalized data.")
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
      des_wo_controls <- design[-which(design[, control_var] == control_val), ]
      h5index <- c(1, 2, which(colnames(final_expression) %in%
                                 des_wo_controls$columnname))
      y1b <- dataManipulate(y = final_expression[, h5index], 
                            x = des_wo_controls, colname = "columnname",
                            format = "Probes", allsamples = TRUE)
      print("Clustering Non-Healthy Samples Median Normalized Heatmap")
      h1b_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(y1b$heatexp))))
      h1b_rowdendro <- as.dendrogram(fastcluster::hclust(dist(y1b$heatexp)))
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  