library(testthat)
library(genBART)

context("cross-sectional data without controls")
test_that("cross-sectional data with no controls outputs correct dendrogram", {
  design.base <- tb.design[which(tb.design$timepoint == 0),]
  expr.base <- tb.expr[1:100,which(tb.design$timepoint == 0)]
  base.means <- apply(expr.base, 1, mean, na.rm = TRUE)
  expr.base.norm <- expr.base - base.means
  h1_test_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(expr.base.norm))))
  h1_test_rowdendro <- as.dendrogram(fastcluster::hclust(dist(expr.base.norm)))
  des.info <- desInfo(y = expr.base, design = design.base, columnname = "columnname")
  dends <- genDendrograms(des.info)
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

context("cross-sectional data with controls")
test_that("cross-sectional data with controls outputs correct dendrograms", {
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
  dends <- genDendrograms(des.info)
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

context("longitudinal data without controls")
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
  dends <- genDendrograms(des.info)
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

context("longitudinal data with controls")
test_that("longitudinal data with controls outputs correct dendrograms", {
  tb.design2 <- tb.design
  tb.design2$clinical_status <- as.character(tb.design2$clinical_status)
  tb.design2$clinical_status[which(tb.design2$clinical_status == "Latent" & tb.design2$timepoint == 0)] <- "HC"
  design.base <- tb.design2[-c(4,19,25,which(tb.design2$timepoint != 0 | tb.design2$clinical_status != "Active")),]
  expr.base <- tb.expr[1:100,-c(4,19,25,which(tb.design2$timepoint != 0 | tb.design2$clinical_status != "Active"))]
  design.hc <- tb.design2[-c(4,19,25,which(tb.design2$clinical_status != "HC")),]
  expr.hc <- tb.expr[1:100,-c(4,19,25,which(tb.design2$clinical_status != "HC"))]
  
  design.all <- tb.design2[-c(4,19,25),]
  expr.all <- tb.expr[1:100,-c(4,19,25)]
  means.all <- apply(expr.all, 1, mean, na.rm = TRUE)
  expr.all.norm <- expr.all - means.all
  h1_test_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(expr.all.norm))))
  h1_test_rowdendro <- as.dendrogram(fastcluster::hclust(dist(expr.all.norm)))
  
  means.hc <- apply(expr.hc, 1, mean, na.rm = TRUE)
  expr.hc.all <- cbind(expr.hc, expr.all[,-which(colnames(expr.all) %in% design.hc$columnname)])
  expr.all.hcnorm <- expr.hc.all - means.hc
  h2_test_coldendro <- as.dendrogram(fastcluster::hclust(dist(t(expr.all.hcnorm))))
  h2_test_rowdendro <- as.dendrogram(fastcluster::hclust(dist(expr.all.hcnorm)))
  
  design.sub <- tb.design2[-c(c(4:6,19:21,25:27), which(tb.design2$clinical_status == "Latent")),]
  expr.sub <- tb.expr[1:100,-c(c(4:6,19:21,25:27), which(tb.design2$clinical_status == "Latent"))]
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
  
  des.info <- desInfo(y = expr.all, design = design.all, 
                      columnname = "columnname",sample_id = "sample_id",patient_id = "monkey_id",
                      long = TRUE,time_var = "timepoint", baseline_var = "timepoint", baseline_val = 0,
                      control_var = "clinical_status", control_val = "HC")
  dends <- genDendrograms(des.info)
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

context("messages")
test_that("correct messages are thrown",{
  suppressWarnings(des.info <- desInfo(y = tb.expr[1:100,], design = tb.design, columnname = "columnname",sample_id = "sample_id",
                      patient_id = "monkey_id",long = TRUE,time_var = "timepoint",control_var = "clinical_status", control_val = "Latent"))
  expect_message(genDendrograms(des.info), "baseline_var and/or baseline_val is unspecified. Cannot produce
baseline median normalized, baseline healthy normalized, or all
samples baseline normalized data.")
  suppressWarnings(des.info <- desInfo(y = tb.expr[1:100,], design = tb.design, columnname = "columnname",sample_id = "sample_id",
                      patient_id = "monkey_id",long = TRUE,time_var = "timepoint"))
  expect_message(genDendrograms(des.info), "baseline_var and/or baseline_val is unspecified. Cannot produce
baseline median normalized or all samples baseline normalized data.")
  des.info <- desInfo(y = tb.expr[1:100,], design = tb.design, columnname = "columnname",sample_id = "sample_id",
                      long = TRUE,time_var = "timepoint",baseline_var = "timepoint",baseline_val = 0,
                      control_var = "clinical_status",control_val = "Latent")
  expect_message(genDendrograms(des.info), "patient_id is not defined. Cannot produce baseline normalized 
data.")
  des.info <- desInfo(y = tb.expr[1:100,], design = tb.design, columnname = "columnname",sample_id = "sample_id",
                      long = TRUE,time_var = "timepoint",baseline_var = "timepoint",baseline_val = 0)
  expect_message(genDendrograms(des.info), "patient_id is not defined. Cannot produce baseline normalized 
data.")
})



  
  
  
  
  
  
  
  
  
  
  
  
  