library(testthat)
library(genBart)

context("cross-sectional data without controls")
test_that("cross-sectional data with no controls outputs correctly normalized data", {
  design.base <- tb.design[which(tb.design$timepoint == 0),]
  expr.base <- tb.expr[1:100,which(tb.design$timepoint == 0)]
  base.means <- apply(expr.base, 1, mean, na.rm = TRUE)
  expr.base.norm <- expr.base - base.means
  rowdend1 <- as.dendrogram(fastcluster::hclust(dist(expr.base.norm)))
  meta <- metaData(y = expr.base, design = design.base, columnname = "columnname")
  normData <- normalizeData(meta)
  clustData <- clusterData(normData)
  expect_equal(normData$y1, expr.base.norm)
  expect_equal(normData$y1b, NULL)
  expect_equal(normData$y2b, NULL)
  expect_equal(normData$y2, NULL)
  expect_equal(normData$norm.method, "mean")
  expect_equal(clustData$rowdend1, rowdend1)
  expect_equal(clustData$rowdend2, NULL)
  expect_equal(clustData$rowdend1b, NULL)
  expect_equal(clustData$rowdend2b, NULL)
  expect_equal(clustData$rowdend3, NULL)
  expect_equal(clustData$norm.method, "mean")
  expect_equal(clustData$dist.method, "euclidean")
  expect_equal(clustData$agg.method, "complete")
})

context("cross-sectional data with controls")
test_that("cross-sectional data with controls outputs correct dendrograms", {
  design.base <- tb.design[which(tb.design$timepoint == 0 & tb.design$clinical_status == "Active"),]
  expr.base <- tb.expr[1:100,which(tb.design$timepoint == 0 & tb.design$clinical_status == "Active")]
  design.hc <- tb.design[which(tb.design$timepoint == 0 & tb.design$clinical_status == "Latent"),]
  expr.hc <- tb.expr[1:100,which(tb.design$timepoint == 0 & tb.design$clinical_status == "Latent")]
  design.all <- rbind(design.hc, design.base)
  expr.all <- cbind(expr.hc, expr.base)
  all.means <- apply(expr.all, 1, mean, na.rm = TRUE)
  expr.all.norm <- expr.all - all.means
  rowdend1 <- as.dendrogram(fastcluster::hclust(dist(expr.all.norm)))
  hc.means <- apply(expr.hc, 1, mean, na.rm = TRUE)
  expr.all.hcnorm <- expr.all - hc.means
  rowdend2 <- as.dendrogram(fastcluster::hclust(dist(expr.all.hcnorm)))
  meta <- metaData(y = expr.all, design = design.all, columnname = "columnname",
                      control.var = "clinical_status", control.val = "Latent")
  normData <- normalizeData(meta)
  clustData <- clusterData(normData)
  expect_equal(normData$y1, expr.all.norm)
  expect_equal(normData$y2, expr.all.hcnorm)
  expect_equal(normData$y1b, NULL)
  expect_equal(normData$y2b, NULL)
  expect_equal(normData$y3, NULL)
  expect_equal(normData$norm.method, "mean")
  expect_equal(clustData$rowdend1, rowdend1)
  expect_equal(clustData$rowdend2, rowdend2)
  expect_equal(clustData$rowdend1b, NULL)
  expect_equal(clustData$rowdend2b, NULL)
  expect_equal(clustData$rowdend3, NULL)
  expect_equal(clustData$norm.method, "mean")
  expect_equal(clustData$dist.method, "euclidean")
  expect_equal(clustData$agg.method, "complete")
})

context("longitudinal data without controls")
test_that("longitudinal data with no controls outputs correct dendrograms", {
  design.all <- tb.design[-c(1,13,22),]
  expr.all <- tb.expr[1:100,-c(1,13,22)]
  means.all <- apply(expr.all, 1, mean, na.rm = TRUE)
  expr.all.norm <- expr.all - means.all
  rowdend1 <- as.dendrogram(fastcluster::hclust(dist(expr.all.norm)))
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
  rowdend3 <- as.dendrogram(fastcluster::hclust(dist(expr.sub.norm)))
  means.base <- apply(expr.base, 1, mean, na.rm = TRUE)
  expr.base.norm <- expr.base - means.base
  rowdend1b <- as.dendrogram(fastcluster::hclust(dist(expr.base.norm)))
  meta <- metaData(y = expr.all, design = design.all, 
                      columnname = "columnname",sample.id = "sample_id",subject.id = "monkey_id",
                      long = TRUE,time.var = "timepoint", baseline.var = "timepoint", baseline.val = 0)
  normData <- normalizeData(meta)
  clustData <- clusterData(normData)
  expect_equal(normData$y1, expr.all.norm)
  expect_equal(normData$y2, NULL)
  expect_equal(normData$y1b, expr.base.norm)
  expect_equal(normData$y2b, NULL)
  expect_equal(normData$y3, expr.sub.norm)
  expect_equal(normData$norm.method, "mean")
  expect_equal(clustData$rowdend1, rowdend1)
  expect_equal(clustData$rowdend2, NULL)
  expect_equal(clustData$rowdend1b, rowdend1b)
  expect_equal(clustData$rowdend2b, NULL)
  expect_equal(clustData$rowdend3, rowdend3)
  expect_equal(clustData$norm.method, "mean")
  expect_equal(clustData$dist.method, "euclidean")
  expect_equal(clustData$agg.method, "complete")
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
  rowdend1 <- as.dendrogram(fastcluster::hclust(dist(expr.all.norm)))
  means.hc <- apply(expr.hc, 1, mean, na.rm = TRUE)
  expr.all.hcnorm <- expr.all - means.hc
  rowdend2 <- as.dendrogram(fastcluster::hclust(dist(expr.all.hcnorm)))
  design.sub <- tb.design2[-c(c(4:6,19:21,25:27), which(tb.design2$clinical_status == "Latent")),]
  expr.sub <- tb.expr[1:100,-c(c(4:6,19:21,25:27), which(tb.design2$clinical_status == "Latent"))]
  expr.sub.norm <- expr.sub
  for(i in 1:length(unique(design.sub$monkey_id))){
    index <- which(design.sub$monkey_id %in% unique(design.sub$monkey_id)[i])
    expr.sub.norm[,index] <- expr.sub[,index] - expr.sub[,index[1]]
  }
  expr.sub.norm <- expr.sub.norm[,-which(design.sub$timepoint == 0)]
  rowdend3 <- as.dendrogram(fastcluster::hclust(dist(expr.sub.norm)))
  means.base <- apply(expr.base, 1, mean, na.rm = TRUE)
  expr.base.norm <- expr.base - means.base
  rowdend1b <- as.dendrogram(fastcluster::hclust(dist(expr.base.norm)))
  expr.base.hc <- expr.all[,which(design.all$timepoint == 0)]
  expr.base.hcnorm <- expr.base.hc - means.hc
  rowdend2b <- as.dendrogram(fastcluster::hclust(dist(expr.base.hcnorm)))
  meta <- metaData(y = expr.all, design = design.all, columnname = "columnname",
                   sample.id = "sample_id",subject.id = "monkey_id",long = TRUE,
                   time.var = "timepoint", baseline.var = "timepoint", baseline.val = 0,
                   control.var = "clinical_status", control.val = "HC")
  normData <- normalizeData(meta)
  clustData <- clusterData(normData)
  expect_equal(normData$y1, expr.all.norm)
  expect_equal(normData$y2, expr.all.hcnorm)
  expect_equal(normData$y1b, expr.base.norm)
  expect_equal(normData$y2b, expr.base.hcnorm)
  expect_equal(normData$y3, expr.sub.norm)
  expect_equal(normData$norm.method, "mean")
  expect_equal(clustData$rowdend1, rowdend1)
  expect_equal(clustData$rowdend2, rowdend2)
  expect_equal(clustData$rowdend1b, rowdend1b)
  expect_equal(clustData$rowdend2b, rowdend2b)
  expect_equal(clustData$rowdend3, rowdend3)
  expect_equal(clustData$norm.method, "mean")
  expect_equal(clustData$dist.method, "euclidean")
  expect_equal(clustData$agg.method, "complete")
})

context("messages")
test_that("correct messages are thrown",{
  suppressWarnings(meta <- metaData(y = tb.expr[1:100,], design = tb.design, columnname = "columnname",sample.id = "sample_id",
                      subject.id = "monkey_id",long = TRUE,time.var = "timepoint",control.var = "clinical_status", control.val = "Latent"))
  expect_message(normalizeData(meta), "baseline.var and/or baseline.val is unspecified. Cannot produce 
baseline mean normalized, baseline healthy normalized, or all 
samples baseline normalized data.")
  suppressWarnings(meta <- metaData(y = tb.expr[1:100,], design = tb.design, columnname = "columnname",sample.id = "sample_id",
                      subject.id = "monkey_id",long = TRUE,time.var = "timepoint"))
  expect_message(normalizeData(meta), "baseline.var and/or baseline.val is unspecified. Cannot produce 
baseline mean normalized or all samples baseline normalized data.")
  meta <- metaData(y = tb.expr[1:100,], design = tb.design, columnname = "columnname",sample.id = "sample_id",
                      long = TRUE,time.var = "timepoint",baseline.var = "timepoint",baseline.val = 0,
                      control.var = "clinical_status",control.val = "Latent")
  expect_message(normalizeData(meta), "subject.id is not defined. Cannot produce baseline normalized 
data.")
  meta <- metaData(y = tb.expr[1:100,], design = tb.design, columnname = "columnname",sample.id = "sample_id",
                      long = TRUE,time.var = "timepoint",baseline.var = "timepoint",baseline.val = 0)
  expect_message(normalizeData(meta), "subject.id is not defined. Cannot produce baseline normalized 
data.")
})



  
  
  
  
  
  
  
  
  
  
  
  
  