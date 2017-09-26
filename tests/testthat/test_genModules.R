library(testthat)
library(genBART)

genesets <- list(s1 = rownames(tb.expr)[1:10], s2 = rownames(tb.expr)[11:20])
y <- tb.expr[1:50,]
x <- tb.design

context("cross-sectional data without controls")
test_that("cross-sectional data with no controls outputs NULL", {
  x.base <- x[which(x$timepoint == 0),]
  y.base <- y[,which(x$timepoint == 0)]
  des.info <- desInfo(y = y.base, design = x.base, columnname = "columnname", long = FALSE)
  mods <- genModules(des.info, genesets)
  expect_equal(mods$base_ctrl, NULL)
  expect_equal(mods$long_base, NULL)
  expect_equal(mods$long_ctrl, NULL)
})

context("cross-sectional data with controls")
test_that("cross-sectional data with controls outputs correctly", {
  x.base <- x[which(x$timepoint == 0),]
  y.base <- y[,which(x$timepoint == 0)]
  x.case <- x.base[which(x.base$clinical_status == "Active"),]
  y.case <- y.base[,which(x.base$clinical_status == "Active")]
  c.mean <- apply(y.base[,which(x.base$clinical_status == "Latent")], 1, mean)
  c.sd <- apply(y.base[,which(x.base$clinical_status == "Latent")], 1, sd)
  test.mat <- as.matrix(y.case)
  test.mat[is.numeric(test.mat)] <- 0
  test.mat[y.case < c.mean-2*c.sd] <- -1
  test.mat[y.case > c.mean+2*c.sd] <- 1
  counts <- list()
  for(i in 1:length(geneset)){
    counts[[i]] <- apply(test.mat[which(rownames(test.mat) %in% geneset[[i]]),], 2, sum)
  }
  counts <- do.call(rbind, counts)
  rownames(counts) <- names(genesets)
  colnames(counts) <- x.case$sample_id
  sets.length <- c()
  for(i in 1:length(genesets)){
    sets.length[i] <- length(genesets[[i]])
  }
  percents <- counts/sets.length + 1
  des.info <- desInfo(y = y.base, design = x.base, columnname = "columnname", long = FALSE,
                      control_var = "clinical_status", control_val = "Latent",
                      baseline_var = "timepoint", baseline_val = 0, sample_id = "sample_id",
                      patient_id = "monkey_id")
  mods <- genModules(des.info, genesets)
  expect_equal(mods$base_ctrl, percents)
  expect_equal(mods$long_base, NULL)
  expect_equal(mods$long_ctrl, NULL)
})

context("cross-sectional data with controls")
test_that("cross-sectional data with controls outputs correctly", {
  x.base <- x[which(x$timepoint == 0),]
  y.base <- y[,which(x$timepoint == 0)]
  x.case <- x.base[which(x.base$clinical_status == "Active"),]
  y.case <- y.base[,which(x.base$clinical_status == "Active")]
  c.mean <- apply(y.base[,which(x.base$clinical_status == "Latent")], 1, mean)
  c.sd <- apply(y.base[,which(x.base$clinical_status == "Latent")], 1, sd)
  test.mat <- as.matrix(y.case)
  test.mat[is.numeric(test.mat)] <- 0
  test.mat[y.case < c.mean-2*c.sd] <- -1
  test.mat[y.case > c.mean+2*c.sd] <- 1
  counts <- list()
  for(i in 1:length(geneset)){
    counts[[i]] <- apply(test.mat[which(rownames(test.mat) %in% geneset[[i]]),], 2, sum)
  }
  counts <- do.call(rbind, counts)
  rownames(counts) <- names(genesets)
  colnames(counts) <- x.case$sample_id
  sets.length <- c()
  for(i in 1:length(genesets)){
    sets.length[i] <- length(genesets[[i]])
  }
  percents <- counts/sets.length + 1
  des.info <- desInfo(y = y.base, design = x.base, columnname = "columnname", long = FALSE,
                      control_var = "clinical_status", control_val = "Latent",
                      baseline_var = "timepoint", baseline_val = 0, sample_id = "sample_id",
                      patient_id = "monkey_id")
  mods <- genModules(des.info, genesets)
  expect_equal(mods$base_ctrl, percents)
  expect_equal(mods$long_base, NULL)
  expect_equal(mods$long_ctrl, NULL)
})

context("longitudinal data without controls")
test_that("longitudinal data without controls outputs correctly", {
  x.nobase <- x[-which(x$timepoint == 0),]
  y.nobase <- y[,-which(x$timepoint == 0)]
  y.nobase <- y.nobase[,order(x.nobase$timepoint)]
  x.nobase <- x.nobase[order(x.nobase$timepoint),]
  base.mean <- apply(y[,which(x$timepoint == 0)], 1, mean)
  base.sd <- apply(y[,which(x$timepoint == 0)], 1, sd)
  test.mat <- as.matrix(y.nobase)
  test.mat[is.numeric(test.mat)] <- 0
  test.mat[y.nobase < base.mean-2*base.sd] <- -1
  test.mat[y.nobase > base.mean+2*base.sd] <- 1
  counts <- list()
  for(i in 1:length(geneset)){
    counts[[i]] <- apply(test.mat[which(rownames(test.mat) %in% geneset[[i]]),], 2, sum)
  }
  counts <- do.call(rbind, counts)
  rownames(counts) <- names(genesets)
  colnames(counts) <- x.nobase$sample_id
  sets.length <- c()
  for(i in 1:length(genesets)){
    sets.length[i] <- length(genesets[[i]])
  }
  percents <- counts/sets.length + 1
  des.info <- desInfo(y = y, design = x, columnname = "columnname", long = TRUE, time_var = "timepoint",
                      baseline_var = "timepoint", baseline_val = 0, sample_id = "sample_id",
                      patient_id = "monkey_id")
  mods <- genModules(des.info, genesets)
  expect_equal(mods$base_ctrl, NULL)
  expect_equal(mods$long_base, percents)
  expect_equal(mods$long_ctrl, NULL)
})

context("longitudinal data with controls")
test_that("longitudinal data with controls outputs correctly", {
  x.base <- x[which(x$timepoint == 0 & x$clinical_status == "Active"),]
  y.base <- y[,which(x$timepoint == 0 & x$clinical_status == "Active")]
  x.ctrl <- x[which(x$timepoint == 0 & x$clinical_status == "Latent"),]
  y.ctrl <- y[,which(x$timepoint == 0 & x$clinical_status == "Latent")]
  x.noctrl <- x[which(x$timepoint != 0 & x$clinical_status == "Active"),]
  y.noctrl <- y[,which(x$timepoint != 0 & x$clinical_status == "Active")]
  y.noctrl <- y.noctrl[,order(x.noctrl$timepoint)]
  x.noctrl <- x.noctrl[order(x.noctrl$timepoint),]
  x.long <- x[which(x$clinical_status == "Active"),]
  y.long <- y[,which(x$clinical_status == "Active")]
  y.long <- y.long[,order(x.long$timepoint)]
  x.long <- x.long[order(x.long$timepoint),]
  base.mean <- apply(y.base, 1, mean)
  base.sd <- apply(y.base, 1, sd)
  c.mean <- apply(y.ctrl, 1, mean)
  c.sd <- apply(y.ctrl, 1, sd)
  test.mat <- as.matrix(y.noctrl)
  test.mat[is.numeric(test.mat)] <- 0
  test.mat[y.noctrl < base.mean-2*base.sd] <- -1
  test.mat[y.noctrl > base.mean+2*base.sd] <- 1
  counts <- list()
  for(i in 1:length(geneset)){
    counts[[i]] <- apply(test.mat[which(rownames(test.mat) %in% geneset[[i]]),], 2, sum)
  }
  counts <- do.call(rbind, counts)
  rownames(counts) <- names(genesets)
  colnames(counts) <- x.noctrl$sample_id
  sets.length <- c()
  for(i in 1:length(genesets)){
    sets.length[i] <- length(genesets[[i]])
  }
  long_base <- counts/sets.length + 1
  
  test.mat <- as.matrix(y.long)
  test.mat[is.numeric(test.mat)] <- 0
  test.mat[y.long < c.mean-2*c.sd] <- -1
  test.mat[y.long > c.mean+2*c.sd] <- 1
  counts <- list()
  for(i in 1:length(geneset)){
    counts[[i]] <- apply(test.mat[which(rownames(test.mat) %in% geneset[[i]]),], 2, sum)
  }
  counts <- do.call(rbind, counts)
  rownames(counts) <- names(genesets)
  colnames(counts) <- x.long$sample_id
  sets.length <- c()
  for(i in 1:length(genesets)){
    sets.length[i] <- length(genesets[[i]])
  }
  long_ctrl <- counts/sets.length + 1
  
  test.mat <- as.matrix(y.base)
  test.mat[is.numeric(test.mat)] <- 0
  test.mat[y.base < c.mean-2*c.sd] <- -1
  test.mat[y.base > c.mean+2*c.sd] <- 1
  counts <- list()
  for(i in 1:length(geneset)){
    counts[[i]] <- apply(test.mat[which(rownames(test.mat) %in% geneset[[i]]),], 2, sum)
  }
  counts <- do.call(rbind, counts)
  rownames(counts) <- names(genesets)
  colnames(counts) <- x.base$sample_id
  sets.length <- c()
  for(i in 1:length(genesets)){
    sets.length[i] <- length(genesets[[i]])
  }
  base_ctrl <- counts/sets.length + 1
  
  des.info <- desInfo(y = y, design = x, columnname = "columnname", long = TRUE, time_var = "timepoint",
                      control_var = "clinical_status", control_val = "Latent",
                      baseline_var = "timepoint", baseline_val = 0, sample_id = "sample_id",
                      patient_id = "monkey_id")
  mods <- genModules(des.info, genesets)
  expect_equal(mods$base_ctrl, base_ctrl)
  expect_equal(mods$long_base, long_base)
  expect_equal(mods$long_ctrl, long_ctrl)
})





