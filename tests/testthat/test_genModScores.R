library(testthat)
library(genBart)

genesets <- list(s1 = rownames(tb.expr)[1:10], s2 = rownames(tb.expr)[11:20])
y <- tb.expr[1:50,]
x <- tb.design

context("cross-sectional data without controls")
test_that("cross-sectional data with no controls outputs NULL", {
  x.base <- x[which(x$timepoint == 0),]
  y.base <- y[,which(x$timepoint == 0)]
  meta <- metaData(y = y.base, design = x.base, columnname = "columnname", long = FALSE)
  mods <- genModScores(meta, genesets)
  expect_equal(mods$scores.ctrl, NULL)
  expect_equal(mods$scores.base, NULL)
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
  for(i in 1:length(genesets)){
    counts[[i]] <- apply(test.mat[which(rownames(test.mat) %in% genesets[[i]]),], 2, sum)
  }
  counts <- do.call(rbind, counts)
  rownames(counts) <- names(genesets)
  colnames(counts) <- x.case$sample_id
  sets.length <- c()
  for(i in 1:length(genesets)){
    sets.length[i] <- length(genesets[[i]])
  }
  percents <- counts/sets.length + 1
  meta <- metaData(y = y.base, design = x.base, columnname = "columnname", long = FALSE,
                      control.var = "clinical_status", control.val = "Latent",
                      baseline.var = "timepoint", baseline.val = 0, sample.id = "sample_id",
                      subject.id = "monkey_id")
  mods <- genModScores(meta, genesets)
  expect_equal(mods$scores.ctrl, percents)
  expect_equal(mods$scores.base, NULL)
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
  for(i in 1:length(genesets)){
    counts[[i]] <- apply(test.mat[which(rownames(test.mat) %in% genesets[[i]]),], 2, sum)
  }
  counts <- do.call(rbind, counts)
  rownames(counts) <- names(genesets)
  colnames(counts) <- x.nobase$sample_id
  sets.length <- c()
  for(i in 1:length(genesets)){
    sets.length[i] <- length(genesets[[i]])
  }
  percents <- counts/sets.length + 1
  meta <- metaData(y = y, design = x, columnname = "columnname", long = TRUE, time.var = "timepoint",
                      baseline.var = "timepoint", baseline.val = 0, sample.id = "sample_id",
                      subject.id = "monkey_id")
  mods <- genModScores(meta, genesets)
  expect_equal(mods$scores.base, percents)
  expect_equal(mods$scores.ctrl, NULL)
})

context("longitudinal data with controls")
test_that("longitudinal data with controls outputs correctly", {
  x$clinical_status <- as.character(x$clinical_status)
  x$clinical_status[which(x$clinical_status == "Latent" & x$timepoint == 0)] <- "control"
  x.base <- x[which(x$timepoint == 0 & x$clinical_status == "Active"),]
  y.base <- y[,which(x$timepoint == 0 & x$clinical_status == "Active")]
  x.ctrl <- x[which(x$clinical_status == "control"),]
  y.ctrl <- y[,which(x$clinical_status == "control")]
  x.noctrl <- x[which(x$clinical_status != "control"),]
  y.noctrl <- y[,which(x$clinical_status != "control")]
  x.long <- x[which(x$clinical_status != "control" & x$timepoint != 0),]
  y.long <- y[,which(x$clinical_status != "control" & x$timepoint != 0)]
  base.mean <- apply(y.base, 1, mean)
  base.sd <- apply(y.base, 1, sd)
  c.mean <- apply(y.ctrl, 1, mean)
  c.sd <- apply(y.ctrl, 1, sd)
  test.mat <- as.matrix(y.long)
  test.mat[is.numeric(test.mat)] <- 0
  test.mat[y.long < base.mean-2*base.sd] <- -1
  test.mat[y.long > base.mean+2*base.sd] <- 1
  counts <- list()
  for(i in 1:length(genesets)){
    counts[[i]] <- apply(test.mat[which(rownames(test.mat) %in% genesets[[i]]),], 2, sum)
  }
  counts <- do.call(rbind, counts)
  rownames(counts) <- names(genesets)
  colnames(counts) <- x.long$sample_id
  sets.length <- c()
  for(i in 1:length(genesets)){
    sets.length[i] <- length(genesets[[i]])
  }
  scores.base <- counts/sets.length + 1
  scores.base <- scores.base[,order(x.long$timepoint, x.long$monkey_id)]
  test.mat <- as.matrix(y.noctrl)
  test.mat[is.numeric(test.mat)] <- 0
  test.mat[y.noctrl < c.mean-2*c.sd] <- -1
  test.mat[y.noctrl > c.mean+2*c.sd] <- 1
  counts <- list()
  for(i in 1:length(genesets)){
    counts[[i]] <- apply(test.mat[which(rownames(test.mat) %in% genesets[[i]]),], 2, sum)
  }
  counts <- do.call(rbind, counts)
  rownames(counts) <- names(genesets)
  colnames(counts) <- x.noctrl$sample_id
  sets.length <- c()
  for(i in 1:length(genesets)){
    sets.length[i] <- length(genesets[[i]])
  }
  scores.ctrl <- counts/sets.length + 1
  
  meta <- metaData(y = y, design = x, columnname = "columnname", long = TRUE, time.var = "timepoint",
                      control.var = "clinical_status", control.val = "control",
                      baseline.var = "timepoint", baseline.val = 0, sample.id = "sample_id",
                      subject.id = "monkey_id")
  mods <- genModScores(meta, genesets)
  expect_equal(mods$scores.ctrl, scores.ctrl)
  expect_equal(mods$scores.base, scores.base)
})





