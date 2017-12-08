library(testthat)
library(genBart)
library(edgeR)

context("limma")
test_that("genModelResults outputs correctly for one comparison", {
  dat <- tb.expr[1:100,]
  
  # Generate lmFit and eBayes (limma) objects needed for genModelResults
  tb.design$Group <- paste(tb.design$clinical_status,tb.design$timepoint,sep = "")
  grp <- factor(tb.design$Group)
  design2 <- model.matrix(~0+grp)
  colnames(design2) <- levels(grp)
  dupcor <- limma::duplicateCorrelation(dat, design2, block = tb.design$monkey_id)
  fit <- limma::lmFit(dat, design2, block = tb.design$monkey_id, correlation = dupcor$consensus.correlation)
  contrasts <- limma::makeContrasts(A_20vsPre = Active20-Active0, levels=design2)
  fit2 <- limma::contrasts.fit(fit, contrasts)
  fit2 <- limma::eBayes(fit2, trend = FALSE)
  
  # Format results
  model.results <- genModelResults(y = dat, data.type = "microarray", object = fit2,
                                   lm.Fit = fit, method = "limma")
  expect_equal(model.results$data.type, "microarray")
  expect_equal(model.results$gene.sets, NULL)
  expect_equal(model.results$annotations, NULL)
  expect_equal(ncol(model.results$results), 7)
  expect_equal(nrow(model.results$results), 100)
  expect_equal(dim(model.results$resids), dim(dat))
})

test_that("genModelResults outputs correctly for multiple comparisons", {
  dat <- tb.expr[1:100,]
  
  # Generate lmFit and eBayes (limma) objects needed for genModelResults
  tb.design$Group <- paste(tb.design$clinical_status,tb.design$timepoint,sep = "")
  grp <- factor(tb.design$Group)
  design2 <- model.matrix(~0+grp)
  colnames(design2) <- levels(grp)
  dupcor <- limma::duplicateCorrelation(dat, design2, block = tb.design$monkey_id)
  fit <- limma::lmFit(dat, design2, block = tb.design$monkey_id, correlation = dupcor$consensus.correlation)
  contrasts <- limma::makeContrasts(A_20vsPre = Active20-Active0, A_42vsPre = Active42-Active0, levels=design2)
  fit2 <- limma::contrasts.fit(fit, contrasts)
  fit2 <- limma::eBayes(fit2, trend = FALSE)
  
  # Format results
  model.results <- genModelResults(y = dat, data.type = "microarray", object = fit2,
                                   lm.Fit = fit, method = "limma")
  expect_equal(model.results$data.type, "microarray")
  expect_equal(model.results$gene.sets, NULL)
  expect_equal(model.results$annotations, NULL)
  expect_equal(ncol(model.results$results), 12)
  expect_equal(nrow(model.results$results), 100)
  expect_equal(dim(model.results$resids), dim(dat))
})

context("DESeq2")
test_that("genModelResults outputs correctly for one comparison", {
  dat <- round(2^tb.expr[1:100,])
  
  # Generate DESeq2 result objects needed for genModelResults
  group <- paste(tb.design$clinical_status, tb.design$timepoint, sep = "")
  design2 <- data.frame(group = group)
  rownames(design2) <- colnames(dat)
  dds <- DESeq2::DESeqDataSetFromMatrix(dat, colData = design2, design = ~ group)
  dds <- DESeq2::DESeq(dds)
  A_20vsPre = DESeq2::results(dds, contrast = c("group", "Active20", "Active0"))
  
  # Format results
  model.results <- genModelResults(y = dat, data.type = "rnaseq", object = list(A_20vsPre),
                                   lm.Fit = dds, method = "deseq2")
  expect_equal(model.results$data.type, "rnaseq")
  expect_equal(model.results$gene.sets, NULL)
  expect_equal(model.results$annotations, NULL)
  expect_equal(ncol(model.results$results), 6)
  expect_equal(nrow(model.results$results), 100)
})

test_that("genModelResults outputs correctly for multiple comparisons", {
  dat <- round(2^tb.expr[1:100,])
  
  # Generate DESeq2 result objects needed for genModelResults
  group <- paste(tb.design$clinical_status, tb.design$timepoint, sep = "")
  design2 <- data.frame(group = group)
  rownames(design2) <- colnames(dat)
  dds <- DESeq2::DESeqDataSetFromMatrix(dat, colData = design2, design = ~ group)
  dds <- DESeq2::DESeq(dds)
  A_20vsPre = DESeq2::results(dds, contrast = c("group", "Active20", "Active0"))
  A_42vsPre = DESeq2::results(dds, contrast = c("group", "Active42", "Active0"))
  
  # Format results
  model.results <- genModelResults(y = dat, data.type = "rnaseq", object = list(A_20vsPre, A_42vsPre),
                                   lm.Fit = dds, method = "deseq2", comp.names = c("A_20vsPre", "A_42vsPre"))
  expect_equal(model.results$data.type, "rnaseq")
  expect_equal(model.results$gene.sets, NULL)
  expect_equal(model.results$annotations, NULL)
  expect_equal(ncol(model.results$results), 10)
  expect_equal(nrow(model.results$results), 100)
})

context("edgeR")
test_that("genModelResults outputs correctly for one comparison", {
  dat <- round(2^tb.expr[1:100,])
  
  # Generate edgeR result objects needed for genModelResults
  group <- paste(tb.design$clinical_status, tb.design$timepoint, sep = "")
  y <- DGEList(counts=dat,group=group)
  y <- calcNormFactors(y)
  design2 <- model.matrix(~0+group)
  y <- estimateDisp(y,design2)
  fit <- glmQLFit(y,design2)
  my.contrasts <- limma::makeContrasts(A_20vsPre=groupActive20-groupActive0, 
                                       A_42vsPre=groupActive42-groupActive0, levels=design2)
  qlf <- glmQLFTest(fit, contrast = my.contrasts[,"A_20vsPre"])
  fit.lrt <- glmFit(y,design2)
  lrt <- glmLRT(fit.lrt,contrast = my.contrasts[,"A_20vsPre"])
  
  # Format results
  model.results <- genModelResults(y = dat, data.type = "rnaseq", object = list(qlf),
                                   lm.Fit = fit, method = "edgeR", comp.names = "A_20vsPre")
  model.results.lrt <- genModelResults(y = dat, data.type = "rnaseq", object = list(lrt), 
                                       lm.Fit = fit.lrt, method = "edgeR", comp.names = "A_20vsPre")
  expect_equal(model.results$data.type, "rnaseq")
  expect_equal(model.results$gene.sets, NULL)
  expect_equal(model.results$annotations, NULL)
  expect_equal(ncol(model.results$results), 6)
  expect_equal(nrow(model.results$results), 100)
  expect_equal(model.results.lrt$data.type, "rnaseq")
  expect_equal(model.results.lrt$gene.sets, NULL)
  expect_equal(model.results.lrt$annotations, NULL)
  expect_equal(ncol(model.results.lrt$results), 6)
  expect_equal(nrow(model.results.lrt$results), 100)
})

test_that("genModelResults outputs correctly for multiple comparisons", {
  dat <- round(2^tb.expr[1:100,])
  
  # Generate edgeR result objects needed for genModelResults
  group <- paste(tb.design$clinical_status, tb.design$timepoint, sep = "")
  y <- DGEList(counts=dat,group=group)
  y <- calcNormFactors(y)
  design2 <- model.matrix(~0+group)
  y <- estimateDisp(y,design2)
  fit <- glmQLFit(y,design2)
  my.contrasts <- limma::makeContrasts(A_20vsPre=groupActive20-groupActive0, 
                                       A_42vsPre=groupActive42-groupActive0, levels=design2)
  A_20vsPre.qlf <- glmQLFTest(fit, contrast = my.contrasts[,"A_20vsPre"])
  A_42vsPre.qlf <- glmQLFTest(fit, contrast = my.contrasts[,"A_42vsPre"])
  fit.lrt <- glmFit(y,design2)
  A_20vsPre.lrt <- glmLRT(fit.lrt,contrast = my.contrasts[,"A_20vsPre"])
  A_42vsPre.lrt <- glmLRT(fit.lrt,contrast = my.contrasts[,"A_42vsPre"])
  # Format results
  model.results <- genModelResults(y = dat, data.type = "rnaseq", object = list(A_20vsPre.qlf, A_42vsPre.qlf),
                                   lm.Fit = fit, method = "edgeR", comp.names = c("A_20vsPre", "A_42vsPre"))
  model.results.lrt <- genModelResults(y = dat, data.type = "rnaseq", object = list(A_20vsPre.lrt, A_42vsPre.lrt), 
                                       lm.Fit = fit.lrt, method = "edgeR", comp.names = c("A_20vsPre", "A_42vsPre"))
  expect_equal(model.results$data.type, "rnaseq")
  expect_equal(model.results$gene.sets, NULL)
  expect_equal(model.results$annotations, NULL)
  expect_equal(ncol(model.results$results), 10)
  expect_equal(nrow(model.results$results), 100)
  expect_equal(model.results.lrt$data.type, "rnaseq")
  expect_equal(model.results.lrt$gene.sets, NULL)
  expect_equal(model.results.lrt$annotations, NULL)
  expect_equal(ncol(model.results.lrt$results), 10)
  expect_equal(nrow(model.results.lrt$results), 100)
})

