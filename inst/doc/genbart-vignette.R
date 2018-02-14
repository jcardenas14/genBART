## ----setup, echo = FALSE, message = FALSE--------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(genBart)
library(limma)

## ------------------------------------------------------------------------
head(tb.design)

## ------------------------------------------------------------------------
meta <- metaData(y = tb.expr, design = tb.design, data.type = "microarray", 
                 columnname = "columnname", long = TRUE, sample.id = "sample_id", 
                 subject.id = "monkey_id", time.var = "timepoint", 
                 baseline.var = "timepoint", baseline.val = 0)

## ------------------------------------------------------------------------
meta.flow <- metaData(y = tb.flow, design = tb.flow.des, data.type = "flow", 
                      columnname = "columnname", long = TRUE, sample.id = "sample_id", 
                      subject.id = "monkey_id", time.var = "timepoint", 
                      baseline.var = "timepoint", baseline.val = 0)

## ------------------------------------------------------------------------
norm.data <- normalizeData(meta = meta, norm.method = "mean")
cluster.data <- clusterData(norm.data = norm.data, dist.method = "euclidean", 
                            agg.method = "complete")

## ------------------------------------------------------------------------
# Form clusters
tb.expr.scale <- data.frame(t(scale(t(tb.expr)))) # center and scale probes
hc <- fastcluster::hclust(dist(tb.expr.scale))
cls <- cutree(hc, 10)
clusters <- list()
for(i in 1:10){
  clusters[[i]] <- names(cls)[which(cls == i)]
  names(clusters)[i] <- paste0("C",i)
}

# Generate module scores
mod.scores <- genModScores(meta = meta, gene.sets = clusters, sd.lim = 2)

## ------------------------------------------------------------------------
tb.design$Group <- paste(tb.design$clinical_status, tb.design$timepoint, sep = "")
grp <- factor(tb.design$Group)


## ------------------------------------------------------------------------
design2 <- model.matrix(~0+grp)
colnames(design2) <- levels(grp)

## ------------------------------------------------------------------------
dupcor <- duplicateCorrelation(tb.expr, design2, block = tb.design$monkey_id)

## ------------------------------------------------------------------------
fit <- lmFit(tb.expr, design2, block = tb.design$monkey_id, 
             correlation = dupcor$consensus.correlation)
contrasts <- makeContrasts(A_20vsPre = Active20-Active0, A_42vsPre = Active42-Active0,
                           L_20vsPre = Latent20-Latent0, L_42vsPre = Latent42-Latent0,
                           levels=design2)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend = FALSE)

## ------------------------------------------------------------------------
tb.flow.des$Group <- paste(tb.flow.des$clinical_status, tb.flow.des$timepoint, sep = "")
grp <- factor(tb.flow.des$Group)

design2 <- model.matrix(~0+grp)
colnames(design2) <- levels(grp)
tb.flow.l2 <- log2(tb.flow + 1) # log2 transform
dupcor <- duplicateCorrelation(tb.flow.l2, design2, block = tb.flow.des$monkey_id)
fit.flow <- lmFit(tb.flow.l2, design2, block = tb.flow.des$monkey_id,
correlation = dupcor$consensus.correlation)

# Additional contrasts since all timepoints available. 
contrasts <- makeContrasts(A_20vsPre = Active20-Active0, A_30vsPre = Active30-Active0,
                           A_42vsPre = Active42-Active0, A_56vsPre = Active56-Active0,
                           L_20vsPre = Latent20-Latent0, L_30vsPre = Latent30-Latent0,
                           L_42vsPre = Latent42-Latent0, L_46vsPre = Latent56-Latent0,
                           levels=design2)

fit2.flow <- contrasts.fit(fit.flow, contrasts)
fit2.flow <- eBayes(fit2.flow, trend = FALSE)

## ------------------------------------------------------------------------
data(gene.symbols) # gene symbols for Illumina probes
mod.results <- genModelResults(y = tb.expr, data.type = "microarray", object = fit2, 
                               lm.Fit = fit, method = "limma", var.symbols = gene.symbols)
mod.results.flow <- genModelResults(y = tb.flow, data.type = "flow", object = fit2.flow, 
                                    lm.Fit = fit.flow, method = "limma")

## ------------------------------------------------------------------------
qus <- runQgen(model.results = mod.results, gene.sets = modules)

## ------------------------------------------------------------------------
# Format expression data to align with flow data
gene.data <- data.frame(t(tb.expr))
rownames(gene.data) <- paste0(tb.design$monkey_id, "_", tb.design$timepoint)
flow.data <- data.frame(t(tb.flow))
flow.data <- flow.data[match(rownames(gene.data), rownames(flow.data), nomatch = 0), ]
gene.data <- gene.data[match(rownames(flow.data), rownames(gene.data), nomatch = 0), ]

# Create time variable
time <- tb.flow.des$timepoint[match(rownames(flow.data),tb.flow.des$columnname,nomatch = 0)]

# Run correlations formatted for BART
corrs <- crossCorr(x = gene.data, y = flow.data, by = time, by.name = "days", 
                   description = "Genes vs Flow", x.var = "Genes", 
                   y.var = "Flow", method = "spearman")


## ---- eval = FALSE-------------------------------------------------------
#  genFile(meta = list(meta, meta.flow), module.scores = mod.scores,
#          dendrograms = cluster.data, model.results = list(mod.results, mod.results.flow),
#          project.name = "BART example")

## ---- eval = FALSE-------------------------------------------------------
#  path <- paste0(getwd(), "/", "BART example")
#  updateFile(load.path = path, qusage.results = qus, corr.results = list(corrs))

## ---- eval = FALSE-------------------------------------------------------
#  runBart()

