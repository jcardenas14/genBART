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
cluster.data <- clusterData(norm.data = norm.data, dist.method = "euclidean", agg.method = "complete")

## ------------------------------------------------------------------------
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
fit <- lmFit(tb.expr, design2, block = tb.design$monkey_id, correlation = dupcor$consensus.correlation)
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
dupcor <- duplicateCorrelation(tb.flow, design2, block = tb.flow.des$monkey_id)

fit.flow <- lmFit(tb.flow, design2, block = tb.flow.des$monkey_id, correlation = dupcor$consensus.correlation)

# Add more contrasts since we have all timepoints available. Additional contrasts could be added as well.
contrasts <- makeContrasts(A_20vsPre = Active20-Active0, A_30vsPre = Active30-Active0,
                           A_42vsPre = Active42-Active0, A_56vsPre = Active56-Active0,
                           L_20vsPre = Latent20-Latent0, L_30vsPre = Latent30-Latent0,
                           L_42vsPre = Latent42-Latent0, L_46vsPre = Latent56-Latent0,
                           levels=design2)

fit2.flow <- contrasts.fit(fit.flow, contrasts)
fit2.flow <- eBayes(fit2.flow, trend = FALSE)

## ------------------------------------------------------------------------
mod.results <- genModelResults(y = tb.expr, data.type = "microarray", object = fit2, lm.Fit = fit, method = "limma")
mod.results.flow <- genModelResults(y = tb.flow, data.type = "flow", object = fit2.flow, lm.Fit = fit.flow, method = "limma")

## ------------------------------------------------------------------------
qus <- qBart(model.results = mod.results, gene.sets = modules)

## ------------------------------------------------------------------------
# Create time variable
time <- module.as$time
module.as$time <- NULL

# Format flow data to run correlations and match flow samples with module.as
flow <- data.frame(t(tb.flow))
flow <- flow[match(rownames(module.as), rownames(flow), nomatch = 0), ]

# Run correlations formatted for BART
corrs <- crossCorr(x = module.as, y = flow, by = time, by.name = "days", 
                   description = "Mod.Act.Score vs Flow", x.var = "Mod.Act.Score", 
                   y.var = "Flow", method = "spearman")


## ------------------------------------------------------------------------
genFile(meta = list(meta, meta.flow), module.scores = mod.scores, dendrograms = cluster.data, 
        model.results = list(mod.results, mod.results.flow), project.name = "BART example")

## ---- eval = FALSE-------------------------------------------------------
#  path <- paste(getwd(), "/", "BART example", sep = "")
#  updateFile(load.path = path, qusage.results = qus, corr.results = list(corrs))

## ---- eval = FALSE-------------------------------------------------------
#  runBart()

