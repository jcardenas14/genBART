## ----setup, echo = FALSE, message = FALSE--------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(genBART)
library(limma)

## ------------------------------------------------------------------------
head(tb.design)

## ------------------------------------------------------------------------
des.info <- desInfo(y = tb.expr, design = tb.design, data_type = "micro", 
                    columnname = "columnname", long = TRUE, 
                    sample_id = "sample_id", patient_id = "monkey_id", 
                    time_var = "timepoint", baseline_var = "timepoint", baseline_val = 0, 
                    control_var = "clinical_status", control_val = "Latent", 
                    summary_var = NULL, responder_var = "clinical_status",  
                    project_name = "TB")

## ------------------------------------------------------------------------
des.info.flow <- desInfo(y = tb.flow, design = tb.flow.des, data_type = "flow", 
                         columnname = "columnname", long = TRUE, patient_id = "monkey_id", 
                         baseline_var = "timepoint", baseline_val = 0, 
                         control_var = "clinical_status", control_val = "Latent", 
                         time_var = "timepoint", summary_var = NULL, 
                         responder_var = "clinical_status", sample_id = "sample_id", 
                         project_name = "TB")

## ------------------------------------------------------------------------
dendros <- genDendrograms(design_info = des.info)

## ------------------------------------------------------------------------
mods <- genModules(design_info = des.info, gene_sets = clusters)

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
mod_results <- genModelResults(design_info = des.info, object = fit2, lm_Fit = fit, method = "limma")
mod_results.flow <- genModelResults(design_info = des.info.flow, object = fit2.flow, lm_Fit = fit.flow,
                                    method = "limma")

## ------------------------------------------------------------------------
qus <- qBart(object = mod_results, gene_sets = modules)

## ------------------------------------------------------------------------
# Create time variable
time <- module.as$time
module.as$time <- NULL

# Format flow data to run correlations and match flow samples with module.as
flow <- data.frame(t(tb.flow))
flow <- flow[match(rownames(module.as), rownames(flow), nomatch = 0), ]

# Run correlations formatted for BART
corrs <- crossCorr(x = module.as, y = flow, by = time, by_name = "days", 
                   description = "Mod.Act.Score vs Flow", x_var = "Mod.Act.Score", 
                   y_var = "Flow", method = "spearman")


## ------------------------------------------------------------------------
genFile(design_info = list(des.info, des.info.flow), module_maps = mods, dendros = dendros, 
        model_results = list(mod_results, mod_results.flow))

## ---- eval = FALSE-------------------------------------------------------
#  path <- paste(getwd(), "/", des.info$project_name, " Pipeline", sep = "")
#  updateFile(load.path = path, qusage_results = qus, corr_results = list(corrs))

## ---- eval = FALSE-------------------------------------------------------
#  runBart()

