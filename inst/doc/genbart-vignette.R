## ----setup, echo = FALSE, message = FALSE--------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(genBART)
library(limma)

## ------------------------------------------------------------------------
#library(genBART)
#library(limma)

## ------------------------------------------------------------------------
tb.design$Group <- paste(tb.design$clinical_status,tb.design$timepoint,sep = "")
grp <- factor(tb.design$Group)


## ------------------------------------------------------------------------
design2 <- model.matrix(~0+grp)
colnames(design2) <- levels(grp)

## ------------------------------------------------------------------------
dupcor <- duplicateCorrelation(tb.expr, design2, block = tb.design$monkey_id)

## ------------------------------------------------------------------------
fit <- lmFit(tb.expr, design2, block = tb.design$monkey_id, correlation = dupcor$consensus.correlation)
contrasts <- makeContrasts(A_20vsPre = Active20-Active0,A_30vsPre = Active30-Active0,
                           A_42vsPre = Active42-Active0,A_56vsPre = Active56-Active0,
                           L_20vsPre = Latent20-Latent0,L_30vsPre = Latent30-Latent0,
                           L_42vsPre = Latent42-Latent0,L_56vsPre = Latent56-Latent0,
                           levels=design2)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend = FALSE)

## ------------------------------------------------------------------------
head(tb.design)

## ------------------------------------------------------------------------
des.info <- desInfo(y = tb.expr, design = tb.design, data_type = "micro", columnname = "columnname", long = TRUE, patient_id = "monkey_id", 
                       baseline_var = "timepoint", baseline_val = 0, time_var = "timepoint", responder_var = "clinical_status", 
                       sample_id = "sample_id", project_name = "TB")

## ------------------------------------------------------------------------
mod_results <- genModelResults(design_info = des.info, object = fit2, lm_Fit = fit, method = "limma")

