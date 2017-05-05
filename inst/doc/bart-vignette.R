## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ------------------------------------------------------------------------
library(BART)

## ------------------------------------------------------------------------
design$Group <- paste(design$Responder.Status, "_", design$Hours, sep = "")
grp <- factor(design$Group)


## ------------------------------------------------------------------------
design2 <- model.matrix(~0+grp)
colnames(design2) <- levels(grp)

## ------------------------------------------------------------------------
dupcor <- duplicateCorrelation(flu_data, design2, block=design$Subject)

## ------------------------------------------------------------------------
fit <- lmFit(flu_data, design2, block = design$Subject, correlation = dupcor$consensus.correlation)
contrasts <- makeContrasts(H0_SvsA = Symp_0-Asymp_0,H5_SvsA = Symp_5-Asymp_5, H12_SvsA = Symp_12-Asymp_12, H21_SvsA = Symp_21-Asymp_21,
                           H29_SvsA = Symp_29-Asymp_29, H36_SvsA = Symp_36-Asymp_36, H45_SvsA = Symp_45-Asymp_45, H53_SvsA = Symp_53-Asymp_53,
                           H60_SvsA = Symp_60-Asymp_60, H69_SvsA = Symp_69-Asymp_69, H77_SvsA = Symp_77-Asymp_77, H84_SvsA = Symp_84-Asymp_84,
                           H93_SvsA = Symp_93-Asymp_93, H101_SvsA = Symp_101-Asymp_101, H108_SvsA = Symp_108-Asymp_108, S_H5vs0 = Symp_5-Symp_0,
                           S_H12vs0 = Symp_12-Symp_0, S_H21vs0 = Symp_21-Symp_0, S_H29vs0 = Symp_29-Symp_0, S_H36vs0 = Symp_36-Symp_0,
                           S_H45vs0 = Symp_45-Symp_0, S_H53vs0 = Symp_53-Symp_0, S_H60vs0 = Symp_60-Symp_0, S_H69vs0 = Symp_69-Symp_0,
                           S_H77vs0 = Symp_77-Symp_0, S_H84vs0 = Symp_84-Symp_0, S_H93vs0 = Symp_93-Symp_0, S_H101vs0 = Symp_101-Symp_0,
                           S_H108vs0 = Symp_108-Symp_0, A_H5vs0 = Asymp_5-Asymp_0,
                           A_H12vs0 = Asymp_12-Asymp_0, A_H21vs0 = Asymp_21-Asymp_0, A_H29vs0 = Asymp_29-Asymp_0, A_H36vs0 = Asymp_36-Asymp_0,
                           A_H45vs0 = Asymp_45-Asymp_0, A_H53vs0 = Asymp_53-Asymp_0, A_H60vs0 = Asymp_60-Asymp_0, A_H69vs0 = Asymp_69-Asymp_0,
                           A_H77vs0 = Asymp_77-Asymp_0, A_H84vs0 = Asymp_84-Asymp_0, A_H93vs0 = Asymp_93-Asymp_0, A_H101vs0 = Asymp_101-Asymp_0,
                           A_H108vs0 = Asymp_108-Asymp_0,levels=design2)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend = FALSE)


## ------------------------------------------------------------------------
head(design)

## ------------------------------------------------------------------------
des.info <- desInfo(y = flu_data, design = design, data_type = "micro",columnname = "Columnname", long = TRUE, patient_id = "Subject", baseline_var = "Hours", baseline_val = 0, time_var = "Hours", responder_var = "Responder.Status", summary_var = "Age", sample_id = "Sample.ID")

## ------------------------------------------------------------------------
mod_results <- genModelResults(design_info = des.info, data_type = "micro",object = fit2, lm_Fit = fit, method = "limma")
head(mod_results$results[,c(1:3,46,89)])
head(mod_results$resids[,c(1:3,46,89)])

