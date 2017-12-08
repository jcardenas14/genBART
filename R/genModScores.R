#' Generate module (gene set) maps for plotting
#' 
#' @param meta list generated from \code{metaData}
#' @param gene.sets list of gene sets
#' @param sd.lim number of standard deviations away from the mean of the 
#'   reference samples
#' @details This function calculates module scores for individual samples. In 
#'   cross sectional studies with controls, the control samples are used to 
#'   determine an upper and lower threshold (mean of controls +/- 2 sd). The 
#'   module proportion for each sample is then calculated based on the 
#'   percentage of genes within a module that are above or below this threshold. 
#'   For example, if 40\% of the genes within a module are above the threshold 
#'   and 15\% are below it, the final module score would be 25\% up (40-15). In 
#'   longitudinal settings, module scores are calculated with respect to 
#'   controls and baseline samples. In cross sectional studies without controls, 
#'   \code{genModules} cannot be used, since there are no reference samples with 
#'   which to calculate a threshold.
#' @return \code{scores.ctrl} data frame of module scores for all samples 
#'   with respect to controls.
#' @return \code{scores.base} data frame of module scores for all time point 
#'   samples with respect to their baseline.
#' @return \code{gene.sets} List of gene sets provided through \code{gene.sets}
#' @examples
#' # Example data
#' data(tb.expr)
#' data(tb.design)
#' data(modules)
#' 
#' # Demonstrate on first 100 probes
#' dat <- tb.expr[1:100, ]
#' 
#' # Create desInfo object
#' meta.data <- metaData(y = tb.expr, design = tb.design, data.type = "microarray", 
#'                     columnname = "columnname", long = TRUE, sample.id = "sample_id",
#'                     subject.id = "monkey_id", time.var = "timepoint",
#'                     baseline.var = "timepoint", baseline.val = 0)
#' 
#' # Generate module maps                                      
#' mods <- genModScores(meta.data, modules)
#' @export
genModScores <- function(meta, gene.sets, sd.lim = 2) {
  y <- meta$y
  design <- meta$design
  Transcript.ID <- rownames(y)
  exprs <- cbind(Transcript.ID = Transcript.ID, y)
  baseline.var <- meta$baseline.var
  baseline.val <- meta$baseline.val
  control.val <- meta$control.val
  control.var <- meta$control.var
  sample.id <- meta$sample.id
  subject.id <- meta$subject.id
  time.var <- meta$time.var
  long <- meta$long
  genes <- unlist(gene.sets)
  gset.names <- rep(names(gene.sets), times = lapply(gene.sets, length))
  gsets <- data.frame(Transcript.ID = genes, Module = gset.names)
  rownames(gsets) <- NULL
  modExp <- merge(exprs, gsets, by = c("Transcript.ID"))
  modOrdNum <- table(factor(gsets$Module))
  modNames <- names(modOrdNum)
  if (is.null(control.var)) {
    if (long) {
      if (is.null(baseline.var) || is.null(baseline.val)) {
        warning("Baseline variable unspecified. Must define baseline variable 
when long = TRUE. Returning NULL.")
        z <- list(scores.ctrl = NULL, scores.base = NULL, gene.sets = NULL)
        return(z)
      }
      desOrd <- design
      if (!is.null(time.var) & !is.null(subject.id)) {
        desOrd <- design[order(design[, time.var], design[, subject.id]), ]
      } 
      expOrd <- modExp[, match(desOrd[, "columnname"], colnames(modExp), 
                               nomatch = 0)]
      if (!is.null(sample.id)) {
        colnames(expOrd) <- desOrd[, sample.id]
      }
      baseData <- expOrd[, which(desOrd[, baseline.var] == baseline.val)]
      baseMean <- apply(baseData, 1, mean, na.rm = TRUE)
      baseSd <- apply(baseData, 1, sd, na.rm = TRUE)
      finalData <- as.matrix(expOrd[, -which(colnames(expOrd) %in% 
                                               colnames(baseData))])
      rownames(finalData) <- modExp$Module
      signMat <- finalData
      signMat[is.numeric(finalData)] <- 0
      signMat[finalData < (baseMean - sd.lim * baseSd)] <- -1
      signMat[finalData > (baseMean + sd.lim * baseSd)] <- 1
      countMat <- c()
      n <- length(modOrdNum)
      for (i in 1:n) {
        subset <- signMat[which(rownames(signMat) == modNames[i]), ]
        if (is.vector(subset)) {
          countMat <- rbind(countMat, subset)
        }
        if (is.matrix(subset)) {
          if (dim(subset)[1] == 0) {
            countMat <- rbind(countMat, rep(0, dim(signMat)[2]))
          }
          if (dim(subset)[1] > 1) {
            countMat <- rbind(countMat, apply(subset, 2, sum, na.rm = TRUE))
          }
        }
      }
      rownames(countMat) <- modNames
      percentMat <- countMat / as.vector(modOrdNum)
      scores.base <- percentMat + 1
      scores.ctrl <- NULL
    }  else {
      warning("long = FALSE and no control variable specified.")
      z <- list(scores.ctrl = NULL, scores.base = NULL, gene.sets = NULL)
      return(z)
    }
  } else {
    desCtrl <- design[design[, control.var] == control.val, ]
    expCtrl <- modExp[, match(desCtrl[, "columnname"], colnames(modExp), 
                              nomatch = 0)]
    desCase <- design[design[, control.var] != control.val, ]
    expCase <- as.matrix(modExp[, match(desCase[, "columnname"], 
                                        colnames(modExp), nomatch = 0)])
    rownames(expCase) <- modExp$Module
    if (!is.null(sample.id)) {
      colnames(expCase) <- desCase[, sample.id]
    }
    ctrlMean <- apply(expCtrl, 1, mean, na.rm = TRUE)
    ctrlSd <- apply(expCtrl, 1, sd, na.rm = TRUE)
    signMat <- expCase
    signMat[is.numeric(expCase)] <- 0
    signMat[expCase < (ctrlMean - sd.lim * ctrlSd)] <- -1
    signMat[expCase > (ctrlMean + sd.lim * ctrlSd)] <- 1
    countMat <- c()
    for (i in 1:length(modOrdNum)) {
      subset <- signMat[which(rownames(signMat) == modNames[i]), ]
      if (is.vector(subset)) {
        countMat <- rbind(countMat, subset)
      }
      if (is.matrix(subset)) {
        if (dim(subset)[1] == 0) {
          countMat <- rbind(countMat, rep(0, dim(signMat)[2]))
        }
        if (dim(subset)[1] > 1) {
          countMat <- rbind(countMat, apply(subset, 2, sum, na.rm = TRUE))
        }
      }
    }
    rownames(countMat) <- modNames
    percentMat <- countMat / as.vector(modOrdNum)
    scores.ctrl <- percentMat + 1
    if (long) {
      if (is.null(baseline.var) || is.null(baseline.val)) {
        warning("Baseline variable unspecified. Must specify baseline variable 
when long = TRUE to calculate scores.base.")
        z <- list(scores.ctrl = scores.ctrl, scores.base = NULL, gene.sets = 
                    gene.sets)
        return(z)
      }
      desCase <- design[design[, control.var] != control.val, ]
      desOrd <- desCase
      if (!is.null(time.var) & !is.null(subject.id)) {
        desOrd <- desCase[order(desCase[, time.var], desCase[, subject.id]), ]
      } 
      expOrd <- modExp[, match(desOrd[, "columnname"], colnames(modExp),
                               nomatch = 0)]
      if (!is.null(sample.id)) {
        colnames(expOrd) <- desOrd[, sample.id]
      }
      baseData <- expOrd[, desOrd[, baseline.var] == baseline.val]
      baseMean <- apply(baseData, 1, mean, na.rm = TRUE)
      baseSd <- apply(baseData, 1, sd, na.rm = TRUE)
      finalData <- as.matrix(expOrd[, !colnames(expOrd) %in% 
                                           colnames(baseData)])
      rownames(finalData) <- modExp$Module
      signMat <- finalData
      signMat[is.numeric(finalData)] <- 0
      signMat[finalData < (baseMean - sd.lim * baseSd)] <- -1
      signMat[finalData > (baseMean + sd.lim * baseSd)] <- 1
      countMat <- c()
      for (i in 1:length(modOrdNum)) {
        subset <- signMat[which(rownames(signMat) == modNames[i]), ]
        if (is.vector(subset)) {
          countMat <- rbind(countMat, subset)
        }
        if (is.matrix(subset)) {
          if (dim(subset)[1] == 0) {
            countMat <- rbind(countMat, rep(0, dim(signMat)[2]))
          }
          if (dim(subset)[1] > 1) {
            countMat <- rbind(countMat, apply(subset, 2, sum, na.rm = TRUE))
          }
        }
      }
      rownames(countMat) <- modNames
      percentMat <- countMat / as.vector(modOrdNum)
      scores.base <- percentMat + 1
    } else {
      scores.base <- NULL
    }
  }
  z <- list(scores.ctrl = scores.ctrl, scores.base = scores.base, gene.sets = 
              gene.sets)
  return(z)
}
