#' Generate and Update BART Result files
#' 
#' Generate/Update BART file that can be directly uploaded into the BART app
#' @param meta A list of a single or multiple objects returned by 
#'   \code{metaData}. Default is NULL. 
#' @param model.results A list of a single or multiple objects returned by 
#'   \code{genModelResults}. Default is NULL. 
#' @param module.scores Object returned by \code{genModScores}. Default is 
#'   NULL. 
#' @param dendrograms Object returned by \code{clusterData}. Default is NULL. 
#' @param qusage.results Object returned by \code{runQgen}. Default is NULL.
#' @param roast.results Object returned by \code{rBart}. Default is NULL. 
#' @param corr.results A list of a single or multiple objects returned by 
#'   \code{crossCorr}. Default is NULL. 
#' @param project.name String giving the name project name. Default is 
#'   "BART Project".
#' @param folder.path Path to create folder that BART file will be saved to. If 
#'   NULL (default), BART folder will be created in the current working 
#'   directory.
#' @details \code{genFile} generates a formatted R data file (bartResults.rda) 
#'   that can be uploaded into the BART web application. The file is created 
#'   based on result objects returned by the various functions in 
#'   \code{genBart}. Since BART can store and display results across multiple 
#'   platforms at one time (e.g. RNA-seq, flow cytometry, metabolomics), some of 
#'   the parameters (i.e. meta, model.results, corr.results) require the input 
#'   objects to be wrapped in a list. The function saves the file in a folder 
#'   whose name is given by \code{project.name}. The folder is created in the 
#'   current working directory.
#' @examples
#' # Example data
#' data(tb.expr)
#' data(tb.design)
#' 
#' # Use first 100 probes to demonstrate
#' dat <- tb.expr[1:100,]
#' 
#' # Create desInfo object
#' meta.data <- metaData(y = dat, design = tb.design, data.type = "microarray", 
#'                     columnname = "columnname", long = TRUE, subject.id = "monkey_id",
#'                     baseline.var = "timepoint", baseline.val = 0, time.var = "timepoint", 
#'                     sample.id = "sample_id")
#'                     
#' # create BART file (minimal example)
#' genFile(meta = list(meta.data), folder.path = tempdir())                     
#' 
#' # generate module scores, normalize and cluster genes
#' mods <- genModScores(meta.data, modules)
#' data.norm <- normalizeData(meta = meta.data)
#' dendros <- clusterData(norm.data = data.norm)
#' 
#' # Update BART file with module scores and clustered genes
#' path <- paste0(tempdir(), "/", "BART Project")
#' updateFile(load.path = path, module.scores = mods, dendrograms = dendros)
#' @export
genFile <- function(meta, model.results = NULL, module.scores = NULL,
                    dendrograms = NULL, qusage.results = NULL,
                    roast.results = NULL, corr.results = NULL, 
                    project.name = "BART Project", folder.path = NULL) {
  exprs <- design <- long <- sample.id <- subject.id <- control.var <- 
    control.val <- baseline.var <- baseline.val <- time.var <- scores.base <- 
    scores.ctrl <- modules <- rowdend1b <- rowdend2b <- rowdend1 <- rowdend2 <-
    rowdend3 <- norm.method <- dist.method <- agg.method <- results.file <- 
    dge.gsets <- dge.annots <- lower.ci <- upper.ci <- gene.sets <- annots <- 
    corrs <- corr.files <- corr.num <- corr.names <- x.var <- y.var <- 
    corr.method <- flow.results <- flow.data <- metab.results <- metab.data <- 
    data.type <- NULL
  if (!is.null(module.scores)) {
    scores.base <- module.scores$scores.base
    scores.ctrl <- module.scores$scores.ctrl
  }
  if (!is.null(dendrograms)) {
    rowdend1b <- dendrograms$rowdend1b
    rowdend2b <- dendrograms$rowdend2b
    rowdend1 <- dendrograms$rowdend1
    rowdend2 <- dendrograms$rowdend2
    rowdend3 <- dendrograms$rowdend3
    norm.method <- dendrograms$norm.method
    dist.method <- dendrograms$dist.method
    agg.method <- dendrograms$agg.method
  }
  if (!is.null(qusage.results)) {
    lower.ci <- qusage.results$lower.ci
    upper.ci <- qusage.results$upper.ci
    gene.sets <- qusage.results$gene.sets
    annots <- qusage.results$annotations
    qusage.results <- qusage.results$qusage.results
  }
  if (!is.null(corr.results)) {
    if (is.list(corr.results)) {
      corr.names <- y.var <- x.var <- corr.method <- c()
      corrs <- corr.files <- list()
      for (i in 1:length(corr.results)) {
        corr.names[i] <- corr.results[[i]]$corr.names
        y.var[i] <- corr.results[[i]]$y.var
        x.var[i] <- corr.results[[i]]$x.var
        corr.method[i] <- corr.results[[i]]$corr.method
        corrs[[i]] <- corr.results[[i]]$corrs
        corr.files[[i]] <- corr.results[[i]]$corr.files
      }
    } else {
      corr.names <- corr.results$corr.names
      y.var <- corr.results$y.var
      x.var <- corr.results$x.var
      corr.method <- corr.results$corr.method
      corrs <- list(corr.results$corrs)
      corr.files <- list(corr.results$corr.files)
    }
    names(corrs) <- names(corr.files) <- corr.names
    corr.num <- length(corrs)
  }
  if (!is.null(meta)) {
    design <- long <- sample.id <- subject.id <- control.var <- control.val <- 
      baseline.var <- baseline.val <- time.var <- data.type <- list()
    for (i in 1:length(meta)) {
      if (meta[[i]]$data.type %in% c("microarray", "rnaseq")) {
        exprs <- meta[[i]]$y
      }
      design[[i]] <- meta[[i]]$design
      if (is.null(meta[[i]]$time.var)) {
        time.var[i] <- list(NULL)
      } else {
        time.var[[i]] <- meta[[i]]$time.var
      }
      if (is.null(meta[[i]]$control.var)) {
        control.var[i] <- list(NULL)
      } else {
        control.var[[i]] <- meta[[i]]$control.var
      }
      if (is.null(meta[[i]]$control.val)) {
        control.val[i] <- list(NULL)
      } else {
        control.val[[i]] <- meta[[i]]$control.val
      }
      if (is.null(meta[[i]]$baseline.var)) {
        baseline.var[i] <- list(NULL)
      } else {
        baseline.var[[i]] <- meta[[i]]$baseline.var
      }
      if (is.null(meta[[i]]$baseline.val)) {
        baseline.val[i] <- list(NULL)
      } else {
        baseline.val[[i]] <- meta[[i]]$baseline.val
      }
      if (is.null(meta[[i]]$sample.id)) {
        sample.id[i] <- list(NULL)
      } else {
        sample.id[[i]] <- meta[[i]]$sample.id
      }
      if (is.null(meta[[i]]$subject.id)) {
        subject.id[i] <- list(NULL)
      } else {
        subject.id[[i]] <- meta[[i]]$subject.id
      }
      data.type[[i]] <- meta[[i]]$data.type
      long[[i]] <- meta[[i]]$long
    }
    names(design) <- names(long) <- names(sample.id) <- names(subject.id) <- 
      names(control.var) <- names(control.val) <- names(baseline.var) <- 
      names(baseline.val) <- names(time.var) <- unlist(data.type)
    if (all(! names(design) %in% c("microarray", "rnaseq"))) {
      design$rnaseq <- long$rnaseq <- sample.id$rnaseq <- subject.id$rnaseq <- 
        control.var$rnaseq <- control.val$rnaseq <- baseline.var$rnaseq <- 
        baseline.val$rnaseq <- time.var$rnaseq <- list(NULL)
    }
    if (all(! names(design) %in% "flow")) {
      design$flow <- long$flow <- sample.id$flow <- subject.id$flow <- 
        control.var$flow <- control.val$flow <- baseline.var$flow <- 
        baseline.val$flow <- time.var$flow <- list(NULL)
    }
    if(all(! names(design) %in% "metab")) {
      design$metab <- long$metab <- sample.id$metab <- subject.id$metab <- 
        control.var$metab <- control.val$metab <- baseline.var$metab <- 
        baseline.val$metab <- time.var$metab <- list(NULL)
    }
  }
  if (!is.null(model.results)) {
    for (i in 1:length(model.results)) {
      if (model.results[[i]]$data.type %in% c("microarray", "rnaseq")) {
        index <- grep("DF.", colnames(model.results[[i]]$results))
        if (length(index) > 0){
          results.file <- model.results[[i]]$results[, -index]
        } else {
          results.file <- model.results[[i]]$results
        }
        dge.gsets <- model.results[[i]]$gene.sets
        dge.annots <- model.results[[i]]$annotations
      }
      if (model.results[[i]]$data.type == "flow") {
        flow.results <- model.results[[i]]$results
        dat <- data.frame(t(model.results[[i]]$y))
        rownames(dat) <- NULL
        flow.data <- cbind(design[[which(names(design) == "flow")]], dat)
      }
      if (model.results[[i]]$data.type == "metab") {
        metab.results <- model.results[[i]]$results
        dat <- data.frame(t(model.results[[i]]$y))
        rownames(dat) <- NULL
        metab.data <- cbind(design[[which(names(design) == "metab")]], dat)
      }
    }
  }
  if (!is.null(folder.path)) {
    if (!dir.exists(paste0(folder.path, "/", project.name, "/"))) {
      dir.create(paste0(folder.path, "/", project.name, "/"))
    }
    path <- file.path(paste0(folder.path, "/", project.name, "/"))
  } else {
    if (!dir.exists(paste(getwd(), "/", project.name, "/", sep = ""))) {
      dir.create(paste(getwd(), "/", project.name, "/", sep = ""))
    } 
    path <- file.path(paste0(getwd(), "/", project.name, "/"))
  }
  save(exprs, design, scores.base, scores.ctrl, modules, rowdend1b, rowdend2b, 
       rowdend1, rowdend2, rowdend3, norm.method, dist.method, agg.method, 
       time.var, control.var, control.val, baseline.var, baseline.val, 
       sample.id, subject.id, results.file, dge.gsets, dge.annots, 
       qusage.results, lower.ci, upper.ci, gene.sets, annots, roast.results, 
       flow.results, flow.data, metab.results, metab.data, corr.num, corr.names, 
       x.var, y.var, corr.method, corrs, corr.files, project.name, 
       file = paste0(path, "/bartResults.rda"))
}

#' @rdname genFile
#' @param load.path Folder path of BART file that needs to be updated.
#' @param output.path Folder path to save updated BART file. If NULL (default),
#'   updated BART file will override the old one.
#' @details \code{updateFile} takes an existing BART file and allows the user to
#'   easily add and/or revise components. The user can override the old BART
#'   file by leaving \code{output.path} as NULL or give a new path in which to
#'   save the updated file.
#' @export
updateFile <- function(load.path = NULL, output.path = NULL, meta = NULL, 
                       model.results = NULL, module.scores = NULL, 
                       dendrograms = NULL, qusage.results = NULL, 
                       roast.results = NULL, corr.results = NULL, 
                       project.name = NULL) {
  if (!is.null(load.path)) {
    if (file.exists(load.path)) {
      qresults <- qusage.results
      rst.results <- roast.results
      modules <- norm.method <- dist.method <- agg.method <- gene.sets <- 
        annots <- NULL
      data <- load(paste(load.path,"/bartResults.rda",sep = ""))
      if (!is.null(module.scores)) {
        scores.base <- module.scores$scores.base
        scores.ctrl <- module.scores$scores.ctrl
      }
      if (!is.null(dendrograms)) {
        rowdend1b <- dendrograms$rowdend1b
        rowdend2b <- dendrograms$rowdend2b
        rowdend1 <- dendrograms$rowdend1
        rowdend2 <- dendrograms$rowdend2
        rowdend3 <- dendrograms$rowdend3
        norm.method <- dendrograms$norm.method
        dist.method <- dendrograms$dist.method
        agg.method <- dendrograms$agg.method
      }
      if (!is.null(qresults)) {
        lower.ci <- qresults$lower.ci
        upper.ci <- qresults$upper.ci
        qusage.results <- qresults$qusage.results
        gene.sets <- qresults$gene.sets
        annots <- qresults$annotations
      }
      if (!is.null(rst.results)) {
        roast.results <- rst.results
      }
      if (!is.null(corr.results)) {
        if (is.list(corr.results)) {
          corr.names <- y.var <- x.var <- corr.method <- c()
          corrs <- corr.files <- list()
          for (i in 1:length(corr.results)) {
            corr.names[i] <- corr.results[[i]]$corr.names
            y.var[i] <- corr.results[[i]]$y.var
            x.var[i] <- corr.results[[i]]$x.var
            corr.method[i] <- corr.results[[i]]$corr.method
            corrs[[i]] <- corr.results[[i]]$corrs
            corr.files[[i]] <- corr.results[[i]]$corr.files
          }
        } else {
          corr.names <- corr.results$corr.names
          y.var <- corr.results$y.var
          x.var <- corr.results$x.var
          corr.method <- corr.results$corr.method
          corrs <- list(corr.results$corrs)
          corr.files <- list(corr.results$corr.files)
        }
        names(corrs) <- names(corr.files) <- corr.names
        corr.num <- length(corrs)
      }
      if (!is.null(model.results)) {
        for (i in 1:length(model.results)) {
          if (model.results[[i]]$data.type %in% c("microarray", "rnaseq")) {
            index <- grep("DF.", colnames(model.results[[i]]$results))
            if (length(index) > 0){
              results.file <- model.results[[i]]$results[, -index]
            } else {
              results.file <- model.results[[i]]$results
            }
            dge.gsets <- model.results[[i]]$gene.sets
            dge.annots <- model.results[[i]]$annotations
          }
          if (model.results[[i]]$data.type == "flow") {
            flow.results <- model.results[[i]]$results
            dat <- data.frame(t(model.results[[i]]$y))
            rownames(dat) <- NULL
            flow.data <- cbind(design[[which(names(design) == "flow")]], dat)
          }
          if (model.results[[i]]$data.type == "metab") {
            metab.results <- model.results[[i]]$results
            dat <- data.frame(t(model.results[[i]]$y))
            rownames(dat) <- NULL
            metab.data <- cbind(design[[which(names(design) == "metab")]], dat)
          }
        }
      }
      if (!is.null(meta)) {
        design <- long <- sample.id <- subject.id <- control.var <- 
          control.val <- baseline.var <- baseline.val <- time.var <- 
          data.type <- list()
        for (i in 1:length(meta)) {
          if (meta[[i]]$data.type %in% c("microarray", "rnaseq")) {
            exprs <- meta[[i]]$y
          }
          design[[i]] <- meta[[i]]$design
          if (is.null(meta[[i]]$time.var)) {
            time.var[i] <- list(NULL)
          } else {
            time.var[[i]] <- meta[[i]]$time.var
          }
          if (is.null(meta[[i]]$control.var)) {
            control.var[i] <- list(NULL)
          } else {
            control.var[[i]] <- meta[[i]]$control.var
          }
          if (is.null(meta[[i]]$control.val)) {
            control.val[i] <- list(NULL)
          } else {
            control.val[[i]] <- meta[[i]]$control.val
          }
          if (is.null(meta[[i]]$baseline.var)) {
            baseline.var[i] <- list(NULL)
          } else {
            baseline.var[[i]] <- meta[[i]]$baseline.var
          }
          if (is.null(meta[[i]]$baseline.val)) {
            baseline.val[i] <- list(NULL)
          } else {
            baseline.val[[i]] <- meta[[i]]$baseline.val
          }
          if (is.null(meta[[i]]$sample.id)) {
            sample.id[i] <- list(NULL)
          } else {
            sample.id[[i]] <- meta[[i]]$sample.id
          }
          if (is.null(meta[[i]]$subject.id)) {
            subject.id[i] <- list(NULL)
          } else {
            subject.id[[i]] <- meta[[i]]$subject.id
          }
          data.type[[i]] <- meta[[i]]$data.type
          long[[i]] <- meta[[i]]$long
        }
        names(design) <- names(long) <- names(sample.id) <- names(subject.id) <- 
          names(control.var) <- names(control.val) <- names(baseline.var) <- 
          names(baseline.val) <- names(time.var) <- unlist(data.type)
        if (all(! names(design) %in% c("microarray", "rnaseq"))) {
          design$rnaseq <- long$rnaseq <- sample.id$rnaseq <- 
            subject.id$rnaseq <- control.var$rnaseq <- control.val$rnaseq <- 
            baseline.var$rnaseq <- baseline.val$rnaseq <- time.var$rnaseq <- 
            list(NULL)
        }
        if (all(! names(design) %in% "flow")) {
          design$flow <- long$flow <- sample.id$flow <- subject.id$flow <- 
            control.var$flow <- control.val$flow <- baseline.var$flow <- 
            baseline.val$flow <- time.var$flow <- list(NULL)
        }
        if(all(! names(design) %in% "metab")) {
          design$metab <- long$metab <- sample.id$metab <- subject.id$metab <- 
            control.var$metab <- control.val$metab <- baseline.var$metab <- 
            baseline.val$metab <- time.var$metab <- list(NULL)
        }
      }
      if(is.null(output.path)){
        path <- file.path(paste0(load.path, "/bartResults.rda"))
      } else {
        dir.create(paste(output.path, "/", project.name, "/", sep = ""))
        path <- file.path(paste0(output.path, "/", project.name, 
                                 "/bartResults.rda")
        )
      }
      save(exprs, design, scores.base, scores.ctrl, modules, rowdend1b, 
           rowdend2b, rowdend1, rowdend2, rowdend3, norm.method, dist.method, 
           agg.method, time.var, control.var, control.val, baseline.var, 
           baseline.val, sample.id, subject.id, results.file, dge.gsets, 
           dge.annots, qusage.results, lower.ci, upper.ci, gene.sets, annots, 
           roast.results, flow.results, flow.data, metab.results, metab.data, 
           corr.num, corr.names, x.var, y.var, corr.method, corrs, corr.files, 
           project.name, file = path)
    } else {
      return(print("File path or data file does not exist"))
    }
  } else {
    return(print("Please enter a valid existing BART file path, or create a new 
                 BART file using genFile"))
  }
}