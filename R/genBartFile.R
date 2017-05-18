#' Generate and Update BART files
#'
#' Generate and/or update BART file that can be directly uploaded into the BART
#' app
#' @param design_info list generated from \code{desInfo}
#' @param model_results Default is NULL. Object generated from function
#'   \code{genModelResults}
#' @param module_maps Default is NULL. Object generated from function
#'   \code{genModules}
#' @param dendros Default is NULL. Object generated from function
#'   \code{genDendrograms}
#' @param qusage_results Default is NULL. Object generated from function
#'   \code{qBart}
#' @param roast_results Default is NULL. Object generated from function
#'   \code{rBart}
#' @param corr_results Default is NULL. Object generated from function
#'   \code{crossCorr}
#' @param illumina Logical. Are the probe IDs Illumina? Default is TRUE.
#' @details \code{genFile} creates a formatted BART file (Unsupervised.RData)
#'   based on input components that the user provides. The function then saves
#'   the file in a folder named after the project name that was defined in the
#'   function \code{\link{desInfo}}. The folder is created in the current
#'   workign directory. All function inputs are the objects generated from
#'   running the other available functions in \code{genBART}.
#' @examples
#' # Example data
#' data(tb.expr)
#' data(tb.design)
#' 
#' # Use first 100 probes to demonstrate
#' dat <- tb.expr[1:100,]
#' 
#' # Create desInfo object
#' des.info <- desInfo(y = dat, design = tb.design, data_type = "micro", 
#'                     columnname = "columnname", long = TRUE, patient_id = "monkey_id",
#'                     baseline_var = "timepoint", baseline_val = 0, time_var = "timepoint", 
#'                     responder_var = "clinical_status", sample_id = "sample_id", 
#'                     project_name = "TB")
#'                     
#' # create BART file (minimal example)
#' genFile(design_info = list(des.info))                     
#'
#' # generate module maps and cluster probes
#' mods <- genModules(des.info, modules)
#' dendros <- genDendrograms(des.info)
#'
#' # Update BART file with module maps and heat maps
#' path <- paste(getwd(), "/", des.info$project_name, " Pipeline", sep = "")
#' updateFile(load.path = path, module_maps = mods, dendros = dendros)
#' @export
genFile <- function(design_info, model_results = NULL, module_maps = NULL,
                    dendros = NULL, qusage_results = NULL,
                    roast_results = NULL, corr_results = NULL,
                    illumina = TRUE) {
  if (!is.null(model_results)) {
    for (i in 1:length(model_results)) {
      if (model_results[[i]]$data_type %in% c("micro", "rna")) {
        index <- grep("DF.", colnames(model_results[[i]]$results))
        if (length(index) > 0){
          results_file <- model_results[[i]]$results[, -index]
        } else {
          results_file <- model_results[[i]]$results
        }
      }
    }
  }
  design <- NULL
  long <- NULL
  mod1 <- NULL
  mod2 <- NULL
  base_mod <- NULL
  long_mod <- NULL
  long_mod2 <- NULL
  h1b_coldendro <- NULL
  h1b_rowdendro <- NULL
  h2b_coldendro <- NULL
  h2b_rowdendro <- NULL
  h1_coldendro <- NULL
  h1_rowdendro <- NULL
  h2_coldendro <- NULL
  h2_rowdendro <- NULL
  h3_coldendro <- NULL
  h3_rowdendro <- NULL
  lowerCI <- NULL
  upperCI <- NULL
  correlation_files_count <- NULL
  correlation_names <- NULL
  y_var <- NULL
  x_var <- NULL
  correlation_method <- NULL
  correlations <- NULL
  correlation_files <- NULL
  f_design <- NULL
  f_summary_var <- NULL
  f_time_var <- NULL
  f_responder_var <- NULL
  f_control_var <- NULL
  f_control_val <- NULL
  f_baseline_var <- NULL
  f_baseline_val <- NULL
  f_sample_id <- NULL
  f_patient_id <- NULL
  f_hc <- NULL
  m_design <- NULL
  m_summary_var <- NULL
  m_time_var <- NULL
  m_responder_var <- NULL
  m_control_var <- NULL
  m_control_val <- NULL
  m_baseline_var <- NULL
  m_baseline_val <- NULL
  m_sample_id <- NULL
  m_patient_id <- NULL
  m_hc <- NULL
  summary_var <- NULL
  time_var <- NULL
  responder_var <- NULL
  control_var <- NULL
  control_val <- NULL
  baseline_var <- NULL
  baseline_val <- NULL
  sample_id <- NULL
  patient_id <- NULL
  hc <- NULL
  flow_results <- NULL
  flow_data <- NULL
  metab_results <- NULL
  metab_data <- NULL
  project_name <- NULL
  data_type <- NULL
  if (!is.null(module_maps)) {
    base_mod <- module_maps$base_mod
    long_mod <- module_maps$long_mod
    long_mod2 <- module_maps$long_mod2
  }
  if (!is.null(dendros)) {
    h1b_coldendro <- dendros$h1b_coldendro
    h1b_rowdendro <- dendros$h1b_rowdendro
    h2b_coldendro <- dendros$h2b_coldendro
    h2b_rowdendro <- dendros$h2b_rowdendro
    h1_coldendro <- dendros$h1_coldendro
    h1_rowdendro <- dendros$h1_rowdendro
    h2_coldendro <- dendros$h2_coldendro
    h2_rowdendro <- dendros$h2_rowdendro
    h3_coldendro <- dendros$h3_coldendro
    h3_rowdendro <- dendros$h3_rowdendro
  }
  if (!is.null(qusage_results)) {
    lowerCI <- qusage_results$lowerCI
    upperCI <- qusage_results$upperCI
    qusage_results <- qusage_results$qusage_results
  }
  if (!is.null(corr_results)) {
    if (is.list(corr_results)) {
      correlation_names <- y_var <- x_var <- correlation_method <- c()
      correlations <- correlation_files <- list()
      for (i in 1:length(corr_results)) {
        correlation_names[i] <- corr_results[[i]]$correlation_names
        y_var[i] <- corr_results[[i]]$y_var
        x_var[i] <- corr_results[[i]]$x_var
        correlation_method[i] <- corr_results[[i]]$correlation_method
        correlations[[i]] <- corr_results[[i]]$correlations
        correlation_files[[i]] <- corr_results[[i]]$correlation_files
      }
    } else {
      correlation_names <- corr_results$correlation_names
      y_var <- corr_results$y_var
      x_var <- corr_results$x_var
      correlation_method <- corr_results$correlation_method
      correlations <- list(corr_results$correlations)
      correlation_files <- list(corr_results$correlation_files)
    }
    names(correlations) <- names(correlation_files) <- correlation_names
    correlation_files_count <- length(correlations)
  }
  if (!is.null(design_info)) {
    for (i in 1:length(design_info)) {
      if (design_info[[i]]$data_type %in% c("micro", "rna")) {
        PROBE_ID <- SYMBOL <- rownames(design_info[[i]]$y)
        final_expression <- cbind(PROBE_ID = PROBE_ID, SYMBOL = SYMBOL,
                                  design_info[[i]]$y)
        design <- design_info[[i]]$design
        summary_var <- design_info[[i]]$summary_var
        time_var <- design_info[[i]]$time_var
        responder_var <- design_info[[i]]$responder_var
        control_var <- design_info[[i]]$control_var
        control_val <- design_info[[i]]$control_val
        baseline_var <- design_info[[i]]$baseline_var
        baseline_val <- design_info[[i]]$baseline_val
        sample_id <- design_info[[i]]$sample_id
        patient_id <- design_info[[i]]$patient_id
        hc <- design_info[[i]]$hc
        project_name <- design_info[[i]]$project_name
        data_type <- design_info[[i]]$data_type
        long <- design_info[[i]]$long
      }
      if (design_info[[i]]$data_type == "flow") {
        f_design <- design_info[[i]]$design
        f_summary_var <- design_info[[i]]$summary_var
        f_time_var <- design_info[[i]]$time_var
        f_responder_var <- design_info[[i]]$responder_var
        f_control_var <- design_info[[i]]$control_var
        f_control_val <- design_info[[i]]$control_val
        f_baseline_var <- design_info[[i]]$baseline_var
        f_baseline_val <- design_info[[i]]$baseline_val
        f_sample_id <- design_info[[i]]$sample_id
        f_patient_id <- design_info[[i]]$patient_id
        f_hc <- design_info[[i]]$hc
      }
      if (design_info[[i]]$data_type == "metab") {
        m_design <- design_info[[i]]$design
        m_summary_var <- design_info[[i]]$summary_var
        m_time_var <- design_info[[i]]$time_var
        m_responder_var <- design_info[[i]]$responder_var
        m_control_var <- design_info[[i]]$control_var
        m_control_val <- design_info[[i]]$control_val
        m_baseline_var <- design_info[[i]]$baseline_var
        m_baseline_val <- design_info[[i]]$baseline_val
        m_sample_id <- design_info[[i]]$sample_id
        m_patient_id <- design_info[[i]]$patient_id
        m_hc <- design_info[[i]]$hc
      }
    }
  }
  if (!is.null(model_results)) {
    for (i in 1:length(model_results)) {
      if (model_results[[i]]$data_type == "flow") {
        flow_results <- model_results[[i]]$results
        dat <- data.frame(t(model_results[[i]]$y))
        rownames(dat) <- NULL
        flow_data <- cbind(f_design, dat)
      }
      if (model_results[[i]]$data_type == "metab") {
        metab_results <- model_results[[i]]$results
        dat <- data.frame(t(model_results[[i]]$y))
        rownames(dat) <- NULL
        metab_data <- cbind(m_design, dat)
      }
    }
  }
  if (!exists("final_expression")) {
    final_expression <- NULL
  }
  if (!exists("results_file")) {
    results_file <- NULL
  }
  illumina <- illumina
  dir.create(paste(getwd(), "/", project_name, " Pipeline/", sep = ""))
  save(mod1, mod2, final_expression, design, base_mod, long, long_mod, long_mod2,
       h1b_rowdendro, h1b_coldendro, h2b_rowdendro, h2b_coldendro, h1_rowdendro,
       h2_rowdendro, h3_rowdendro, h1_coldendro, h2_coldendro, h3_coldendro,
       summary_var, time_var, responder_var, control_var, control_val,
       baseline_var, baseline_val, sample_id, patient_id, f_design,
       f_summary_var, f_time_var, f_responder_var, f_control_var, f_control_val,
       f_baseline_var, f_baseline_val, f_sample_id, f_patient_id, f_hc,
       m_design, m_summary_var, m_time_var, m_responder_var, m_control_var,
       m_control_val, m_baseline_var, m_baseline_val, m_sample_id, m_patient_id,
       m_hc, hc, results_file, project_name, lowerCI, upperCI, qusage_results,
       roast_results, flow_results, flow_data, metab_results, metab_data,
       illumina, correlation_files_count, correlation_names, y_var, x_var,
       correlation_method, correlations, correlation_files, data_type,
       file = file.path(paste(getwd(), "/", project_name, " Pipeline/",
                              "Unsupervised.Rdata", sep = "")))
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
updateFile <- function(load.path = NULL, output.path = NULL, design_info = NULL, model_results = NULL, module_maps = NULL,
                       dendros = NULL, qusage_results = NULL,
                       roast_results = NULL, corr_results = NULL,
                       illumina = NULL) {
  if (!is.null(load.path)) {
    if (file.exists(load.path)) {
      qresults <- qusage_results
      mod1 <- mod2 <- long <- data_type <- NULL
      data <- load(paste(load.path,"/Unsupervised.RData",sep = ""))
      if (!is.null(module_maps)) {
        base_mod <- module_maps$base_mod
        long_mod <- module_maps$long_mod
        long_mod2 <- module_maps$long_mod2
      }
      if (!is.null(dendros)) {
        h1b_coldendro <- dendros$h1b_coldendro
        h1b_rowdendro <- dendros$h1b_rowdendro
        h2b_coldendro <- dendros$h2b_coldendro
        h2b_rowdendro <- dendros$h2b_rowdendro
        h1_coldendro <- dendros$h1_coldendro
        h1_rowdendro <- dendros$h1_rowdendro
        h2_coldendro <- dendros$h2_coldendro
        h2_rowdendro <- dendros$h2_rowdendro
        h3_coldendro <- dendros$h3_coldendro
        h3_rowdendro <- dendros$h3_rowdendro
      }
      if (!is.null(qresults)) {
        lowerCI <- qresults$lowerCI
        upperCI <- qresults$upperCI
        qusage_results <- qresults$qusage_results
      }
      if (!is.null(corr_results)) {
        if (is.list(corr_results)) {
          correlation_names <- y_var <- x_var <- correlation_method <- c()
          correlations <- correlation_files <- list()
          for (i in 1:length(corr_results)) {
            correlation_names[i] <- corr_results[[i]]$correlation_names
            y_var[i] <- corr_results[[i]]$y_var
            x_var[i] <- corr_results[[i]]$x_var
            correlation_method[i] <- corr_results[[i]]$correlation_method
            correlations[[i]] <- corr_results[[i]]$correlations
            correlation_files[[i]] <- corr_results[[i]]$correlation_files
          }
        } else {
          correlation_names <- corr_results$correlation_names
          y_var <- corr_results$y_var
          x_var <- corr_results$x_var
          correlation_method <- corr_results$correlation_method
          correlations <- list(corr_results$correlations)
          correlation_files <- list(corr_results$correlation_files)
        }
        names(correlations) <- names(correlation_files) <- correlation_names
        correlation_files_count <- length(correlations)
      }
      if (!is.null(model_results)) {
        for (i in 1:length(model_results)) {
          if (model_results[[i]]$data_type %in% c("micro", "rna")) {
            index <- grep("DF.", colnames(model_results[[i]]$results))
            if (length(index) > 0){
              results_file <- model_results[[i]]$results[, -index]
            } else {
              results_file <- model_results[[i]]$results
            }
          }
          if (model_results[[i]]$data_type == "flow") {
            flow_results <- model_results[[i]]$results
            dat <- data.frame(t(model_results[[i]]$y))
            rownames(dat) <- NULL
            flow_data <- cbind(f_design, dat)
          }
          if (model_results[[i]]$data_type == "metab") {
            metab_results <- model_results[[i]]$results
            dat <- data.frame(t(model_results[[i]]$y))
            rownames(dat) <- NULL
            metab_data <- cbind(m_design, dat)
          }
        }
      }
      if (!is.null(design_info)) {
        for (i in 1:length(design_info)) {
          if (design_info[[i]]$data_type %in% c("micro", "rna")) {
            PROBE_ID <- SYMBOL <- rownames(design_info[[i]]$y)
            final_expression <- cbind(PROBE_ID = PROBE_ID, SYMBOL = SYMBOL,
                                      design_info[[i]]$y)
            design <- design_info[[i]]$design
            summary_var <- design_info[[i]]$summary_var
            time_var <- design_info[[i]]$time_var
            responder_var <- design_info[[i]]$responder_var
            control_var <- design_info[[i]]$control_var
            control_val <- design_info[[i]]$control_val
            baseline_var <- design_info[[i]]$baseline_var
            baseline_val <- design_info[[i]]$baseline_val
            sample_id <- design_info[[i]]$sample_id
            patient_id <- design_info[[i]]$patient_id
            hc <- design_info[[i]]$hc
            project_name <- design_info[[i]]$project_name
            data_type <- design_info[[i]]$data_type
            long <- design_info[[i]]$long
          }
          if (design_info[[i]]$data_type == "flow") {
            f_design <- design_info[[i]]$design
            f_summary_var <- design_info[[i]]$summary_var
            f_time_var <- design_info[[i]]$time_var
            f_responder_var <- design_info[[i]]$responder_var
            f_control_var <- design_info[[i]]$control_var
            f_control_val <- design_info[[i]]$control_val
            f_baseline_var <- design_info[[i]]$baseline_var
            f_baseline_val <- design_info[[i]]$baseline_val
            f_sample_id <- design_info[[i]]$sample_id
            f_patient_id <- design_info[[i]]$patient_id
            f_hc <- design_info[[i]]$hc
          }
          if (design_info[[i]]$data_type == "metab") {
            m_design <- design_info[[i]]$design
            m_summary_var <- design_info[[i]]$summary_var
            m_time_var <- design_info[[i]]$time_var
            m_responder_var <- design_info[[i]]$responder_var
            m_control_var <- design_info[[i]]$control_var
            m_control_val <- design_info[[i]]$control_val
            m_baseline_var <- design_info[[i]]$baseline_var
            m_baseline_val <- design_info[[i]]$baseline_val
            m_sample_id <- design_info[[i]]$sample_id
            m_patient_id <- design_info[[i]]$patient_id
            m_hc <- design_info[[i]]$hc
          }
        }
      }
      if(is.null(output.path)){
        setwd(load.path)
        save(mod1, mod2, final_expression, design, base_mod, long_mod, long, long_mod2,
             h1b_rowdendro, h1b_coldendro, h2b_rowdendro, h2b_coldendro, h1_rowdendro,
             h2_rowdendro, h3_rowdendro, h1_coldendro, h2_coldendro, h3_coldendro,
             summary_var, time_var, responder_var, control_var, control_val,
             baseline_var, baseline_val, sample_id, patient_id, f_design,
             f_summary_var, f_time_var, f_responder_var, f_control_var, f_control_val,
             f_baseline_var, f_baseline_val, f_sample_id, f_patient_id, f_hc,
             m_design, m_summary_var, m_time_var, m_responder_var, m_control_var,
             m_control_val, m_baseline_var, m_baseline_val, m_sample_id, m_patient_id,
             m_hc, hc, results_file, project_name, lowerCI, upperCI, qusage_results,
             roast_results, flow_results, flow_data, metab_results, metab_data,
             illumina, correlation_files_count, correlation_names, y_var, x_var,
             correlation_method, correlations, correlation_files, data_type,
             file = "Unsupervised.RData")
      } else {
        setwd(output.path)
        dir.create(paste(getwd(), "/", project_name, " Pipeline/", sep = ""))
        save(mod1, mod2, final_expression, design, base_mod, long_mod, long, long_mod2,
             h1b_rowdendro, h1b_coldendro, h2b_rowdendro, h2b_coldendro, h1_rowdendro,
             h2_rowdendro, h3_rowdendro, h1_coldendro, h2_coldendro, h3_coldendro,
             summary_var, time_var, responder_var, control_var, control_val,
             baseline_var, baseline_val, sample_id, patient_id, f_design,
             f_summary_var, f_time_var, f_responder_var, f_control_var, f_control_val,
             f_baseline_var, f_baseline_val, f_sample_id, f_patient_id, f_hc,
             m_design, m_summary_var, m_time_var, m_responder_var, m_control_var,
             m_control_val, m_baseline_var, m_baseline_val, m_sample_id, m_patient_id,
             m_hc, hc, results_file, project_name, lowerCI, upperCI, qusage_results,
             roast_results, flow_results, flow_data, metab_results, metab_data,
             illumina, correlation_files_count, correlation_names, y_var, x_var,
             correlation_method, correlations, correlation_files, data_type,
             file = file.path(paste(getwd(), "/", project_name, " Pipeline/",
                                    "Unsupervised.Rdata", sep = "")))
      }
    } else {
      return(print("File path or data file does not exist"))
    }
  } else {
    return(print("Please enter a valid existing BART file path, or create a new BART file using genFile"))
  }
}






