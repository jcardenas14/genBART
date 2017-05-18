#' Microarray design file
#' 
#' A microarray design file dataset for cynomolgus macaques infected with M. 
#' tuberculosis. Contains sample information such as time point and clinical 
#' status.
#' 
#' @format A data frame with 30 rows and 7 variables. \describe{ 
#'   \item{columnname}{names matching column names in expression dataset} 
#'   \item{monkey_id}{id assigned to annotate each monkey} 
#'   \item{timepoint}{times at which samples were drawn, in days} 
#'   \item{timerange}{groups time points into stages} \item{sample_id}{id 
#'   assigned to annotate each sample} \item{clinical_status}{denotes if disease
#'   is active or latent} }
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84152}
"tb.design"


#' Microarray data on cynomolgus macaques infected with M. tuberculosis
#' 
#' A microarray dataset on cynomolgus macaques infected with M. tuberculosis. 
#' The study was conducted on 38 monkeys at 11 unequally spaced time points. The
#' dataset presented is a subset of the full dataset. 
#' 
#' @format A data frame with 4000 rows and 30 columns, where each row is a probe
#'   and each column a sample.
#'   
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84152}
"tb.expr"


#' Flow data on cynomolgus macaques infected with M. tuberculosis
#' 
#' A dataset containing flow variables of absolute counts and percentages.
#' 
#' @format A data frame with 13 rows of flow variables and 213 columns of
#'   samples
#'   
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84152}
"tb.flow"


#' Flow design file
#' 
#' A flow design file dataset for cynomolgus macaques infected with M.
#' tuberculosis. Contains sample information such as time point and clinical
#' status.
#' 
#' @format A data frame with 213 rows of samples and 7 variables. \describe{  
#'   \item{columnname}{names matching column names in expression dataset} 
#'   \item{monkey_id}{id assigned to annotate each monkey} 
#'   \item{timepoint}{times at which samples were drawn, in days} 
#'   \item{timerange}{groups time points into stages} \item{sample_id}{id 
#'   assigned to annotate each sample} \item{clinical_status}{denotes if disease
#'   is active or latent} }
#'   
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84152}
"tb.flow.des"

#' Baylor Modules
#' 
#' A list containing the baylor modules (gene sets)
#' 
#' @format A large list containing 260 gene sets.
#' 
#' @source baylor modules
"modules"

#' Cluster Gene Sets
#' 
#' A list of clusters
#' 
#' @format A list of 10 clusters formed through hierarchical clustering of the 
#'   4000 probes described in \code{\link{tb.expr}}. The probes were clustered
#'   after centering and scaling each of the rows.
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84152}
"clusters"


#' Module Activity Scores
#' 
#' This is a dataset containing module activity scores for cynomolgus macaques 
#' infected with M. tuberculosis. For a given module, monkey, and time point, 
#' each score was calculated by summing the difference in probe-level 
#' post-infection expression values and the average basline expression value for
#' all probes within the module. This sum was then divided by the number of 
#' probes within each module plus the number of samples within each time point.
#' 
#' @format A data frame with 194 rows of samples and 260 columns of baylor
#'   modules.
#'   
#' @source baylor modules and 
#'   \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84152}
"module.as"


#' Gene Symbols
#' 
#' A vector of gene symbols to match the probes ids in expression
#' 
#' @format A vector of 1000 gene symbols
#' 
#' @source Illumina
"gene.symbols"

