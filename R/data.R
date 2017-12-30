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


#' Gene Symbols
#' 
#' A vector of gene symbols to match the probes ids in expression
#' 
#' @format A vector of 1000 gene symbols
#' 
#' @source Illumina
"gene.symbols"

