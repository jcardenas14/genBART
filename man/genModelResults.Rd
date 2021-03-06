% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genModelResults.R
\name{genModelResults}
\alias{genModelResults}
\title{Generate formatted results file}
\usage{
genModelResults(y = NULL, data.type = "rnaseq", method = "limma", object,
  lm.Fit = NULL, comp.names = NULL, var.symbols = NULL,
  gene.sets = NULL, annotations = NULL)
}
\arguments{
\item{y}{Expression data frame the model was run on. Default is NULL. See 
Details for further description.}

\item{data.type}{Type of data being analyzed ("rnaseq", "microarray", "flow", 
"metab"). Default is data.type="rnaseq".}

\item{method}{string denoting modeling method used ("limma","deseq2","edgeR")}

\item{object}{results object generated from \code{\link[limma]{eBayes}} 
(limma), \code{\link[DESeq2]{results}} (DESeq2), 
\code{\link[edgeR]{glmLRT}} or \code{\link[edgeR]{glmQLFTest}} (edgeR). 
Objects generated from \code{results}, \code{glmLRT}, or \code{glmlQLFTest} 
must be wrapped in a list. See Details for further description.}

\item{lm.Fit}{linear model fit object generated from 
\code{\link[limma]{lmFit}} (limma). This parameter can can be left as NULL 
when \code{method} equals "deseq2" or "edgeR".}

\item{comp.names}{Optional vector of comparison (contrast) names. Default 
is NULL.}

\item{var.symbols}{Optional vector of additional annotations for the 
variables. Otherwise, the rownames of the expression data are used. Default 
is NULL.}

\item{gene.sets}{A list of gene sets. Default is NULL. See Details for 
further description.}

\item{annotations}{A data frame of additional annotations for the gene sets. 
Default is NULL. See Details for further description.}
}
\value{
\code{data.type} string denoting the type of data that was analyzed

\code{results} the formatted results returned as a data frame

\code{resids} data frame of residuals. Returned only if 
  data.type="microarray" or "rnaseq" and method="limma". Used to estimate the 
  VIFs when running the Qusage algorithm in \code{\link{runQgen}}.

\code{gene.sets} list of gene sets provided by the user. NULL if no 
  list provided.

\code{annotations} data frame of gene set annotations provided by the
  user. Null if no annotations are provided.
}
\description{
Generate formatted results file from result objects returned by limma, 
DESeq2, and edgeR pipelines
}
\details{
This function takes results obtained from differential analysis 
  pipelines found in \code{limma}, \code{DESeq2}, or \code{edgeR} and formats 
  them for the BART app. 
  
  The expression data \code{y} and \code{lm.Fit} objects are used to obtain 
  the residual matrix from the fitted model. These parameters are only needed 
  when method = "limma" and data.type = "rnaseq" or "microarray" and can 
  otherwise be left as NULL. It is important to remember that \code{y} should 
  be the expression data used for modeling (e.g. voom transformed data). The 
  residual matrix is stored as an element of the returned list and can be 
  used in downstream gene set analysis using \code{\link{runQgen}} (Please 
  visit for more details).
  
  The \code{object} parameter takes as input model result objects returned by 
  functions in limma, DESeq2, or edgeR. When method = "limma", the expected 
  input is the single object returned by \code{\link[limma]{eBayes}} since it 
  is able to store results across multiple comparisons. When method = 
  "deseq2" or "edgeR", the result object(s) returned by 
  \code{\link[DESeq2]{results}}, \code{\link[edgeR]{glmLRT}}, or 
  \code{\link[edgeR]{glmQLFTest}} must be wrapped in a list in which each 
  element is an object containing the results for a single comparison. 
  
  The \code{comp.names} parameter is a character vector of comparison names 
  that is particularly useful when method = "deseq2" or "edgeR" since 
  comparison names are not extracted from the result objects generated by 
  either of those pipelines. When using limma, the comparison names can also 
  be defined in \code{\link[limma]{makeContrasts}}. It is important that the 
  names are written in the correct order. For example, if object = list(AvsB, 
  CvsD), where AvsB and CvsD are result objects for the comparisons "group A 
  vs group B" and "group C vs group D" respectively, then comp.names = 
  c("CvsD", "AvsB") would incorrectly assign the name "CvsD" to the 
  comparison "group A vs group B" and vice versa. The \code{var.symbols} 
  parameter is typically used to provide a character vector of gene symbols. 
  The vector provided must be the same length and in the same order as the 
  row names of the data used for modeling.
  
  The \code{gene.sets} parameter is a list in which each element is a 
  character vector of gene names comprising a gene set. The gene names must 
  match the rownames of the data used for modeling. The gene sets are used to 
  create modular maps for each comparison in the DGE section of BART. The 
  \code{annotations} parameter is a data frame consisting of two columns. The 
  first column consists of gene set names and the second column consists of 
  additional descriptions for the gene sets.
}
\examples{
# Example data
data(tb.expr)
data(tb.design)

# Only use first 100 genes to demonstrate
dat <- tb.expr[1:100,]

# Generate lmFit and eBayes (limma) objects needed for genModelResults
tb.design$Group <- paste(tb.design$clinical_status,tb.design$timepoint,sep = "")
grp <- factor(tb.design$Group)
design2 <- model.matrix(~0+grp)
colnames(design2) <- levels(grp)

dupcor <- limma::duplicateCorrelation(dat, design2, block = tb.design$monkey_id)
fit <- limma::lmFit(dat, design2, block = tb.design$monkey_id, 
             correlation = dupcor$consensus.correlation)
contrasts <- limma::makeContrasts(A_20vsPre = Active20-Active0, A_42vsPre = Active42-Active0, 
                                  levels=design2)
fit2 <- limma::contrasts.fit(fit, contrasts)
fit2 <- limma::eBayes(fit2, trend = FALSE)

# Format results
model.results <- genModelResults(y = dat, data.type = "microarray", object = fit2,
                                 lm.Fit = fit, method = "limma")
}
