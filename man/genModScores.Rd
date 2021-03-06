% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genModScores.R
\name{genModScores}
\alias{genModScores}
\title{Generate modular (gene set) maps for plotting}
\usage{
genModScores(meta, gene.sets, sd.lim = 2, annotations = NULL)
}
\arguments{
\item{meta}{list returned by \code{metaData}}

\item{gene.sets}{list of gene sets. See \code{\link{genModelResults}} for 
additional formatting details.}

\item{sd.lim}{number of standard deviations away from the mean of the 
reference samples. Default is 2.}

\item{annotations}{A data frame of additional annotations for the gene sets. 
Default is NULL. See \code{\link{genModelResults}} for additional details.}
}
\value{
\code{scores.ctrl} data frame of module scores for all samples 
  with respect to controls.

\code{scores.base} data frame of module scores for all time point 
  samples with respect to their baseline.

\code{gene.sets} List of gene sets provided through \code{gene.sets}
}
\description{
Generate modular (gene set) maps for plotting
}
\details{
This function calculates module scores for individual samples. In 
  cross sectional studies with controls, the control samples are used to 
  determine an upper and lower threshold (mean of controls +/- 2 sd). The 
  module proportion for each sample is then calculated based on the 
  percentage of genes within a module that are above or below this threshold. 
  For example, if 40\% of the genes within a module are above the threshold 
  and 15\% are below it, the final module score would be 25\% up (40-15). In 
  longitudinal settings, module scores are calculated with respect to 
  controls and baseline samples. In cross sectional studies without controls, 
  \code{genModules} cannot be used, since there are no reference samples with 
  which to calculate a threshold.
}
\examples{
# Example data
data(tb.expr)
data(tb.design)
data(modules)

# Demonstrate on first 100 probes
dat <- tb.expr[1:100, ]

# Create desInfo object
meta.data <- metaData(y = tb.expr, design = tb.design, data.type = "microarray", 
                    columnname = "columnname", long = TRUE, sample.id = "sample_id",
                    subject.id = "monkey_id", time.var = "timepoint",
                    baseline.var = "timepoint", baseline.val = 0)

# Generate module maps                                      
mods <- genModScores(meta.data, modules)
}
