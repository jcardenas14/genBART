#' Hierarchical clustering of normalized expression data
#' 
#' Perform hierarchical clustering on normalized data
#' @param norm.data list of normalized expression data returned by 
#'   \code{normalizeData}
#' @param dist.method The distance measure to be used. This must be one of 
#'   "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#'   See \code{\link[stats]{dist}} for more details.
#' @param agg.method The agglomeration method to be used. This must be one of 
#'   "single", "complete", "average", "mcquitty", "ward.D", "ward.D2", 
#'   "centroid" or "median". \code{\link[fastcluster]{hclust}} for more details.
#' @details This function performs hierarchical clustering on the rows of the 
#'   normalized expression data contained in \code{norm.data}. 
#' @return \code{rowdend1b} dendrogram from hierarchical clustering of genes on
#'   baseline samples normalized according to \code{norm.method} specified in 
#'   \code{norm.data}. NULL if \code{y1b} in \code{norm.data} is NULL.
#' @return \code{rowdend2b} dendrogram from hierarchical clustering of genes on
#'   baseline samples normalized to controls according to \code{norm.method} 
#'   specified in \code{norm.data}. NULL if \code{y2b} in \code{norm.data} is 
#'   NULL.
#' @return \code{rowdend1} dendrogram from hierarchical clustering of genes on
#'   all samples normalized according to \code{norm.method} specified in 
#'   \code{norm.data}. NULL if \code{y1} in \code{norm.data} is NULL.
#' @return \code{rowdend2} dendrogram from hierarchical clustering of genes on
#'   all samples normalized to controls according to \code{norm.method} 
#'   specified in \code{norm.data}. NULL if \code{y2} in \code{norm.data} is 
#'   NULL.
#' @return \code{rowdend3} dendrogram from hierarchical clustering of genes on
#'   all samples normalized to their baseline. NULL if \code{y3} in 
#'   \code{norm.data} is NULL.
#' @return \code{norm.method} string describing the normalization method used in
#'   \code{\link{normalizeData}}
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
#'                     columnname = "columnname", long = TRUE, sample.id = "sample_id",
#'                     subject.id = "monkey_id", time.var = "timepoint",
#'                     baseline.var = "timepoint", baseline.val = 0)
#' 
#' # Normalize data
#' data.norm <- normalizeData(meta = meta.data)
#' 
#' # Cluster data
#' dendros <- clusterData(norm.data = data.norm)
#' @importFrom stats as.dendrogram dist qt sd
#' @export
clusterData <- function(norm.data, dist.method = "euclidean", 
                        agg.method = "complete") {
  names(norm.data) <- c("rowdend1b", "rowdend2b", "rowdend1", "rowdend2", 
                        "rowdend3", "norm.method")
  rowDendros <- vector("list", (length(norm.data)-1))
  names(rowDendros) <- names(norm.data)[1:(length(norm.data)-1)]
  for (i in 1:(length(norm.data)-1)) {
    if (is.null(norm.data[[i]])) {
      rowDendros[i] <- list(NULL)
    } else {
      if (names(rowDendros)[i] %in% "rowdend1b") {
        print(paste0("clustering genes from baseline samples normalized to ", 
                     norm.data$norm.method, "..."))
      } else if (names(rowDendros)[i] %in% "rowdend2b") {
        print("clustering genes from baseline samples normalized to controls..."
        )
      } else if (names(rowDendros)[i] %in% "rowdend1") {
        print(paste0("clustering genes from all samples normalized to ", 
                     norm.data$norm.method, "..."))
      } else if (names(rowDendros)[i] %in% "rowdend2") {
        print("clustering genes from all samples normalized to controls...")
      } else if (names(rowDendros)[i] %in% "rowdend3") {
        print("clustering genes from all samples normalized to baseline...")
      }
      rowDendros[[i]] <- as.dendrogram(
        fastcluster::hclust(
          dist(norm.data[[i]], method = dist.method), method = agg.method
        )
      )
    }
  }
  rowDendros$norm.method <- norm.data$norm.method
  rowDendros$dist.method <- dist.method
  rowDendros$agg.method <- agg.method
  return(rowDendros)
}



