\name{plot_tree}
\alias{plot_tree}
\title{plot_tree}
\usage{
plot_tree(hierarchical_clustering,ids,k_clusters, colors=rainbow(k_clusters),main,consensus_labels=NULL)
}
\description{
Plots a hierarchical clustering coloring labels by ids.
}
\arguments{
  \item{hierarchical_clustering}{Dendrogram returned by hclust.}
  \item{ids}{Numeric or character vector, each positition indicating to which cluster an element belongs, as the one returned by cutree. }
  \item{k_clusters}{Numeric. Number of clusters.}
  \item{colors}{Colors vector (pallette to be used)}
  \item{main}{Character. Title for plot.}
  \item{consensus_labels}{String vector. Labels for each cluster to be plotted below the dendrogram.}
}
\value{
Plots and returns 0.
}
\examples{
require(dendextend)
distance_matrix <- matrix(runif(100),ncol=10,nrow=10)
distance_matrix <- accomodate_dist(distance_matrix, names=NULL)
distance_matrix <- as.dist(distance_matrix)
hc <- hclust(distance_matrix)
clusters <- cutree(hc, k=2)

plot_tree(hc,clusters,2, main="example")
}
