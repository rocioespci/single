# SINGLe 
SNIPs In Nanopore reads of Gene Libraries

Full code for manuscript Espada et al 2020

This git includes:
- R package for reproducing or adapting the analysis proposed (single.Rproj)
- The main code in R (using the package) to repeat the analysis or adapt it to your data (pipelines/main_SINGLe.r)
- Clustering and consensus downstream can be done by editing and running pipelines/main_clustering_and_consensus.r.

## Installation
The easiest way to install this package is to run in an R console:
require(devtools)
install_github("rocioespci/single")

## Run
Open Pipelines/main_SINGLe.r, edit the inputs and run as any R code. 
Open Pipelines/main_Clustering_Consensus.r, edit the inputs and run as any R code. 

