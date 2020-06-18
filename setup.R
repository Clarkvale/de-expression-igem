

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install()
BiocManager::install("GEOquery")
BiocManager::install("Biobase")
BiocManager::install("limma")
BiocManager::install("edgeR")


install.packages("devtools")
library(devtools)
devtools::install_github("aryoda/R_enumerations@v0.3.0-beta")




