## Install packages required for local instance of Strongyloides RNA-seq Browser
# Note, if you're on a Windows computer, you may need to change the folder your R packages are stored in from Read-only.
setRepositories(ind = c(1,2,3,4,5,6))
if (!requireNamespace("pacman", quietly = TRUE))
    install.packages("pacman")
library(pacman)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
library(BiocManager)

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

BiocManager::install("DO.db", ask=FALSE)
BiocManager::install("HDO.db", ask=FALSE)
BiocManager::install("igraph", ask=FALSE)
BiocManager::install("clusterProfiler", ask=FALSE)
BiocManager::install(c("limma","edgeR"),ask=FALSE)

pacman::p_load(shiny,shinyjs,shinyWidgets,htmltools,shinythemes,DT,Cairo,tidyverse,ggforce,gt,plotly,magrittr,ggthemes,gplots,svglite,Cairo,heatmaply,RColorBrewer,rcartocolor,openxlsx,egg,dendextend,vctrs,markdown,config)
# Note, on mac you may need to download XQuartz for some of these packages, specially Cairo to work correctly. 

libPath <- .libPaths()
remove.packages("patchwork", lib=libPath)
devtools::install_github("thomasp85/patchwork")