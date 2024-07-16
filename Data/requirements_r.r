if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
packages <- c("GSEABase", "GSVA", "limma", "msigdbr", "ggpubr", 
              "stringr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer", 
              "ggnewscale", "cowplot", "ComplexHeatmap", "patchwork", 
              "cluster", "BiocManager")

# Checking and installing required packages
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  install.packages(new_packages)
}

# Bioconductor packages
bioconductor_packages <- c("GSEABase", "GSVA", "limma", "msigdbr", "ComplexHeatmap")
installed_bioconductor <- BiocManager::installed()
new_bioconductor_packages <- bioconductor_packages[!(bioconductor_packages %in% installed_bioconductor)]
if(length(new_bioconductor_packages)) {
  BiocManager::install(new_bioconductor_packages)
}

lapply(packages, library, character.only = TRUE)
