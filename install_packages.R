# Install R dependencies (CRAN + Bioconductor)
cran <- c("shiny", "flexdashboard", "plotly", "DT", "ggplot2", "survival", "survminer",
          "data.table", "msigdbr", "rmarkdown", "testthat", "BiocManager")
install.packages(setdiff(cran, rownames(installed.packages())), repos = "https://cloud.r-project.org")
BiocManager::install("clusterProfiler", ask = FALSE, update = FALSE)
