library(devtools)
library(BiocManager)

devtools::install_github("mtmorgan/dirichletmultinomial",
  ref = "master", repos = BiocManager::repositories()
)
