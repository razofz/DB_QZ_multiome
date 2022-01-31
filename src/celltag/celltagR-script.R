install_non_conda_package <- function(package) {
        if (!requireNamespace(package, quietly=T)) {
	    library(devtools)
	    # BiocManager::install(package)
	    devtools::install(package)
    }
}
install_non_conda_package('GenomeInfoDbData')
install_non_conda_package('EnsDb.Hsapiens.v86')

