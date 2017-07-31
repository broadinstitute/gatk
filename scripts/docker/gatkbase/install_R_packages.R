require(devtools)
source("http://bioconductor.org/biocLite.R")
biocLite("DNAcopy")

if (!("getopt" %in% rownames(installed.packages()))) {
install_version("getopt", version = "1.20.0", repos = "http://cran.us.r-project.org")
}

if (!("optparse" %in% rownames(installed.packages()))) {
install_version("optparse", version = "1.3.2", repos = "http://cran.us.r-project.org")
}

dependencies = c("naturalsort","ggplot2","gplots","reshape","gsalib")
if (!all(dependencies %in% rownames(installed.packages()))) {
  install.packages(setdiff(dependencies, rownames(installed.packages())), repos="http://cran.cnr.Berkeley.edu")
}

# HMM is only required for testing and not in production:
hmmUrl = "http://cran.r-project.org/src/contrib/HMM_1.0.tar.gz"
if (!("HMM" %in% rownames(installed.packages()))) {
install_version("HMM", version = "1.0", repos = "http://cran.us.r-project.org")
}
q(save="no")
