source("http://bioconductor.org/biocLite.R")
biocLite("DNAcopy")
#Make sure to use http not https as this will give an "unsupported URL scheme" error
getoptUrl="http://cran.r-project.org/src/contrib/getopt_1.20.0.tar.gz"
if (!("getopt" %in% rownames(installed.packages()))) {
  install.packages(getoptUrl, repos=NULL, type="source")
}
optparseUrl="http://cran.r-project.org/src/contrib/Archive/optparse/optparse_1.3.2.tar.gz"
if (!("optparse" %in% rownames(installed.packages()))) {
  install.packages(optparseUrl, repos=NULL, type="source")
}
datatableUrl="http://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.10.4-2.tar.gz"
if (!("data.table" %in% rownames(installed.packages()))) {
    install.packages(datatableUrl, repos=NULL, type="source")
}
dependencies = c("naturalsort","ggplot2","gplots","reshape","gsalib")
repos <- c("http://cran.cnr.Berkeley.edu",
           "https://cran.mtu.edu",
           "http://lib.stat.cmu.edu/R/CRAN/")
missing <- which(!(dependencies %in% rownames(installed.packages())))
try <- 1
while(length(missing)!=0 & try <= length(repos)) {
    install.packages(dependencies[missing], repos = repos[try])
    missing <- which(!(dependencies %in% rownames(installed.packages())))
    try <- try + 1
}

# HMM is only required for testing and not in production:
hmmUrl = "http://cran.r-project.org/src/contrib/HMM_1.0.tar.gz"
if (!("HMM" %in% rownames(installed.packages()))) {
  install.packages(hmmUrl, repos=NULL, type="source")
}
q(save="no")
