source("http://bioconductor.org/biocLite.R")
biocLite("DNAcopy")

getoptUrl="http://lib.stat.cmu.edu/R/CRAN/bin/macosx/contrib/3.1/getopt_1.20.0.tgz"
if (!("getopt" %in% rownames(installed.packages))) {
  install.packages(getoptUrl, repos=NULL, type="source")
}
optparseUrl="http://cran.r-project.org/src/contrib/optparse_1.3.0.tar.gz"
if (!("optparse" %in% rownames(installed.packages))) {
  install.packages(optparseUrl, repos=NULL, type="source")
}
dependencies = c("naturalsort")
if (!all(dependencies %in% rownames(installed.packages()))) {
  install.packages(setdiff(dependencies, rownames(installed.packages())), repos="http://cran.cnr.Berkeley.edu")
}
q(save="no")
