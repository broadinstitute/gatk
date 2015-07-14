source("http://bioconductor.org/biocLite.R")
biocLite("DNAcopy")
argparserurl="http://cran.r-project.org/src/contrib/Archive/argparser/argparser_0.1.tar.gz"
if (length(setdiff("argparser", rownames(installed.packages()))) > 0) {
  install.packages(argparserurl, repos=NULL, type="source")
}
dependencies = c("naturalsort")
if (length(setdiff(dependencies, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(dependencies, rownames(installed.packages())), repos="http://cran.cnr.Berkeley.edu")
}
q(save="no")
