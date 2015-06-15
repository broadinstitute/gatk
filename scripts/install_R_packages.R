source("http://bioconductor.org/biocLite.R")
biocLite("DNAcopy")
dependencies = c("argparser", "naturalsort")
if (length(setdiff(dependencies, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(dependencies, rownames(installed.packages())), repos="http://cran.cnr.Berkeley.edu")
}
q(save="no")
