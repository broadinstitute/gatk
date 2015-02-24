dependencies = c("ggplot2","gplots","reshape","gsalib")
if (length(setdiff(dependencies, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(dependencies, rownames(installed.packages())), repos="http://cran.cnr.Berkeley.edu")  
}
q(save="no")
