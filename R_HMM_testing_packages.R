###############################################################################
# If you edit this file you MUST release a new version of the gatkbase docker #
# built with the updated r dependencies                                       #
#                                                                             #
# you MUST also manually clear the travis cache for master before running the #
# pull request tests in order to make sure it's still working                 #
###############################################################################

options(warn = 2)     # treat warnings as errors, otherwise script can fail silently if a package fails to install

InstallPackageFromArchive = function(packageName, packageURL) {
    # make sure to use http not https as this will give an "unsupported URL scheme" error
    if (!(packageName %in% rownames(installed.packages()))) {
        install.packages(packageURL, repos = NULL, type = "source", clean = TRUE)
    }
}

# HMM is only required for testing and not in production:
InstallPackageFromArchive("HMM", "http://cran.r-project.org/src/contrib/HMM_1.0.tar.gz")

q(save = "no")
