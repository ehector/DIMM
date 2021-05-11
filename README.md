# About

This repository is an R package implementing methods from the following publication:
- Emily C. Hector and Peter X.-K. Song (2020). A distributed and integrated method of moments for high-dimensional correlated data analysis. 
Journal of the American Statistical Association, pages 1–14. doi: 10.1080/01621459.2020.1736082.

Briefly, the method performs regression analysis of high-dimensional correlated responses. It divides outcomes into blocks according to a supplied indicator, 
analyses blocks using composite likelihood, and combines blocks using a meta-estimator asymptotically equivalent to the optimal generalized method of moments (GMM) equation.

Please email ehector@ncsu.edu with any questions or bug-reports.

# Installation

The DIMM R package can be installed in one of two ways:
- from the downloaded gzipped tarball as R CMD INSTALL DIMM_1.0-1.tar.gz
- from the downloaded and renamed DIMM folder as R CMD build DIMM and R CMD INSTALL DIMM_1.0-1.tar.gz

Please make sure to have all packages listed in the DESCRIPTION file already installed.

# Citation

If you use the DIMM R package, please consider citing the relevant manuscript: Hector & Song (2020).

# References

The posdef.matrix function was written by Ravi Varadhan: https://stat.ethz.ch/pipermail/r-help/2008-February/153708.

Emily C. Hector and Peter X.-K. Song (2020). A distributed and integrated method of moments for high-dimensional correlated data analysis. Journal of the American Statistical Association, pages 1–14. doi: 10.1080/01621459.2020.1736082.

Lars P. Hansen. (1982), Large sample properties of generalized method of moments estimators, Econometrica 50(4), 1029-1054.

Christiano Varin, Nancy Reid & David Firth. (2011), An overview of composite likelihood methods, Statistica Sinica 21(1), 5-42.
