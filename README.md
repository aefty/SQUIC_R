# SQUIC_R

// The following system environment variables need to be set:

export KMP_DUPLICATE_LIB_OK=TRUE
export R_LD_LIBRARY_PATH=$HOME

// Run the following command to install the library:

library(devtools)
install_github("aefty/SQUIC_Pkg_R")

Note to install QUIC run the following comand in R:
install.packages("https://cran.r-project.org/src/contrib/Archive/QUIC/QUIC_1.1.1.tar.gz", type="source")
