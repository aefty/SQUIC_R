# SQUIC_R

// The following system environment variables need to be set:

export KMP_DUPLICATE_LIB_OK=TRUE
export R_LD_LIBRARY_PATH=$HOME

// Run the following command to install the library:

library(devtools)
install_github("aefty/SQUIC_Pkg_R")

Note to install QUIC run the following comand in R:
install.packages("https://cran.r-project.org/src/contrib/Archive/QUIC/QUIC_1.1.1.tar.gz", type="source")

library(SQUIC)

p=10
n=130
lambda=.5
max_iter=10
tol=1e-3

iC_star <- Matrix::Diagonal(p);
Data_train<-SQUIC_generate_data(iC_star,n,normalized=TRUE);
Data_test<-SQUIC_generate_data(iC_star,3,normalized=TRUE);

a<-SQUIC(Data_train=Data_train,lambda=lambda, max_iter=max_iter, drop_tol=tol, term_tol=tol, M=NULL, X0=NULL, W0=NULL, Data_test=Data_test)

a<-SQUIC(Data_train=Data_train,lambda=lambda, max_iter=0, drop_tol=tol, term_tol=tol)

S_test=cov(t(as.matrix(Data_test))) \*(ncol(Data_test)-1)/ncol(Data_test);

sum(diag(a$X %\*% S_test))

compileAttributes(verbose=TRUE)
