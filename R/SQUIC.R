usethis::use_package("Matrix") 
usethis::use_package("MASS") # Use for multivriate data generation
usethis::use_package("MLmetrics") # Use for multivriate data generation

# For validation 
usethis::use_package("BigQuic") # Use for multivriate data generation


.onLoad <- function(libname, pkgname){
	print("## SQUIC init start ##");
	print("1) Setting envirnment variable KMP_DUPLICATE_LIB_OK=TRUE")
	Sys.setenv(KMP_DUPLICATE_LIB_OK =TRUE)
	print("## SQUIC init finished ##");
}

SQUIC <- function(data_train, lambda, max_iter, drop_tol, term_tol,verbose=1, mode=0, M=NULL, X0=NULL, W0=NULL, data_test=NULL) {
  
	mode = c(1,1,1);
  p<-nrow(data_train);

   # M is provided <=> mode[0] > 0;
   # Both X0 and W0 ae provided <=> mode[1] > 0;
   # Test data provided <=> mode[2] > 0;

  if(is.null(M)){
	  # Make empty sparse matrix.
	  M<-Matrix::sparseMatrix(dims = dim(p,p), i={}, j={})
  }else{
    if(!isSymmetric(M)){
      stop('M must be symmetric.')
    }
  }

  if(is.null(X0) || is.null(W0)){
	  # Make empty sparse matrix.
	  X0<-Matrix::sparseMatrix(dims = dim(p,p), i={}, j={})
	  X0<-Matrix::sparseMatrix(dims = dim(p,p), i={}, j={})
  }else{
     if(!isSymmetric(X0) || isSymmetric(W0)){
      stop('X0 and W0 must be symmetric.')
    }
  }
  
  if(is.null(data_test)){
	  # Make 1x1 matrix for data_test.
	  data_test<-matrix(1, nrow = 1, ncol = 1)
  }else{
     if(nrow(data_test) !=p ){
      stop('data_test must the same number of rows as data_train.')
    }
  }

   out<-.Call(`_SQUIC_SQUIC_C`, data_train , lambda , max_iter , drop_tol , term_tol , verbose , mode ,M , X0 , W0 , data_test );
   
   return(out);
}