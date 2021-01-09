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

SQUIC <- function(data_train, lambda, max_iter, drop_tol, term_tol, M=NULL, X0=NULL, W0=NULL, data_test=NULL,verbose=1) {
  
	mode = c(1,1,1)

   # M is provided <=> mode[0] > 0;
   # Both X0 and W0 ae provided <=> mode[1] > 0;
   # Test data provided <=> mode[2] > 0;

  if(is.null(M)){
	  # Make dummy sparse matrix. Dimension dosent matter.
	  M<-Matrix::Diagonal(1);
      mode[1]=0;
    }

  if(is.null(X0) || is.null(W0)){
	  # Make dummy sparse matrix. Dimension dosent matter.
	  X0<-Matrix::Diagonal(1);
	  W0<-Matrix::Diagonal(1);
	  mode[2]=0;
  }
  
  if(is.null(data_test)){
	# Make dummy dense matrix. Dimension dosent matter.
	data_test<-matrix(1:4, nrow = 2, ncol = 2)
    mode[3]=0;
  }

   out<-.Call(`_SQUIC_SQUIC_BASE`, mode , data_train , lambda , max_iter , drop_tol , term_tol , M , X0 , W0 , data_test , verbose);
   
   if(max_iter>0)
   {
	   out$X = Matrix::drop0(out$X);
	   out$W = Matrix::drop0(out$W);  
   }

   return(out);
}