usethis::use_package("Matrix") 
usethis::use_package("MASS") # Use for multivriate data generation
usethis::use_package("MLmetrics") # Use for multivriate data generation

# For validation 
usethis::use_package("BigQuic") # Use for multivriate data generation

# Here we set the envioerment variable KMP_DUPLICATE_LIB_OK=TRUE
# This is needed for Mac due to potential conflict with other OMP versions
.onLoad <- function(libname, pkgname){
	print("## SQUIC init start ##");
	print("1) Setting envirnment variable KMP_DUPLICATE_LIB_OK=TRUE")
	Sys.setenv(KMP_DUPLICATE_LIB_OK =TRUE)
	print("## SQUIC init finished ##");
}

# Main function
SQUIC <- function(data_train, lambda, max_iter, drop_tol, term_tol,verbose=1, mode=0, M=NULL, X0=NULL, W0=NULL, data_test=NULL) {
  
  p<-nrow(data_train);

  if(is.null(M)){
	  # Make empty sparse matrix.
	  M<-Matrix::sparseMatrix(dims = c(p,p), i={}, j={})
  }else{
    if(!isSymmetric(M)){
      stop('M must be symmetric.')
    }
  }

  if(is.null(X0) || is.null(W0)){
	  # Make empty sparse matrix.
	  X0<-Matrix::sparseMatrix(dims = c(p,p), i={}, j={})
	  X0<-Matrix::sparseMatrix(dims = c(p,p), i={}, j={})
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

# Cross validation
SQUIC_CV<-function(data_full , lambda_set,K=4, drop_tol=1e-2,term_tol=1e-3 , max_iter=3 , criterion="LL" , M=NULL , X0=NULL , W0=NULL)
{

	p=nrow(data_full);
	n_full=ncol(data_full);

	#Construct active sample sets
	lambda_set<-sort(lambda_set,decreasing =TRUE)
	nlambda=length(lambda_set)
	active_set<-split(1:n_full, rep(1:K, length = n_full));

	# Cross Validation Matrix
	CV<-matrix(0,nrow=K,ncol=nlambda);

	U_set=list();
	for (k in 1:K) # For each active set of sample indicies ...
	{
		# Construct the test and train data
		data_train<-data_full[,-active_set[[k]]];
		data_test<-data_full[,active_set[[k]]];
		n_train<-ncol(data_train);
		n_test<-ncol(data_test);

		temp_list=list();
		for (l in 1:nlambda){ # for each lambda in the lambda set

			# Set of lambda 
			lambda=lambda_set[l];
	
			# run a rough (low tolerences and iterations)
			squic_output<-SQUIC::SQUIC(data_train=data_train, lambda=lambda, max_iter=max_iter, drop_tol=drop_tol , term_tol=term_tol , verbose=0 , M=M , X0=X0, W0=W0 , data_test=data_test );

			#Extract the results form SQUIC
			X<-squic_output$X;
			logdetX<-squic_output$info_logdetX;
			trXS_test<-squic_output$info_trXS_test;
			nnzX<-Matrix::nnzero(X);

			#logliklihood (Not negative logliklihood!!!)
			logliklihood<-( p*log(2*3.14) + logdetX - trXS_test )*n_test/2;

			# For all criterion smaller is better
			if (criterion == "LL") # Defualt criterion is logliklihood (negative)
			{
				CV[k,l]=-logliklihood;
           	} 
		   	else if (criterion == "AIC") # AIC Criterion 
		   	{
                CV[k,l] = 2* (p+(nnzX-p)/2) - 2*logliklihood ;
            }
            else if (criterion == "BIC") # BIC Criterion
			{
                CV[k,l] = (p+(nnzX-p)/2) *log(n_test) - 2*logliklihood; 
            }
			else{
				stop("Criterion specfied is not valid.")
			}
		}
	}

	# Find the smallest value in the CV matrix and select the corresponding lambda
	CV_mean<-colMeans(CV);
	l_inx<-which.min(CV_mean);
	lambda_opt<-lambda_set[l_inx];

	output <- list(
		"lambda_opt"    = lambda_opt,
    "CV_mean"	      = CV_mean
	);

	return(output);
}