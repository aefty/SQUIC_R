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
SQUIC <- function(Y1, lambda, max_iter, drop_tol, term_tol,verbose=1, mode=0, M=NULL, X0=NULL, W0=NULL, Y2=NULL) {
  
  p<-nrow(Y1);

  if(is.null(M)){
	  # Make empty sparse matrix of type dgCMatrix.
	  M<-as(Matrix::sparseMatrix(dims = c(p,p), i={}, j={}),"dgCMatrix");
  }else{
    if(!isSymmetric(M)){
      stop('M must be symmetric.');
    }
  }

  if(is.null(X0) || is.null(W0)){
	  # Make empty sparse matrix of type dgCMatrix.
	  X0<-as(Matrix::sparseMatrix(dims = c(p,p), i={}, j={}),"dgCMatrix");
	  W0<-as(Matrix::sparseMatrix(dims = c(p,p), i={}, j={}),"dgCMatrix");
  }else{
     if(!isSymmetric(X0) || isSymmetric(W0)){
      stop('X0 and W0 must be symmetric.')
    }
  }
  
  if(is.null(Y2)){
	  # Make 1x1 matrix for Y2. This input must be provided , 
    # if p_test!=p_train ... thatn the test data is ignored.
	  Y2<-matrix(1.0, nrow = 1, ncol = 1)
  }else{
     if(nrow(Y2) !=p ){
      stop('Y2 must the same number of rows as Y1.')
    }
  }

   out<-SQUIC::SQUIC_R(Y1 , lambda , max_iter , drop_tol , term_tol , verbose , mode ,M , X0 , W0 , Y2 );
   return(out);
}

SQUIC_S<-function(data, lambda_sample=.5,lambda_set_length=10 , M=NULL){

	p=nrow(data);	
	
	# Get sample covarinace matrix by running SQUIC with max_iter=0;
	squic_output<-SQUIC::SQUIC(Y1=data,lambda=lambda_sample, max_iter=0, drop_tol=0, term_tol=0,verbose=0, M=M );

	# Get absolute value of the max and mean of nonzeros in S
	S<-squic_output$S;
	S_abs_vec = abs(S@x);

	S_abs_max     <- max(S_abs_vec);
	S_abs_min     <- min(S_abs_vec);
	S_abs_mean    <- mean(S_abs_vec);
	S_abs_sd      <- sd(S_abs_vec);
	S_nnz_per_row <- length(S_abs_vec)/p;

	up  <- (lambda_sample + S_abs_max)/2;
	low <- (lambda_sample + S_abs_min)/2;

	delta      <- (low-up)/(lambda_set_length-1);
	lambda_set <- seq(up , low , delta);
	
	output <- list(
		"S_nnz_per_row" = S_nnz_per_row, 
		"S_abs_max"		= S_abs_max, 
		"S_abs_min" 	= S_abs_min, 
		"S_abs_mean"	= S_abs_mean, 
		"S_abs_sd" 		= S_abs_sd,
		"lambda_set"	= lambda_set 					
	);

	return(output);
}



# Cross validation
SQUIC_CV<-function(data , lambda_set, K=5 , criterion="AIC", drop_tol=1e-3,term_tol=1e-2 , max_iter=5  , M=NULL , X0=NULL , W0=NULL)
{

	p=nrow(data);
	n_full=ncol(data);

	#Construct active sample sets
	lambda_set<-sort(lambda_set,decreasing =TRUE)
	nlambda <- length(lambda_set)
	active_set<-split(1:n_full, rep(1:K, length = n_full));

	# Cross Validation Matrix
	CV<-matrix(0,nrow=K,ncol=nlambda);

	U_set=list();
	for (k in 1:K) # For each active set of sample indicies ...
	{
		# Construct the test and train data
		data_train<-data[,-active_set[[k]]];
		data_test<-data[,active_set[[k]]];
		n_train<-ncol(data_train);
		n_test<-ncol(data_test);

		temp_list=list();
		for (l in 1:nlambda){ # for each lambda in the lambda set

			# Set of lambda 
			lambda=lambda_set[l];
	
			# run a rough (low tolerences and iterations)
			out<-SQUIC::SQUIC( Y1=data_train , lambda=lambda , max_iter=max_iter , drop_tol=drop_tol , term_tol=term_tol , verbose=0 , M=M , X0=X0 , W0=W0 , Y2=data_test );

			#Extract the results form SQUIC
			X<-out$X;
			logdetX<-out$info_logdetX_Y1;
			trXS_test<-out$info_trXS_Y2;
			nnzX<-Matrix::nnzero(X);

			#logliklihood
			logliklihood<- (-p*log(2*3.14)  +  logdetX - trXS_test )*n_test/2;

			# The inverse covariance is symmetric matrix
			num_params = (p+(nnzX-p)/2);

			if(criterion=="AIC")
			{
				CV[k,l] <-  (2*num_params - 2*logliklihood) ;
			}else if(criterion=="AICc")
			{
				CV[k,l] <-  (2*num_params - 2*logliklihood) + 2*(num_params^2+num_params)/(n_train-num_params-1);

			}else if(criterion=="BIC")
			{
				CV[k,l] <-  (log(n_train)*num_params - 2*logliklihood) ;
			}else{
				stop("Unkown criterion")
			}
		}
	}

	# Find the smallest value in the CV matrix and select the corresponding lambda
	CV_mean    <-colMeans(CV);
	l_inx	   <-which.min(CV_mean);
	lambda_opt <-lambda_set[l_inx];

	output <- list(
		"lambda_opt"   = lambda_opt,
    	"CV_mean"	   = CV_mean
	);

	return(output);
}