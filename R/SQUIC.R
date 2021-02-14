usethis::use_package("Matrix") 
usethis::use_package("MLmetrics") # Use for multivriate data generation

# For validation 
usethis::use_package("BigQuic") 

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
  
  p	= nrow(Y1);


  if(lambda<0){
	  stop('#SQUIC: lambda cannot be negative.');
  }
  if(max_iter<0){
	  stop('#SQUIC: max_iter cannot be negative.');
  }
  
  if(drop_tol<0){
	  stop('#SQUIC: drop_tol cannot be negative.');
  }
  
  if(term_tol<0){
	  stop('#SQUIC: term_tol cannot be negative.');
  } 

  if(mode < 0 || mode > 9){
	  stop('#SQUIC: mode must be in [0,9].');
  } 


  if(is.null(M)){
	  # Make empty sparse matrix of type dgCMatrix.
	  M	= as(Matrix::sparseMatrix(dims = c(p,p), i={}, j={}),"dgCMatrix");
  }else{
    if(!isSymmetric(M)){
      stop('#SQUIC: M must be symmetric.');
    }
  }

  if(is.null(X0) || is.null(W0)){
	  # Make empty sparse matrix of type dgCMatrix.
	  X0 = as(Matrix::sparseMatrix(dims = c(p,p), i={}, j={}),"dgCMatrix");
	  W0 = as(Matrix::sparseMatrix(dims = c(p,p), i={}, j={}),"dgCMatrix");
  }else{
     if(!isSymmetric(X0) || isSymmetric(W0)){
      stop('#SQUIC: X0 and W0 must be symmetric.')
    }
  }
  
  if(is.null(Y2)){
	  # Make 1x1 matrix for Y2. This input must be provided , 
    # if p_test!=p_train ... thatn the test data is ignored.
	  Y2 = matrix(1.0, nrow = 1, ncol = 1)
  }else{
     if(nrow(Y2) != p ){
      stop('#SQUIC: Y2 must the same number of rows as Y1.')
    }
  }

   out = SQUIC::SQUIC_R(Y1 , lambda , max_iter , drop_tol , term_tol , verbose , mode ,M , X0 , W0 , Y2 );
   return(out);
}

# Sparse Sample covariance matrix
SQUIC_S<-function(data, lambda_sample=.2, M=NULL){

	# Get sample covarinace matrix by running SQUIC with max_iter=0;
	squic_output = SQUIC::SQUIC(Y1=data,lambda=lambda_sample, max_iter=0, drop_tol=0, term_tol=0,verbose=0, M=M );

	output	= list(
		"S"             = squic_output$S					
	);

	return(output);
}

# Cross validation
SQUIC_CV<-function(data , lambda_set , K=5 , drop_tol=1e-3 , term_tol=1e-2 , max_iter=3  , M=NULL , X0=NULL , W0=NULL)
{
	#Dimensions
	p		= nrow(data);
	n_full	= ncol(data);

	#Construct active sample sets
	lambda_set	= sort(lambda_set,decreasing =TRUE)
	nlambda		= length(lambda_set)
	active_set	= split(1:n_full, rep(1:K, length = n_full));

	# Cross Validation Matrix
	CV_AIC	= matrix(0,nrow=K,ncol=nlambda);
	CV_BIC	= matrix(0,nrow=K,ncol=nlambda);

	# For each active set of sample indicies ...
	for (k in 1:K) 
	{
		# Construct the test and train data
		data_train	= data[,-active_set[[k]]];
		data_test	= data[,active_set[[k]]];
		n_train		= ncol(data_train);
		n_test		= ncol(data_test);

		# for each lambda in the lambda set
		for (l in 1:nlambda) 
		{ 
			# Set of lambda 
			lambda	=lambda_set[l];
	
			# run a rough (low tolerences and iterations)
			out	= SQUIC::SQUIC( Y1=data_train , lambda=lambda , max_iter=max_iter , drop_tol=drop_tol , term_tol=term_tol , verbose=0 , M=M , X0=X0 , W0=W0 , Y2=data_test );

			#Extract the results form SQUIC
			X			= out$X;
			logdetX_Y1	= out$info_logdetX_Y1;
			trXS_Y2		= out$info_trXS_Y2;
			nnzX		= Matrix::nnzero(X);

			#logliklihood
			logliklihood = (-p*log(2*3.14)  +  logdetX_Y1 - trXS_Y2 )*n_test/2;

			# The inverse covariance is symmetric matrix
			num_params = (p+(nnzX-p)/2);

			CV_AIC[k,l]	 = (2*num_params - 2*logliklihood);
			CV_BIC[k,l]	 = (log(n_train)*num_params - 2*logliklihood) ;
			
		}
	}

	# Find the smallest value in the CV matrix and select the corresponding lambda
	CV_AIC_mean	 = colMeans(CV_AIC);
	CV_BIC_mean	 = colMeans(CV_BIC);

	lambda_opt_AIC	= lambda_set[which.min(CV_AIC_mean)];
	lambda_opt_BIC	= lambda_set[which.min(CV_BIC_mean)];

	output	= list(
		"CV_AIC_mean"		= CV_AIC_mean,
		"CV_BIC_mean"		= CV_BIC_mean,	
		"lambda_opt_AIC"	= lambda_opt_AIC,		
		"lambda_opt_BIC"	= lambda_opt_BIC				
	);

	return(output);
}

# Cross validation auto search
SQUIC_CVX<-function(data , lambda=1 , alpha=1, nnzX_per_row_max=100 , R=20 , K=5 , criterion="AIC" , drop_tol=1e-3 , term_tol=1e-2 , max_iter=3 , M=NULL , X0=NULL , W0=NULL)
{
	if(nnzX_per_row_max<=0){
		stop("#SQUIC: nnzX_per_row_max nonzero postive number.")
	}

	#Dimensions
	p		= nrow(data);
	n_full	= ncol(data);

	#Construct active sample sets
	active_set	= split(1:n_full, rep(1:K, length = n_full));

	# Cross Validation Matrix
	CV	= vector(length=K);


	lambda_set = c();
	f0_set = c();
	f1_set = c();
	f2_set = c();
	alpah_set=c();

	for (r in 1:R) # Line Search iteration
	{

		nnzX_per_row_sum=0;

		for (k in 1:K) # For each active set of sample indicies ...
		{
			# Construct the test and train data
			data_train	= data[,-active_set[[k]]];
			data_test	= data[,active_set[[k]]];
			n_train		= ncol(data_train);
			n_test		= ncol(data_test);

			# run a rough (low tolerences and iterations)
			out	= SQUIC::SQUIC( Y1=data_train , lambda=lambda , max_iter=max_iter , drop_tol=drop_tol , term_tol=term_tol , verbose=0 , M=M , X0=X0 , W0=W0 , Y2=data_test );

			#Extract the results form SQUIC
			X				 = out$X;
			logdetX_Y1		 = out$info_logdetX_Y1;
			trXS_Y2			 = out$info_trXS_Y2;
			nnzX			 = Matrix::nnzero(X);
			nnzX_per_row_sum = nnzX/p + nnzX_per_row_sum;

			#logliklihood
			logliklihood = (-p*log(2*3.14)  +  logdetX_Y1 - trXS_Y2 )*n_test/2;

			# The inverse covariance is symmetric matrix
			num_params = (p+(nnzX-p)/2);

			if(criterion=="AIC")
			{
				CV[k] = (2*num_params - 2*logliklihood);
			}else if(criterion=="BIC")
			{
				CV[k] = (log(n_train)*num_params - 2*logliklihood) ;
			}else{	
				stop("Unkown Cirterion");
			}
		}

		# Average Number nonzers in X
		nnzX_per_row_avg = nnzX_per_row_sum/K;

		# Steph size
		alpha = max(nnzX_per_row_max - nnzX_per_row_avg,0)/nnzX_per_row_max;
		alpah_set<-c(alpah_set,alpha);

		# Function Values 
		f0 =	mean(CV);
		f0_set <- c( f0_set , f0 );
		
		# Delta x
		h=lambda_set[r] - lambda_set[r-1];
		
		# Delta f (Gradient)
		if(r>1)
		{
			f1 = (f0_set[r] - f0_set[r-1]) / h;
		}else{
			f1 = 0;
		}
		f1_set <- c( f1_set , f1 );

		# Delta^2 f (Hessian)
		if(r>2)
		{
			f2 = (f1_set[r] - f1_set[r-1]) / h;
		}else{
			f2 = 0;
		}
		f2_set <- c( f2_set , f2 );


		#Store Lambda
		lambda_set <- c(lambda_set, lambda);


		if(r<3){

			lambda_new =  lambda *0.95;

		}else{

			lambda_new =  lambda - alpha * f1/f2;

		}



		if(lambda_new - lambda < 1e-2){

			break;

		}








		print(sprintf("#Iteration %d Done - lambda=%f %s=%f g=%f", r, lambda, criterion, CV_mean[r],g));

		# Break condition
	
		if( g < 0 || g==0 )
		{
			print(sprintf("#Optimal lambda=%f g=%f", lambda,g));
			#break;
		}

		update= -g*alpha;

	
		lambda = lambda-alpha;

		if(lambda<0){
			break;
		}

		print(sprintf("g=%f alpha=%f update=%f next_lambda=%f",g,alpha,update,lambda));
	
	}

	# Find the smallest value in the CV matrix and select the corresponding lambda
	lambda_opt	= lambda_set[which.min(CV_mean)];

	output	= list(
		"lambda_set"	= lambda_set,
		"CV_mean"		= CV_mean,
		"lambda_opt"	= lambda_opt			
	);

	return(output);
}

# Cross validation auto search
SQUIC_CVX2<-function(data , lambda=1 , nnzX_per_row_max=100 , R=20 , K=5 , criterion="AIC" , drop_tol=1e-3 , term_tol=1e-2 , max_iter=3 , M=NULL , X0=NULL , W0=NULL)
{

	if(nnzX_per_row_max<=0){
		stop("#SQUIC: nnzX_per_row_max nonzero postive number.")
	}

	#Dimensions
	p		= nrow(data);
	n_full	= ncol(data);

	#Construct active sample sets
	active_set	= split(1:n_full, rep(1:K, length = n_full));

	# Cross Validation Matrix
	CV	= vector(length=K);
	CV_mean = c();
	lambda_set = c();
	g_set = c();
	alpah_set=c();



	g=1;
	alpha=1;
	update = 1/2;

	lambda=lambda*2;

	

	for (r in 1:R) # Line Searh iteration
	{
		for (k in 1:K) # For each active set of sample indicies ...
		{
			# Construct the test and train data
			data_train	= data[,-active_set[[k]]];
			data_test	= data[,active_set[[k]]];
			n_train		= ncol(data_train);
			n_test		= ncol(data_test);

			# run a rough (low tolerences and iterations)
			out	= SQUIC::SQUIC( Y1=data_train , lambda=lambda , max_iter=max_iter , drop_tol=drop_tol , term_tol=term_tol , verbose=0 , M=M , X0=X0 , W0=W0 , Y2=data_test );

			#Extract the results form SQUIC
			X				= out$X;
			logdetX_Y1		= out$info_logdetX_Y1;
			trXS_Y2			= out$info_trXS_Y2;
			nnzX			= Matrix::nnzero(X);
			nnzX_per_row	= nnzX/p;

			#logliklihood
			logliklihood = (-p*log(2*3.14)  +  logdetX_Y1 - trXS_Y2 )*n_test/2;

			# The inverse covariance is symmetric matrix
			num_params = (p+(nnzX-p)/2);

			if(criterion=="AIC")
			{
				CV[k] = (2*num_params - 2*logliklihood);
			}else if(criterion=="BIC")
			{
				CV[k] = (log(n_train)*num_params - 2*logliklihood) ;
			}else{	
				stop("Unkown Cirterion");
			}
		
		}

		CV_mean <- c(CV_mean, mean(CV));
		lambda_set <- c(lambda_set, lambda);


		print(sprintf("#Iteration %d Done - lambda=%f %s=%f", r, lambda, criterion, CV_mean[r]));

		# Break condition
		if(r>1){
			g = (CV_mean[r]- CV_mean[r-1]) / (lambda_set[r]- lambda_set[r-1]) ;
			g = g /  (CV_mean[1]- CV_mean[2]) / (lambda_set[1]- lambda_set[2]) ;

			alpha = max(nnzX_per_row_max - nnzX_per_row,0)/nnzX_per_row_max;
			update= g*alpha;



			print(sprintf("g=%f alpha=%f update=%f",g,alpha,update));

			if( update < 0 || update==0 )
			{
				print(sprintf("#Optimal lambda=%f", lambda));
				break;
			}

		}

		g_set <- c(g_set, g);
		alpah_set <- c(alpah_set, alpha);
		

		lambda = lambda*update;

		if(lambda<0){
			break;
		}
		
		print(sprintf("g_set=%f",g_set));
		print(sprintf("alpah_set=%f",alpah_set));
		print(sprintf("lambda_set=%f",lambda_set));
		print(sprintf("CV_mean=%f",CV_mean));

		
		
	}

	# Find the smallest value in the CV matrix and select the corresponding lambda
	lambda_opt	= lambda_set[which.min(CV_mean)];

	output	= list(
		"lambda_set"	= lambda_set,
		"CV_mean"		= CV_mean,
		"lambda_opt"	= lambda_opt			
	);

	return(output);
}