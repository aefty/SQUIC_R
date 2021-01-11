SQUIC_DEMO.load_data<-function(type="trid",p_power=5,n=100,normalized=FALSE)
{

	set.seed(1);

    start_time <- Sys.time()
    matrix_folder=system.file("extdata",package = "SQUIC")

    p=4^p_power;
	
    print(sprintf("# Reading Matrix From file: type=%s p=%d n=%d normalized=%d",type,p,n,normalized));

	iC_file_name=paste(matrix_folder,"/",type,p,"iC",".rmat",sep = "");
	C_file_name=paste(matrix_folder,"/",type,p,"C",".rmat",sep = "");

	iC_star=Matrix::readMM(iC_file_name)
	C_star=Matrix::readMM(C_file_name)
	
	# Generate data
    print(sprintf("# Generating data ...",type,p,n,normalized));
	mu_star <- replicate(p, 0);
 	data <- MASS::mvrnorm(n, mu_star, C_star, tol = 1e-2, empirical = FALSE, EISPACK = FALSE);
	data <- Matrix::t(data);

	if(normalized==TRUE){
		for (i in 1:p) {
			sd_data<-sd(data[i,]);
			data[i,]<-data[i,]/sd_data;
		}
	}
	finish_time <- Sys.time()
	print(sprintf("# Generating data finished: time=%f",finish_time-start_time));

	output <- list(
		"data" = data, 
		"X_star" = iC_star,
		"W_star" = C_star
	);

	return(output);
}



SQUIC_DEMO.cv <- function(type="trid",criterion="LL", p_power=5, n=100, K=4, lambda_sample=.4, max_iter=10,tol=1e-3,lambda_set_length=10,M=NULL) 
{

	time_naive	<- replicate(lambda_set_length, 0);
	f1_naive	<- replicate(lambda_set_length, 0);	
	time_cv = 0;
	f1_cv = 0;

	# Generate data
	out       <- SQUIC::SQUIC_DEMO.load_data( type=type , p_power=p_power , n=n ,normalized=TRUE);
	X_star    <- out$X_star;
	data_full <- out$data;

	#########################
	## Naive Search SQUIC ##
	#########################

	# This is the niave lambda range
	lambda_set <- seq(1 , .1 , -1/lambda_set_length);

	for (i in 1:lambda_set_length)
	{
		lambda=lambda_set[i];
		out           <- SQUIC::SQUIC_DEMO.compare(alg="SQUIC", data_full=data_full , lambda=lambda , tol=tol , max_iter=max_iter , X_star=X_star , M=M , verbose = 0);
		time_naive[i] <- out$time;
		f1_naive[i]   <- out$f1;
	}

	#########################
	## CV SQUIC Search ##
	#########################
	start_time <- Sys.time();

    # First get lambda_set
	print(sprintf("#Search for lambda_set"));
	out        <- SQUIC::SQUIC.lambda_set(data_full=data_full, lambda_sample=lambda_sample , lambda_set_length=lambda_set_length , M=M );
	lambda_set <- out$lambda_set;

	print(sprintf("#lambda_set sample results: S_abs_max%f  S_abs_min=%f  S_abs_mean=%f  S_abs_sd=%f",
		out$S_abs_max, 
		out$S_abs_min, 
		out$S_abs_mean, 
		out$S_abs_sd
		));
	print(lambda_set);	

	# Cross validation to find best lambda
	out  <- SQUIC::SQUIC.cv(data_full=data_full,lambda_set=lambda_set,K=K,tol=1e-3,max_iter=3,criterion=criterion ,M=M ,X0=NULL,W0=NULL);
	lambda_cv <-out$lambda;

	out <- SQUIC::SQUIC(data_train=data_full, lambda=lambda_cv, max_iter=max_iter, drop_tol=tol, term_tol=tol, M=M , X0=NULL, W0=NULL,data_test=NULL,verbose=0)

    start_end <- Sys.time();

	# Convert matrix to labels for accurecy test
	X_label      <-  as.vector( ((out$X)!=0)*1 );
	X_star_label <-  as.vector( ((X_star)!=0)*1 );

	output <- list(
		"time_naive"     =  time_naive, 
		"f1_naive" 	 	 =  f1_naive, 
		"time_cv" 	 	 =  start_end - start_time ,
		"f1_cv" 	 	 =  MLmetrics::F1_Score(X_star_label,X_label, positive = "1"),
		"lambda_cv"	 	 =	lambda_cv,
		"lambda_set"     =  lambda_set,
		"X" = out$X
	);

	return(output);
}


#done
SQUIC_DEMO.performance <- function(type="trid",lambda=0.4,n=100,tol=1e-4,max_iter=10) {

    #Hard coded values
    p_power_max<-6;

	time_squic		<-replicate(p_power_max, 0);
	time_equal		<-replicate(p_power_max, 0);	
	time_bigquic	<-replicate(p_power_max, 0);
	time_quic		<-replicate(p_power_max, 0);

	for (i in 1:p_power_max) {

        p_power=i

        # Generate data
	    out<-SQUIC::SQUIC_DEMO.load_data(type=type ,p_power=p_power ,n=n ,normalized=TRUE);
	    X_star<-out$X_star;
	    data_full<-out$data;

		print(sprintf("Benchmark for p=%d started",4^p_power));

		out<-SQUIC::SQUIC_DEMO.compare(alg="SQUIC" , data_full=data_full , lambda=lambda ,tol=tol , max_iter=max_iter ,  M=NULL,  X_star=NULL, erbose = 0);
		time_squic[i]<-out$time;

		out<-SQUIC::SQUIC_DEMO.compare(alg="EQUAL" , data_full=data_full ,  lambda=lambda ,tol=tol , max_iter=max_iter ,  M=NULL, X_star=NULL, verbose = 0);
		time_equal[i]<-out$time;

		out<-SQUIC::SQUIC_DEMO.compare(alg="BigQuic",data_full=data_full , lambda=lambda , tol=tol , max_iter=max_iter ,  M=NULL, X_star=NULL, verbose = 0);
		time_bigquic[i]<-out$time;

		out<-SQUIC::SQUIC_DEMO.compare(alg="QUIC" ,  data_full=data_full , lambda=lambda , tol=tol , max_iter=max_iter ,  M=NULL, X_star=NULL, verbose = 0);
		time_quic[i]<-out$time;			
	}

	output <- list(
		"time_squic"   			= time_squic, 
		"time_equal" 			= time_equal, 
		"time_bigquic" 			= time_bigquic, 		
		"time_quic" 			= time_quic 			
		)

	return(output);
}


#done
SQUIC_DEMO.drop_tol <- function(type="trid",p_power=5,n=100,lambda=0.3,term_tol=1e-12,max_iter=10,verbose=0) {

    #hard code
    drop_tol_max_power=10

	time		<-replicate(drop_tol_max_power, 0);
	f1			<-replicate(drop_tol_max_power, 0);
	acc			<-replicate(drop_tol_max_power, 0);
	nnz_W		<-replicate(drop_tol_max_power, 0);
	nnz_X		<-replicate(drop_tol_max_power, 0);	
	norm_f_XW	<-replicate(drop_tol_max_power, 0);	

	# Generate data
	out<-SQUIC::SQUIC_DEMO.load_data(type=type ,p_power=p_power ,n=n ,normalized=TRUE);
	iC_star<-out$iC_star;
	data_full<-out$data;

	for (i in 1:drop_tol_max_power) {

		drop_tol=10^(-1*i);

		print(sprintf("Test for drop_tol=%f started",drop_tol));
		time_start <- Sys.time();
		out	 	   <-SQUIC::SQUIC(data_train=data_full,lambda=lambda, max_iter=max_iter, drop_tol=drop_tol, term_tol=term_tol, M=NULL, X0=NULL, W0=NULL, data_test=NULL,verbose=verbose);
		end_start  <- Sys.time();

		X				<- out$X;
		W		 		<- out$W;
		nnz_W[i] 		<- Matrix::nnzero(W);
		nnz_X[i] 		<- Matrix::nnzero(X);
		norm_f_XW[i]	<- Matrix::norm( X%*%W - Matrix::Diagonal(nrow(X)) ,"f");

		# Convert matrix to labels
		X_label 		<-  as.vector( ((X)!=0)*1 );
		X_star_label	<-  as.vector( ((X_star)!=0)*1 );

		time[i]	<- end_start-time_start;
		f1[i]	<- MLmetrics::F1_Score(X_star_label,X_label, positive = "1");
		acc[i]	<- MLmetrics::Accuracy(X_star_label,X_label);	
	}

	output <- list(
		"time"		= time, 
		"f1"		= f1, 
		"acc"		= acc, 		
		"nnz_X"		= nnz_W,	
		"nnz_X"		= nnz_X,	
		"norm_f_XW"	= norm_f_XW
		)

	return(output);
}


SQUIC_DEMO.compare <- function(alg,data_full,lambda=0.5,tol=1e-4,max_iter=10, X_star= NULL, M=NULL, verbose = 0) {
	

		data_full_t <- Matrix::t(data_full);

		time_start <- Sys.time()
		if(alg=="QUIC")
		{
		print("#QUIC")
		# QUIC
		out				<-QUIC::QUIC(S=cov(data_full_t), rho=lambda, path = NULL, tol = tol, msg = verbose, maxIter = max_iter, X.init =NULL, W.init = NULL)
		X			    <-out$X;
		}
		else if(alg=="BigQuic")
		{
		print("#BigQuic")
		# BigQuic
		out				<-BigQuic::BigQuic(X = data_full_t, inputFileName = NULL, outputFileName = NULL, lambda = lambda, numthreads = 8, maxit = max_iter, epsilon = tol, k = 0, memory_size = 8000, verbose = verbose, isnormalized = 1, seed = NULL, use_ram = TRUE);
		X      			<-out$precision_matrices[[1]];
		}
		else if(alg=="SQUIC")
		{
		print("#SQUIC")
		# SQUIC
		out	 		 <-SQUIC::SQUIC(data_train=data_full,lambda=lambda, max_iter=max_iter, drop_tol=tol, term_tol=tol, M=M , X0=NULL, W0=NULL, data_test=NULL,verbose=verbose);
		X      		 <-out$X;
		}
		else if(alg=="EQUAL")
		{
		print("#EQUAL")
		# EQUAL
		out				 <-EQUAL::EQUAL(X=data_full_t,type=TRUE,sdiag=FALSE,lambda=lambda,lambda.min=sqrt(log(ncol(data_full_t))/nrow(data_full_t)),nlambda=1,err=tol,maxIter =max_iter,rho=1)
		X		         <-out$Omega[[1]];
		}
		else
		{
			stop("Algorithem not found");
		};
		time_end 	<- Sys.time()

	if(!is.null(X_star))
	{
		# Convert matrix to labels
        print("#Computing F1 Score & Accuracy")
		X_label <-  as.vector( ((X)!=0)*1 );
		X_star_label <-  as.vector( ((X_star)!=0)*1 );

		output <- list(
			"time" = time_end-time_start,
			"X" = X, 
			"f1" = MLmetrics::F1_Score(X_star_label,X_label, positive = "1"),
			"acc" = MLmetrics::Accuracy(X_star_label,X_label)				
		);
	}else{
		# Convert matrix to labels
		output <- list(
			"X" = X, 
			"time" = time_end-time_start	
		);

	}

	return(output);
}


