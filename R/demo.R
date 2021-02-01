
DEMO.dataset_folder="";

DEMO.set_dataset_folder<-function(folder){
	DEMO.dataset_folder=folder;
}

DEMO.load_data<-function(p_power,normalize=TRUE)
{

	# hard code n=100, all the example synthetic dataset have 100 samples
	n=100;
	p=2^p_power;

    matrix_folder=DEMO.dataset_folder;
	
    print(sprintf("# Reading Matrix From file: p=%d n=%d",p,n));

	filename=paste(matrix_folder , "/" , "p", p , "_n" , n , ".RData", sep = "");
	out<-get(load(filename));

	if(normalize){
		v=apply(data_set$data,1,var);
		v=1/sqrt(v);
		D=Matrix::Diagonal(n=length(v),x=v);
		data_set$data= D%*% data_set$data;
		output <- list(
			"data" = out$data, 
			"X_star" = out$X_star,
			"D" = D
		);
	}else{
		output <- list(
			"data" = out$data, 
			"X_star" = out$X_star
		);
	}

	
	return(output);
}


DEMO.lambda_search<- function(p_power ,lambda_sample=.3, K=5, criterion="AIC"){

  	# Generate data
	out<-SQUIC::DEMO.load_data( p_power = p_power );
	X_star<-out$X_star;
	data<-out$data;

	# Generate a lambda_set
	out<-SQUIC::SQUIC_S(data=data , lambda_sample=lambda_sample ,lambda_set_length=10);
	print("SQUIC::SQUIC_S");
	print(out);
	lambda_set<-out$lambda_set;

	# Do CV on for best lambda
	out<-SQUIC::SQUIC_CV(data=data , lambda_set=lambda_set , K=K , criterion=criterion );
	print("SQUIC::SQUIC_CV AIC");
	print(out);
	lambda_opt <- out$lambda_opt;
	CV_mean    <- out$CV_mean;


	f1_set	    <-replicate(length(lambda_set), 0);
	acc_set     <-replicate(length(lambda_set), 0);
	nnzpr_X_set <-replicate(length(lambda_set), 0);	
	
	for (i in 1:length(lambda_set)) {
		out<-SQUIC::DEMO.compare(alg="SQUIC" , data=data , lambda=lambda_set[i] , tol=1e-3 , max_iter=5 , X_star=X_star );
		f1_set[i]      <-out$f1;
		acc_set[i]     <-out$acc;
		nnzpr_X_set[i] <- (Matrix::nnzero(out$X)/nrow(out$X));
	}

	output <- list(
		"nnzpr_X_set"	 = nnzpr_X_set,
		"f1_set"     	 = f1_set, 
		"acc_set"    	 = acc_set,
		"lambda_opt"     = lambda_opt,				
		"lambda_set" 	 = lambda_set,
		"CV_mean"        = CV_mean
	);

	return(output);
}

DEMO.performance <- function(p_power_max,lambda=0.4,n=100,tol=1e-4,max_iter=10) 
{

	time_squic		<-replicate(p_power_max, 0);
	time_equal		<-replicate(p_power_max, 0);	
	time_quic		<-replicate(p_power_max, 0);

	for (i in 8:p_power_max) {

        p=2^i;

        # Generate data
	    out<-SQUIC::DEMO.load_data(p=p);
	    X_star<-out$X_star;
	    data_full<-out$data;

		print(sprintf("Benchmark for p=%d started",p));

		out<-SQUIC::DEMO.compare(alg="SQUIC"   , data_train=data_full , lambda=lambda , tol=tol , max_iter=max_iter , X_star=NULL);
		time_squic[i]<-out$time;

		out<-SQUIC::DEMO.compare(alg="EQUAL"   , data_train=data_full , lambda=lambda , tol=tol , max_iter=max_iter , X_star=NULL);
		time_equal[i]<-out$time;

		out<-SQUIC::DEMO.compare(alg="QUIC"    ,  data_train=data_full , lambda=lambda , tol=tol , max_iter=max_iter , X_star=NULL);
		time_quic[i]<-out$time;			
	}

	output <- list(
		"time_squic"   			= time_squic, 
		"time_equal" 			= time_equal, 		
		"time_quic" 			= time_quic 			
		)

	return(output);
}

DEMO.compare <- function(alg,data,lambda=0.5,tol=1e-4,max_iter=10, X_star= NULL) 
{
	data_t <- Matrix::t(data);

	verbose = 1;

	time_start <- Sys.time()
	if(alg=="QUIC")
	{
		print("#QUIC")
		# QUIC
		out	<-QUIC::QUIC(S=cov(data_t), rho=lambda, path = NULL, tol = tol, msg = verbose, maxIter = max_iter, X.init =NULL, W.init = NULL)
		X	<-out$X;
	}
	else if(alg=="SQUIC")
	{
		print("#SQUIC")
		# SQUIC
		out	<-SQUIC::SQUIC(Y1=data ,lambda=lambda , max_iter=max_iter , drop_tol=tol/10 , term_tol=tol , verbose=verbose );
		X	<-out$X;
	}
	else if(alg=="EQUAL")
	{
		print("#EQUAL")
		# EQUAL
		out	<-EQUAL::EQUAL(X=data_t,type=TRUE,sdiag=FALSE,lambda=lambda,lambda.min=sqrt(log(ncol(data_t))/nrow(data_t)),nlambda=1,err=tol,maxIter = max_iter,rho=1)
		X	<-out$Omega[[1]];
	}
	else
	{
		stop("Alg not found");
	};
	
	time_end	<- Sys.time()

	if(!is.null(X_star))
	{
		# Convert matrix to labels
        print("#Computing F1 Score & Accuracy")
		X_label <-  as.vector( ((X)!=0)*1 );
		X_star_label <-  as.vector( ((X_star)!=0)*1 );

		output <- list(
			"time" = time_end-time_start,
			"X"    = X, 
			"f1"   = MLmetrics::F1_Score(X_star_label,X_label, positive = "1"),
			"acc"  = MLmetrics::Accuracy(X_star_label,X_label)				
		);

	}else{

		# Convert matrix to labels
		output <- list(
			"X"    = X, 
			"time" = time_end-time_start	
		);
	}

	return(output);
}
