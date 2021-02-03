DEMO.F1 <-function(X,X_star)
{
	X@x=X@x*0+1;
	X_star@x=X_star@x*0+1;

	temp=X-X_star;

	TP=length(which(temp@x==0));
	FP_plus_FN=length(which(temp@x != 0));
		
	F1= TP / ( TP + 0.5*( FP_plus_FN ) );

	return(F1);
}


DEMO.set_dataset_folder<-function(folder)
{
	assign("SQUIC_DEMO_dataset_folder", folder, envir = .GlobalEnv)
}

DEMO.load_data<-function(p , n , normalize=TRUE)
{
    matrix_folder=SQUIC_DEMO_dataset_folder;
	
    print(sprintf("# Reading Matrix From file: p=%d n=%d",p,n));

	filename=paste(matrix_folder , "/" , "p", format(p, scientific = FALSE) , "_n" , format(n, scientific = FALSE) , ".RData", sep = "");
	out<-get(load(filename));

	if(normalize){
		v=apply(out$data,1,var);
		# R uses n-1 as denominator .. here we fore n
		v = v*(n-1)/n;
		v=1/sqrt(v);
		D=Matrix::Diagonal(n=length(v),x=v);
		out$data= as.matrix(D%*% out$data);
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


DEMO.lambda_search<- function(p,n,lambda_sample=.1, K=5){

  	# Generate data
	out=SQUIC::DEMO.load_data( p=p , n=n );
	X_star=out$X_star;
	data=out$data;

	time_CV_S = Sys.time();
	
	# Generate a lambda_set
	out=SQUIC::SQUIC_S(data=data , lambda_sample=lambda_sample ,lambda_set_length=10);
	print("SQUIC::SQUIC_S");
	lambda_set=out$lambda_set;

	# Do CV on for best lambda using AIC
	out=SQUIC::SQUIC_CV(data=data , lambda_set=lambda_set , K=K );
	print("SQUIC::SQUIC_CV");
	lambda_opt_AIC=out$lambda_opt_AIC;
	lambda_opt_BIC=out$lambda_opt_BIC;
	lambda_opt_LL=out$lambda_opt_LL;

	CV_mean_AIC=out$CV_mean_AIC;
	CV_mean_BIC=out$CV_mean_BIC;
	CV_mean_LL=out$CV_mean_LL;

	time_CV_S = as.numeric(Sys.time()- time_CV_S);

	# Compute the entire lambda path
	f1_set=replicate(length(lambda_set), 0);
	acc_set=replicate(length(lambda_set), 0);
	nnzpr_X_set=replicate(length(lambda_set), 0);	
	time_set=replicate(length(lambda_set), 0);	
	
	for (i in 1:length(lambda_set)) {
		time_set[i]= Sys.time();
		out=SQUIC::DEMO.compare(alg="SQUIC" , data=data , lambda=lambda_set[i] , tol=1e-3 , max_iter=5 , X_star=X_star );
		f1_set[i]=out$f1;
		nnzpr_X_set[i]=(Matrix::nnzero(out$X)/nrow(out$X));
		time_set[i]= as.numeric(Sys.time()- time_set[i]);
	}

	output=list(
		"time_set"        = time_set,
		"time_CV_S"       = time_CV_S,
 		"nnzpr_X_set"	  = nnzpr_X_set,
		"f1_set"     	  = f1_set, 
		"lambda_set" 	  = lambda_set,	
		"lambda_opt_AIC"  = lambda_opt_AIC,	
		"lambda_opt_BIC"  = lambda_opt_BIC,	
		"lambda_opt_LL"   = lambda_opt_LL								
		"CV_mean_AIC"	  = CV_mean_AIC,
		"CV_mean_BIC"	  = CV_mean_BIC,
		"CV_mean_LL"	  = CV_mean_LL			
	);

	return(output);
}

DEMO.performance <- function(p_set,lambda=0.4,n=100,tol=1e-4,max_iter=10) 
{

	time_squic		<-replicate(p_set, 0);
	time_equal		<-replicate(p_set, 0);	
	time_quic		<-replicate(p_set, 0);

	for (i in 1:length(p_set)) {

        p=p_set[i];

        # Generate data
	    out<-SQUIC::DEMO.load_data( p=p , n=n );
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
        print("#Computing F1-Score & Accuracy")
		F1=DEMO.F1(X,X_star);

		output <- list(
			"time" = time_end-time_start,
			"X"    = X, 
			"f1"   = F1			
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
