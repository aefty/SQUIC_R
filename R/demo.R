DEMO.make_data<-function(type="trid",p=4^5,n=100)
{

	set.seed(1);

    start_time <- Sys.time()
	
    print(sprintf("# Generating Percision Matrix: type=%s p=%d n=%d",type,p,n));

	if(type=="eye") # Idendity Matrix for iC_star
	{
		iC_star <- Matrix::Diagonal(p);
	}
	else if(type=="trid") # Tridiagiaonl matrix for iC_star 
	{
		iC_star <- Matrix::bandSparse(p, p,
				(-1):1,
				list(rep(-6/9, p-1), 
					rep(10/6, p), 
					rep(-6/9, p-1)));
		iC_star[1,1]=4/3;
		iC_star[p,p]=4/3;
	}
		else if(type=="rand")  # Random matrix for iC_star (averag of 5 nnz per row) 
	{
		nnz_per_row=5;

		# Make PSD symmetric Random Matrix
		iC_star <-Matrix::rsparsematrix(p,p,NULL,nnz_per_row*p/2,symmetric=TRUE);
		x=Matrix::colSums(abs(iC_star))+1;
		D=Matrix::Diagonal(p,x);
		iC_star<- iC_star+D;

	}else{
		stop("Unknown matrix type.")
	}

	# Generate data
	z    <- replicate(n,rnorm(p));
	iC_L <- chol(iC_star);
	data <- matrix(solve(iC_L,z),p,n);

	finish_time <- Sys.time()
	print(sprintf("# Generating Data: time=%f",finish_time-start_time));

	output <- list(
		"data" = data, 
		"X_star" = iC_star
	);

	return(output);
}


DEMO.make_data_and_save<-function(){

	type="rand";
	n=100;

	for (i in 1:6) {
		p=4^i;
		print(sprintf("# Making dataset file: type=%s p=%d n=%d",type,p,n));
		data<-SQUIC::DEMO.make_data(p=p,type=type,n=n);
		file_name=paste("dataset_",type,"_",p,"_",n,".RData",sep = "");
		save(data,file=file_name, version = 2);
	}


	type="eye";
	n=100;

	for (i in 1:6) {
		p=4^i;
		print(sprintf("# Making dataset file: type=%s p=%d n=%d",type,p,n));
		data<-SQUIC::DEMO.make_data(p=p,type=type,n=n);
		file_name=paste("dataset_",type,"_",p,"_",n,".RData",sep = "");
		save(data,file=file_name, version = 2);
	}


	type="trid";
	n=100;
	for (i in 1:6) {
		p=4^i;
		print(sprintf("# Making dataset file: type=%s p=%d n=%d",type,p,n));
		data<-SQUIC::DEMO.make_data(p=p,type=type,n=n);
		file_name=paste("dataset_",type,"_",p,"_",n,".RData",sep = "");
		save(data,file=file_name, version = 2);
	}

}


DEMO.load_data<-function(type="trid",p=4^5,n=100)
{
    matrix_folder=system.file("extdata",package = "SQUIC")
	
    print(sprintf("# Reading Matrix From file: type=%s p=%d n=%d",type,p,n));

	filename=paste(matrix_folder , "/" , "dataset_" , type , "_" , p , "_" , n , ".RData", sep = "");
	out<-get(load(filename));

	output <- list(
		"data" = out$data, 
		"X_star" = out$X_star
	);

	return(output);
}

DEMO.lambda_search<- function(type="trid", p=4^5 , n=100 , lambda_sample=.4){

  	# Generate data
	out<-SQUIC::DEMO.load_data(type=type , p=p ,n=n );
	X_star<-out$X_star;
	data<-out$data;

	# Generate a lambda_set
	out<-SQUIC::SQUIC_S(data=data , lambda_sample=lambda_sample ,lambda_set_length=10);
	print("SQUIC::SQUIC_S");
	print(out);
	lambda_set<-out$lambda_set;

	# Do CV on for best lambda
	out<-SQUIC::SQUIC_CV(data=data , lambda_set=lambda_set , criterion="LL" );
	print("SQUIC::SQUIC_CV LL");
	print(out);
	lambda_opt_LL=out$lambda_opt;

	out<-SQUIC::SQUIC_CV(data=data , lambda_set=lambda_set , criterion="AIC" );
	print("SQUIC::SQUIC_CV AIC");
	print(out);
	lambda_opt_AIC=out$lambda_opt;


	out<-SQUIC::SQUIC_CV(data=data , lambda_set=lambda_set , criterion="BIC" );
	print("SQUIC::SQUIC_CV BIC");
	print(out);
	lambda_opt_BIC=out$lambda_opt;

	f1_set	<-replicate(length(lambda_set), 0);
	acc_set <-replicate(length(lambda_set), 0);
	nnzpr_X_set <-replicate(length(lambda_set), 0);	
	
	for (i in 1:length(lambda_set)) {
		out<-SQUIC::DEMO.compare(alg=alg , data=data , lambda=lambda_set[i] , tol=1e-4 , max_iter=10 , X_star=X_star);
		f1_set[i]<-out$f1;
		acc_set[i]<-out$acc;
		nnzpr_X_set[i]<-Matrix::nnzero(out$X)/p;
	}

	output <- list(
		"nnzpr_X_set"	 = nnzpr_X_set,
		"f1_set"     	 = f1_set, 
		"acc_set"    	 = acc_set,
		"lambda_opt_LL"  = lambda_opt_LL,
		"lambda_opt_AIC" = lambda_opt_AIC,
		"lambda_opt_BIC" = lambda_opt_BIC,				
		"lambda_set" 	 = lambda_set
	);


	return(output);
}

DEMO.performance <- function(type="trid",lambda=0.4,n=100,tol=1e-4,max_iter=10) 
{
    #Hard coded values
    p_power_max		<-6;

	time_squic		<-replicate(p_power_max, 0);
	time_equal		<-replicate(p_power_max, 0);	
	time_quic		<-replicate(p_power_max, 0);

	for (i in 1:p_power_max) {

        p=4^i;

        # Generate data
	    out<-SQUIC::DEMO.load_data(type=type , p=p ,n=n );
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
		out	<-SQUIC::SQUIC(Y1=data ,lambda=lambda , max_iter=max_iter , drop_tol=tol/2 , term_tol=tol , verbose=verbose );
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
