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


DEMO.performance <- function(type="trid",lambda=0.4,n=100,tol=1e-4,max_iter=10) 
{

    #Hard coded values
    p_power_max<-5;

	time_squic		<-replicate(p_power_max, 0);
	time_equal		<-replicate(p_power_max, 0);	
	time_bigquic	<-replicate(p_power_max, 0);
	time_quic		<-replicate(p_power_max, 0);

	for (i in 1:p_power_max) {

        p=4^i;

        # Generate data
	    out<-SQUIC::DEMO.make_data(type=type ,p=p ,n=n ,normalized=TRUE);
	    X_star<-out$X_star;
	    data_full<-out$data;

		print(sprintf("Benchmark for p=%d started",p));

		out<-SQUIC::DEMO.compare(alg="SQUIC"   , data_full=data_full , lambda=lambda , tol=tol , max_iter=max_iter , X_star=NULL);
		time_squic[i]<-out$time;

		out<-SQUIC::DEMO.compare(alg="EQUAL"   , data_full=data_full , lambda=lambda , tol=tol , max_iter=max_iter , X_star=NULL);
		time_equal[i]<-out$time;

		out<-SQUIC::DEMO.compare(alg="QUIC"    ,  data_full=data_full , lambda=lambda , tol=tol , max_iter=max_iter , X_star=NULL);
		time_quic[i]<-out$time;			
	}

	output <- list(
		"time_squic"   			= time_squic, 
		"time_equal" 			= time_equal, 		
		"time_quic" 			= time_quic 			
		)

	return(output);
}

DEMO.compare <- function(alg,data_full,lambda=0.5,tol=1e-4,max_iter=10, X_star= NULL) 
{
	data_full_t <- Matrix::t(data_full);

	 verbose = 1;

	time_start <- Sys.time()
	if(alg=="QUIC")
	{
		print("#QUIC")
		# QUIC
		out	<-QUIC::QUIC(S=cov(data_full_t), rho=lambda, path = NULL, tol = tol, msg = verbose, maxIter = max_iter, X.init =NULL, W.init = NULL)
		X	<-out$X;
	}
	else if(alg=="SQUIC")
	{
		print("#SQUIC")
		# SQUIC
		out	<-SQUIC::SQUIC(data_train=data_full ,lambda=lambda , max_iter=max_iter , drop_tol=tol/2 , term_tol=tol , verbose=verbose );
		X	<-out$X;
	}
	else if(alg=="EQUAL")
	{
		print("#EQUAL")
		# EQUAL
		out	<-EQUAL::EQUAL(X=data_full_t,type=TRUE,sdiag=FALSE,lambda=lambda,lambda.min=sqrt(log(ncol(data_full_t))/nrow(data_full_t)),nlambda=1,err=tol,maxIter =max_iter,rho=1)
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
