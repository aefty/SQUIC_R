SQUIC_DEMO.data<-function(type="trid",p_power=5,n=100)
{

	set.seed(1);

    start_time <- Sys.time()
    p=4^p_power;
	
    print(sprintf("# Generating Percision Matrix: type=%s p=%d n=%d",type,p,n));

	if(type=="eye") # Idendity Matrix for iC_star
	{
			iC_star <- Matrix::Diagonal(p);
			C_star <- solve(iC_star);
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
			C_star <- solve(iC_star);	
	}
		else if(type=="rand")  # Random matrix for iC_star (averag of 5 nnz per row) 
	{
			nnz_per_row=5;

			# Make PSD symmetric Random Matrix
			iC_star <-Matrix::rsparsematrix(p,p,NULL,nnz_per_row*p/2,symmetric=TRUE);
			x=Matrix::colSums(abs(iC_star))+1;
			D=Matrix::Diagonal(p,x);
			iC_star<- iC_star+D;

			# invert it for covariance
			C_star <- solve(iC_star);

			# Extract the diagional to make unit variance for C_star
			# also use the transformation for iC_star
			x=sqrt(diag(C_star));
			ix=1/x;
			
			D=Matrix::Diagonal(p,ix);
			C_star<-D%*%C_star%*%D

			D=Matrix::Diagonal(p,x);
			iC_star<-D%*%iC_star%*%D

	}else{
			stop("Unknown matrix type.")
	}

	# Generate data
    print(sprintf("# Generating data ...",type,p,n,normalized));
	mu_star <- replicate(p, 0);
 	data <- MASS::mvrnorm(n, mu_star, C_star, tol = 1e-2, empirical = FALSE, EISPACK = FALSE);
	data <- Matrix::t(data);

	finish_time <- Sys.time()
	print(sprintf("# Generating data finished: time=%f",finish_time-start_time));

	output <- list(
		"data" = data, 
		"X_star" = iC_star,
		"W_star" = C_star
	);

	return(output);
}
