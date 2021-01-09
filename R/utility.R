make_data<-function(type="trid",p=10,n=5,normalized=FALSE,matrix_folder="")
{

	set.seed(10);
	start_time <- Sys.time();
	print(sprintf("# Generating data started: type=%s p=%d n=%d normalized=%d",type,p,n,normalized));

	if(matrix_folder!="")
	{

		iC_file_name=paste(matrix_folder,type,p,"iC",".rmat",sep = "");
		C_file_name=paste(matrix_folder,type,p,"iC",".rmat",sep = "");

		iC_star=readMM(iC_file_name)
		C_star=readMM(C_file_name)

	}
	else
	{
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
		}

	# Generate data
	mu_star <- replicate(p, 0);
 	data <- MASS::mvrnorm(n, mu_star, C_star, tol = 1e-6, empirical = FALSE, EISPACK = FALSE);
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
		"iC_star" = iC_star,
		"C_star" = C_star
	);

	return(output);
}

