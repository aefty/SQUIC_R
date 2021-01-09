SQUIC_CV_other<-function(data_full,lambda_set,K=4,tol=1e-2,max_iter=3,criterion="LL",M=NULL,X0=NULL,W0=NULL)
{

	p=nrow(data_full);
	n_full=ncol(data_full);

	#Construct active sample sets
	lambda_set<-sort(lambda_set,decreasing =TRUE)
	nlambda=length(lambda_set)
	active_set<-split(1:n_full, rep(1:K, length =n_full));

	
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

			lambda=lambda_set[l];

			if(!is.null(M)){
				M@x = M@x*lambda;
			}
			
			# run a rough (low tolerences and iterations)
			squic_output<-SQUIC::SQUIC(data_train=data_train, lambda=lambda, max_iter=max_iter, drop_tol=tol, term_tol=tol, M=M , X0=X0, W0=W0,data_test=data_test,verbose=0);

			#Extract the results form SQUIC
			X<-squic_output$X;
			logdetX<-squic_output$info_logdetX;
			trXS_test<-squic_output$info_trXS_test;

			#logliklihood (Not negative logliklihood!!!)
			logliklihood<-( p*log(2*3.14) + logdetX - trXS_test )*n_test/2;

			# For all criterion smaller is better
			if (criterion == "LL") # Defualt criterion is logliklihood (negative)
			{
				CV[k,l]=-logliklihood;
           	} 
		   	else if (criterion == "AIC") # AIC Criterion 
		   	{
                CV[k,l] = 2*sum(X!=0) - 2*logliklihood ;
            }
            else if (criterion == "BIC") # BIC Criterion
			{
                CV[k,l] = sum(X!=0)*log(n_test) - 2*logliklihood; 			
            }
			else{
				stop("Criterion specfied is not valid.")
			}

			temp_list[[l]]<-abs(X - Matrix::Diagonal(p,Matrix::diag(X)) );
		}

		U_set[[k]]<-temp_list;
	}

	# Find the smallest value in the CV matrix and select the corresponding lambda
	CV_mean<-apply(CV, 2, mean);
	l_inx<-which.min(CV_mean);
	lambda_opt<-lambda_set[l_inx];

	# Estimate M Matrix
	print("Constructing M");
	# 1) Construct   trace(U_set[[i]] U_set[[j]])
	U<-matrix(0,nrow=K,ncol=K);
	for (i in 1:K){
		for (j in 1:i){
			U[i,j]=sum( U_set[[i]][[l_inx]] * U_set[[j]][[l_inx]] );
		}
	}
	U=Matrix::forceSymmetric(U,uplo = "L");

	# 2) Eigen decompositon of the smaller matrix
	eig<-eigen(U, symmetric=TRUE, only.values = FALSE, EISPACK = FALSE);
	eig_val_inx <- which(abs(eig$values)==max(abs(eig$values)));
	eigvec<-abs(eig$vectors[,eig_val_inx]);

	M = 0;
	for (k in 1:K)
	{
		M = M+U_set[[k]][[l_inx]]* eigvec[k]/K;
	}
	
	#Clean Matrix of zeros
	M  <- drop0(M);

	if(nnzero(M)==0){
		M = NULL;
	}else{
		M@x = 1/M@x;
		M@x = M@x / max(M@x);
	}

	output <- list(
		"CV"   		= CV, 
		"lambda"    = lambda_opt,	
		"M"         = M		
		);

	return(output);
}