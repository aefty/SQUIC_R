SQUIC.lambda_set<-function(data_full, lambda_sample=.5,lambda_set_length=10 , M=NULL){

	# Get sample covarinace matrix by running SQUIC with max_iter=0;
	squic_output<-SQUIC(data=data_full,lambda=lambda_sample, max_iter=0, drop_tol=0, term_tol=0, M=M, X0=NULL, W0=NULL, data_test=NULL,verbose=0);

	# Get absolute value of the max and mean of nonzeros in S
	S<-squic_output$S;
	S_abs_vec = abs(S@x);

	S_abs_max<-max(S_abs_vec);
	S_abs_min<-min(S_abs_vec);
	S_abs_mean<-mean(S_abs_vec);
	S_abs_sd<-sd(S_abs_vec);
	nnz_S<-length(S_abs_vec);

	up  <- (lambda_sample + S_abs_max)/2;
	low <- (lambda_sample + S_abs_min)/2;

	delta     <- (low-up)/lambda_set_length;
	lambda_set<- seq(up , low , delta);
	
	output <- list(
		"nnz_S"			= nnz_S, 
		"S_abs_max"		= S_abs_max, 
		"S_abs_min" 	= S_abs_min, 
		"S_abs_mean"	= S_abs_mean, 
		"S_abs_sd" 		= S_abs_sd,
		"lambda_set"	= lambda_set 					
	);

	return(output);
}