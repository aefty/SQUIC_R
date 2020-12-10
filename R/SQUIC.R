#dyn.load("~/libSQUIC.dylib")
usethis::use_package("Matrix") 
usethis::use_package("MASS") # Use for multivriate data generation
usethis::use_package("BigQuic") # Use for multivriate data generation

#usethis::use_package("https://cran.r-project.org/src/contrib/Archive/QUIC/QUIC_1.1.1.tar.gz", type="source")

SQUIC <- function(Data, M, lambda, max_iter, drop_tol, term_tol, X0, W0) {
  
  mode_M<-1;
  mode_X0W0<-1;
  if(is.null(M)){
      M<-Matrix::Diagonal(1);
      mode_M<-0;
    }
  if(is.null(X0) || is.null(W0)){
	X0<-Matrix::Diagonal(1);
	W0<-Matrix::Diagonal(1);
    mode_X0W0<-0;
  }
  
   .Call(`_SQUIC_SQUIC_BASE`, mode_M, mode_X0W0, Data, M, lambda, max_iter, drop_tol, term_tol, X0, W0);
}


SQUIC_demo_iid <- function(p=100, n=20, lambda=.5,max_iter=10,eps=1e-3) {

	sprintf("Generating data for test started");

	start_time <- Sys.time()
	Y=matrix( rnorm(p*p,mean<-0,sd<-1), p, n);
	end_time=Sys.time()

	sprintf("Generating data for test finished: time=%f",end_time-start_time);

	ouput<-SQUIC(Y, NULL, lambda, max_iter, drop_tol=1e-4, term_tol=eps, NULL, NULL);

	return(ouput);
}


SQUIC_demo_trid <- function(p=100, n=20, lambda=.5, M=NULL, max_iter=10,eps=1e-3) {

	sprintf("Generating data for test started");

	start_time<-Sys.time()
	iC_star <- Matrix::bandSparse(p, p,
                (-1):1,
                list(rep(-0.5, p-1), 
                     rep(1.25, p), 
                     rep(-0.5, p-1)));

	C_star <- solve(iC_star);
	mu_star <- replicate(p, 0);
 	Data <- MASS::mvrnorm(n, mu_star, C_star, tol = 1e-6, empirical = FALSE, EISPACK = FALSE);
 	Data <- t(Data);
 	end_time <- Sys.time()

	sprintf("Generating data for test finished: time=%f",end_time-start_time);

	drop_tol <- eps;
	term_tol <- eps;
	ouput <- SQUIC(Data, M, lambda, max_iter, drop_tol=1e-4, term_tol, NULL, NULL);
	return(ouput);
}


SQUIC_demo_benchmark <- function(p_power_max,lambda=0.5,n=100,eps=1e-3,max_iter=10) {

	time_squic<-replicate(p_power_max, 0);
	time_bigquic<-replicate(p_power_max, 0);

	error_iC_1<-replicate(p_power_max, 0);
	error_C_1<-replicate(p_power_max, 0);

	error_iC_0<-replicate(p_power_max, 0);
	error_C_0<-replicate(p_power_max, 0);

	for (i in 1:p_power_max) {

		p <- 2^i;

		sprintf("Generating data for test started for p=%i",p);

		start_time <- Sys.time()
		iC_star <- Matrix::bandSparse(p, p,
	                (-1):1,
	                list(rep(-0.5, p-1), 
	                     rep(1.25, p), 
	                     rep(-0.5, p-1)));

		C_star <- solve(iC_star);
		mu_star <- replicate(p, 0);
	 	Data <- MASS::mvrnorm(n, mu_star, C_star, tol = 1e-6, empirical = FALSE, EISPACK = FALSE);
	 	Data_t <- t(Data);
	 	end_time <- Sys.time()

		sprintf("Generating data for test finished: time=%f",end_time-start_time);

		bigquic_start_time <- Sys.time()
		#bigquic_output<-BigQuic::BigQuic(X = Data, inputFileName = NULL, outputFileName = NULL, lambda = lambda, numthreads = 8, maxit = max_iter, epsilon = eps, k = 0, memory_size = 8000, verbose = 1, isnormalized = 0, seed = NULL, use_ram = TRUE);
		bigquic_output<-QUIC::QUIC(S=cov(Data), rho=lambda, path = NULL, tol = eps, msg = 1, maxIter = max_iter, X.init =NULL, W.init = NULL)
		bigquic_end_time <- Sys.time()

		squic_start_time <- Sys.time()
		squic_output<-SQUIC(Data_t, NULL, lambda=lambda, max_iter=max_iter, drop_tol=1e-4, term_tol=eps, NULL, NULL);
		squic_end_time <- Sys.time()

		time_squic[i]<-squic_end_time-squic_start_time;
		time_bigquic[i]<-bigquic_end_time-bigquic_start_time;

		# L1 per row error with respect to the QUIC
		error_iC_1[i]<-Matrix::norm(squic_output$iC-bigquic_output$X,"1")/p;
		error_C_1[i]<-Matrix::norm(squic_output$C-bigquic_output$W,"1")/p;

		# Measure of the difference of the nonzero pattern ... someting like the L0 nrom
		error_iC_0[i]<-Matrix::norm( (squic_output$iC!=0)-1*(bigquic_output$X!=0),"1")/p;
		error_C_0[i]<-Matrix::norm( (squic_output$C!=0)-1*(bigquic_output$W!=0),"1")/p;
	}

	output <- list(
		"time_squic"   = time_squic, 
		"time_bigquic" = time_bigquic, 
		"error_iC_1"   = error_iC_1,
		"error_C_1"    = error_C_1,
		"error_iC_0"   = error_iC_0,
		"error_C_0"    = error_C_0		
		)

	return(output);
}
