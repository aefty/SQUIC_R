// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <stdio.h>
#include <assert.h>
#include <algorithm>

// SQUIC is fixed with type: long in (This is a requirement for Cholmod)
#define integer long
// SQUIC Library iterface
extern "C"
{

    void SQUIC_CPP(
        int mode,
        integer p,
        integer n_train, double *Y_train,
        integer n_test, double *Y_test,
        double lambda,
        integer *M_rinx, integer *M_cptr, double *M_val, integer M_nnz,
        int max_iter, double drop_tol, double term_tol, int verbose,
        integer *&X_rinx, integer *&X_cptr, double *&X_val, integer &X_nnz,
        integer *&W_rinx, integer *&W_cptr, double *&W_val, integer &W_nnz,
        int &info_num_iter,
        double *&info_times,     //length must be 6: [time_total,time_impcov,time_optimz,time_factor,time_aprinv,time_updte]
        double *&info_objective, // length must be size max_iter
        double &info_dgap,
        double &info_logdetx,
        double &info_trXS_test)
}

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
List SQUIC_R(arma::mat &data_train, double lambda, int max_iter, double drop_tol, double term_tol, int verbose, int mode, arma::sp_mat &M, arma::sp_mat &X0, arma::sp_mat &W0, arma::mat &data_test)
{

    // Set SQUIC execution style
    bool EXECSTYLE_M_provided = M.n_nonzero > 0;
    bool EXECSTYLE_X0W0_provided = (X0.n_nonzero > 0) && (W0.n_nonzero > 0);
    bool EXECSTYLE_data_test_provided = data_test.n_rows == data_train.n_rows;

    // Get the key size parameters
    integer p = data_train.n_rows;
    integer n_train = data_train.n_cols;
    integer n_test = -1;

    if (p < 2) // Only work iwth matrices
    {
        stop(" The number of random variables 'p' must be larger than 1.");
    }

    // Make matrix pointers
    integer *M_i;
    integer *M_j;
    double *M_val;
    integer M_nnz = 0;

    integer *X_i;
    integer *X_j;
    double *X_val;
    integer X_nnz = 0;

    integer *W_i;
    integer *W_j;
    double *W_val;
    integer W_nnz = 0;

    if (EXECSTYLE_data_test_provided) //We are going to use the test data
    {
        n_test = data_test.n_cols;

        if (data_test.n_rows != data_train.n_rows) // Matricies must be the correct dimension
        {
            stop(" The number of random variables 'p' must equal for both training and testing datasets.");
        }

        if (n_test < 1) // Matricies must be the correct dimension
        {
            stop(" The testing datasets is empty.");
        }
    }

    if (EXECSTYLE_M_provided) //We are given the M matrix
    {
        M.sync();
        M_nnz = M.n_nonzero;

        // ERROR CHECKS
        {
            // 1) Symmetry check done outside in R code.
            // 2) Matricies must be the correct dimension
            if (p != M.n_rows && p != M.n_cols)
            {
                stop("Matrix M must be of size (pxp) p=%d", p);
            }
        }

        //COPY M MATRIX
        auto M_row_indices = arma::access::rwp(M.row_indices);
        auto M_col_ptrs = arma::access::rwp(M.col_ptrs);
        auto M_values = arma::access::rwp(M.values);
        {
            M_i = new integer[M_nnz];
            M_j = new integer[p + 1];
            M_val = new double[M_nnz];

            for (integer i = 0; i < M_nnz; ++i)
            {
                M_i[i] = M_row_indices[i];
                M_val[i] = M_values[i];
            }

            for (integer i = 0; i < p + 1; ++i)
            {
                M_j[i] = M_col_ptrs[i];
            }
        }
    }

    if (EXECSTYLE_X0W0_provided) // We are using BOTH X0 and W0
    {
        X0.sync();
        W0.sync();
        X_nnz = X0.n_nonzero;
        W_nnz = W0.n_nonzero;

        // ERROR CHECKS
        {
            // 1) Symmetry check done outside in R code.
            // 2) Matricies must be the correct dimension
            if (p != X0.n_rows && p != X0.n_cols)
            {

                if (EXECSTYLE_M_provided) // failed error check, delete previous allocated memory
                {
                    delete[] M_i;
                    delete[] M_j;
                    delete[] M_val;
                }
                stop("Matrix X0 must be of size (pxp) p=%d.", p);
            }
            if (p != W0.n_rows && p != W0.n_cols)
            {
                if (EXECSTYLE_M_provided) // failed error check, delete previous allocated memory
                {
                    delete[] M_i;
                    delete[] M_j;
                    delete[] M_val;
                }
                stop("Matrix X0 must be of size (pxp) p=%d.", p);
            }

            if (X_nnz < p)
            {
                if (EXECSTYLE_M_provided)
                {
                    delete[] M_i;
                    delete[] M_j;
                    delete[] M_val;
                }
                stop("Matrix X0 must be at least diagional (X0.nnz >= p).");
            }

            if (W_nnz < p)
            {
                if (EXECSTYLE_M_provided)
                {
                    delete[] M_i;
                    delete[] M_j;
                    delete[] M_val;
                }
                stop("Matrix W0 must be at least diagional (W0.nnz >= p).");
            }
        }

        // COPY X0 MATRIX
        auto X_row_indices = arma::access::rwp(X0.row_indices);
        auto X_col_ptrs = arma::access::rwp(X0.col_ptrs);
        auto X_values = arma::access::rwp(X0.values);
        {
            X_i = new integer[X_nnz];
            X_j = new integer[p + 1];
            X_val = new double[X_nnz];

            for (integer i = 0; i < X_nnz; ++i)
            {
                X_i[i] = X_row_indices[i];
                X_val[i] = X_values[i];
            }

            for (integer i = 0; i < p + 1; ++i)
            {
                X_j[i] = X_col_ptrs[i];
            }
        }

        // COPY W0 MATRIX
        auto W_row_indices = arma::access::rwp(W0.row_indices);
        auto W_col_ptrs = arma::access::rwp(W0.col_ptrs);
        auto W_values = arma::access::rwp(W0.values);
        {
            W_i = new integer[W_nnz];
            W_j = new integer[p + 1];
            W_val = new double[W_nnz];

            for (integer i = 0; i < W_nnz; ++i)
            {
                W_i[i] = W_row_indices[i];
                W_val[i] = W_values[i];
            }

            for (integer i = 0; i < p + 1; ++i)
            {
                W_j[i] = W_col_ptrs[i];
            }
        }
    }

    // Default Result Values
    int info_num_iter = -1;                                            // Number of Newten steps required by SQUIC
    double info_dgap = -1e-12;                                         // Duality Gap between primal and dual
    double info_logdetx = -1e-12;                                      // Can be used for likelilook (AIC or BIC) computation of test data
    double info_trXS_test = -1e-12;                                    // Can be used for likelilook (AIC or BIC) computation of test data
    double *info_times_buffer = new double[6];                         // This need to be of size 6
    double *info_objective_buffer = new double[std::max(1, max_iter)]; // The objective value list, must be of size max(max_iter,1). If max_iter=0, we still keep this with size of 1

    // Run SQUIC
    SQUIC_CPP(
        mode,
        p,
        n_train, data_train.memptr(),
        n_test, data_test.memptr(),
        lambda,
        M_i, M_j, M_val, M_nnz,
        max_iter, drop_tol, term_tol, verbose,
        X_i, X_j, X_val, X_nnz,
        W_i, W_j, W_val, W_nnz,
        info_num_iter,
        info_times_buffer,
        info_objective_buffer,
        info_dgap,
        info_logdetx,
        info_trXS_test);

    // Copy data it standard format
    // In order to access the internal arrays of the SpMat class call .sync()
    arma::SpMat<double> iC(p, p);
    arma::SpMat<double> C(p, p);
    iC.sync();
    C.sync();

    // Making space for the elements
    iC.mem_resize(X_nnz);
    C.mem_resize(W_nnz);

    // Copying elements
    std::copy(X_i, X_i + X_nnz, arma::access::rwp(iC.row_indices));
    std::copy(X_j, X_j + p + 1, arma::access::rwp(iC.col_ptrs));
    std::copy(X_val, X_val + X_nnz, arma::access::rwp(iC.values));
    arma::access::rw(iC.n_rows) = p;
    arma::access::rw(iC.n_cols) = p;
    arma::access::rw(iC.n_nonzero) = X_nnz;

    std::copy(W_i, W_i + W_nnz, arma::access::rwp(C.row_indices));
    std::copy(W_j, W_j + p + 1, arma::access::rwp(C.col_ptrs));
    std::copy(W_val, W_val + W_nnz, arma::access::rwp(C.values));
    arma::access::rw(C.n_rows) = p;
    arma::access::rw(C.n_cols) = p;
    arma::access::rw(C.n_nonzero) = W_nnz;

    Rcpp::List output;

    if (max_iter == 0) // Special max_iter==0: SQUIC only compute the sparse sample covariance S
    {
        output = Rcpp::List::create(
            Named("S") = C,
            Named("info_time_total") = info_times_buffer[0],
            Named("info_time_impcov") = info_times_buffer[1]);
    }
    else // Regular case return all values
    {

        //Copy info_objective_buffer keeping only info_num_iter elements
        arma::Col<double> info_objective(info_objective_buffer, info_num_iter);

        if (EXECSTYLE_data_test_provided)
        {
            output = Rcpp::List::create(
                Named("X") = iC,
                Named("W") = C,
                Named("info_time_total") = info_times_buffer[0],
                Named("info_time_impcov") = info_times_buffer[1],
                Named("info_time_optimz") = info_times_buffer[2],
                Named("info_time_factor") = info_times_buffer[3],
                Named("info_time_aprinv") = info_times_buffer[4],
                Named("info_time_updte") = info_times_buffer[5],
                Named("info_objective") = info_objective,
                Named("info_duality_gap") = info_dgap,
                Named("info_logdetX") = info_logdetx,
                Named("info_trXS_test") = info_trXS_test);
        }
        {
            // no info_trXS_test ouput
            output = Rcpp::List::create(
                Named("X") = iC,
                Named("W") = C,
                Named("info_time_total") = info_times_buffer[0],
                Named("info_time_impcov") = info_times_buffer[1],
                Named("info_time_optimz") = info_times_buffer[2],
                Named("info_time_factor") = info_times_buffer[3],
                Named("info_time_aprinv") = info_times_buffer[4],
                Named("info_time_updte") = info_times_buffer[5],
                Named("info_objective") = info_objective,
                Named("info_duality_gap") = info_dgap,
                Named("info_logdetX") = info_logdetx);
        }
    }

    // Delete Buffers
    if (EXECSTYLE_M_provided)
    {
        delete[] M_i;
        delete[] M_j;
        delete[] M_val;
    }

    if (EXECSTYLE_X0W0_provided)
    {
        delete[] X_i;
        delete[] X_j;
        delete[] X_val;

        delete[] W_i;
        delete[] W_j;
        delete[] W_val;
    }

    delete[] info_times_buffer;
    delete[] info_objective_buffer;

    return output;
}