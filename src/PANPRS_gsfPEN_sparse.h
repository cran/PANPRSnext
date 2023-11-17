#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Main CPP function
//' @param summary_betas matrix of summary statistics
//' @param ld_J vector of indices of SNPs in LD with the current SNP
//' @param index_matrix matrix of indices of SNPs in LD with the current SNP
//' @param index_J vector of indices of SNPs in LD with the current SNP
//' @param ld_vec vector of LD values
//' @param SD_vec matrix of SD values
//' @param tuning_matrix matrix of tuning parameters
//' @param lambda0_vec vector of lambda0 values
//' @param z_matrix matrix of z values
//' @param lambda_vec_func vector of lambda values
//' @param func_lambda matrix of lambda values
//' @param Ifunc_SNP vector of indices of SNPs in LD with the current SNP
//' @param dims vector of dimensions
//' @param params vector of parameters
// [[Rcpp::export]]
Rcpp::List gsfPEN_sparse_cpp(
    arma::Mat<double> summary_betas,
    arma::Col<int> ld_J,
    arma::Mat<int> index_matrix,
    arma::Col<int> index_J,
    arma::Col<double> ld_vec,
    arma::Mat<double> SD_vec,
    arma::Mat<double> tuning_matrix,
    arma::Col<double> lambda0_vec,
    arma::Mat<double> z_matrix,
    arma::Col<double> lambda_vec_func,
    arma::Mat<int> func_lambda,
    arma::Col<int> Ifunc_SNP,
    arma::Col<int> dims,
    arma::Col<double> params);
