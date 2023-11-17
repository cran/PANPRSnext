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
//' @param dims vector of dimensions
//' @param params vector of parameters
// [[Rcpp::export]]
Rcpp::List gsPEN_cpp(
    arma::Mat<double> summary_betas,
    arma::Col<int> ld_J,
    arma::Mat<int> index_matrix,
    arma::Col<int> index_J,
    arma::Col<double> ld_vec,
    arma::Mat<double> SD_vec,
    arma::Mat<double> tuning_matrix,
    arma::Col<int> dims,
    arma::Col<double> params);
