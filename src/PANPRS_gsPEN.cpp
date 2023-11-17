// Prevent bounds checking for increased performance
#define ARMA_NO_DEBUG 1

#include <RcppArmadillo.h>

#include <stdio.h>
#include <stdlib.h>

#include "PANPRS_gsPEN.h"

Rcpp::List gsPEN_cpp(
    arma::Mat<double> summary_betas,
    arma::Col<int> ld_J,
    arma::Mat<int> index_matrix,
    arma::Col<int> index_J,
    arma::Col<double> ld_vec,
    arma::Mat<double> SD_vec,
    arma::Mat<double> tuning_matrix,
    arma::Col<int> dims,
    arma::Col<double> params)
{
  int num_SNP = dims(0);

  int P = dims(1);
  int Q = dims(2);

  int nrow_tuning_matrix = dims(5);

  int nrow_beta_matrix = dims(7);
  int ncol_beta_matrix = dims(8);

  double upper_val = params(0);
  int num_iter = params(1);
  int breaking = params(2);
  int z_scale = params(3);
  int df_max = params(4);
  int num_indices = params(5);

  double epsilon = 0.0001;

  arma::Col<int> num_iter_vec(nrow_tuning_matrix, arma::fill::zeros);
  arma::Col<double> sum_betas(P, arma::fill::zeros);

  arma::Mat<double> beta_matrix(nrow_beta_matrix, ncol_beta_matrix, arma::fill::zeros);
  arma::Mat<double> joint_b_matrix(P, Q, arma::fill::zeros);
  arma::Mat<double> temp_b_matrix(P, Q, arma::fill::zeros);
  arma::Mat<int> skip(P, Q, arma::fill::zeros);

  GetRNGstate();

  for (int tuning_index = 0; tuning_index < nrow_tuning_matrix; tuning_index++)
  {
    double lambda1 = tuning_matrix(tuning_index, 0);
    double lambda2 = tuning_matrix(tuning_index, 2);
    double tau2 = tuning_matrix(tuning_index, 3);

    joint_b_matrix.zeros();
    temp_b_matrix.zeros();
    skip.zeros();
    sum_betas.zeros();

    bool converges = true;
    for (int n = 1; n <= num_iter; n++)
    {

      if (num_indices != 0)
      {
        for (int i = 0; i < num_indices; i++)
        {
          int j = index_J(i);

          for (int q = 0; q < Q; q++)
          {
            double bj_bar = summary_betas(j, q);
            if (bj_bar != 0.0)
            {
              double threshold = ((z_scale == 1) * SD_vec(j, q) + (z_scale == 0)) *
                                 (lambda1 + lambda2 / (sum_betas(j) + tau2));

              joint_b_matrix(j, q) = (bj_bar > threshold) * (bj_bar - threshold) +
                                     (bj_bar < -threshold) * (bj_bar + threshold);
            }

            if (summary_betas(j, q) * joint_b_matrix(j, q) < 0)
              perror("Sign inverse error");
          }
        }
      }

      // Finished initialization
      // Now, iterate over all SNPs in LD
      for (int p = 0; p < num_SNP; p++)
      {
        int j = index_matrix(p, 0);

        for (int q = 0; q < Q; q++)
        {
          if (skip(j, q) == 0)
          {
            if (summary_betas(j, q) != 0.0)
            {
              double bj_bar = summary_betas(j, q);

              for (int i = index_matrix(p, 1); i <= index_matrix(p, 2); i++)
              {
                bj_bar -= ld_vec(i) * joint_b_matrix(ld_J(i), q);
              }

              double threshold = ((z_scale == 1) * SD_vec(j, q) + (z_scale == 0)) *
                                 (lambda1 + lambda2 / (sum_betas(j) + tau2));

              if (fabs(bj_bar) > upper_val)
              {
                if (breaking == 1)
                {
                  num_iter_vec(tuning_index) = -1;
                  break;
                }
                else
                {
                  bj_bar = 0.0;
                  skip(j, q) = 1;
                }
              }

              joint_b_matrix(j, q) = (bj_bar > threshold) * (bj_bar - threshold) +
                                     (bj_bar < -threshold) * (bj_bar + threshold);
            }
            else
            {
              joint_b_matrix(j, q) = 0.0;
            }
          }
          else
          {
            joint_b_matrix(j, q) = 0.0;
          }
        }
      }

      // Update the skip matrix
      skip.elem(
              arma::find(
                  arma::abs(joint_b_matrix) > upper_val))
          .ones();

      // Check for divergence
      int df_q = arma::accu(
          arma::abs(joint_b_matrix(arma::span::all, 0)) != 0.0);
      if (df_q > df_max)
        converges = false;

      if (!converges)
      {
        num_iter_vec(tuning_index) = -2;
        break;
      }

      // Check if the beta matrix has converged
      bool found = arma::any(
          arma::find(
              arma::abs(temp_b_matrix - joint_b_matrix) > epsilon));

      if (!found)
      {
        // Determine which entries need to update
        arma::uvec beta_indices = arma::find(joint_b_matrix != 0.0);
        beta_indices.for_each(
            [&](arma::uword idx)
            {
              arma::uvec beta_subscript = arma::ind2sub(arma::size(joint_b_matrix), idx);
              beta_matrix(tuning_index, idx) = joint_b_matrix(beta_subscript(0), beta_subscript(1));
            });

        num_iter_vec(tuning_index) = n;
        break;
      }

      // Update for next iteration if we haven't converged yet
      temp_b_matrix(arma::span::all, arma::span::all) = joint_b_matrix(arma::span::all, arma::span::all);
      sum_betas(arma::span::all) = arma::sum(arma::abs(joint_b_matrix(arma::span::all, arma::span::all)), 1);

      num_iter_vec(tuning_index) = n;
    }
  }
  PutRNGstate();

  return Rcpp::List::create(
      Rcpp::Named("beta_matrix") = beta_matrix,
      Rcpp::Named("num_iter_vec") = num_iter_vec,
      Rcpp::Named("tuning_matrix") = tuning_matrix);
}
