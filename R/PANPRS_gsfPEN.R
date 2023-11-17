#' Run the gsfPEN algorithm for multiple traits, with functional annotations.
#' @param summary_z A matrix of summary statistics for each SNP and trait.
#' @param n_vec A vector of sample sizes for each of the Q traits corresponding to the Q columns of summary_z.
#' @param plinkLD A matrix of LD values for each pair of SNPs.
#' @param func_index A matrix of functional annotations for each SNP and trait. For the element at i-th row, j-th column, the entry 0 means SNP i without j-th functional annotation; entry 1 means otherwise.
#' @param n_iter The number of iterations to run the algorithm.
#' @param upper_val The upper bound for the tuning parameter.
#' @param breaking The number of iterations to run before checking for convergence.
#' @param z_scale The scaling factor for the summary statistics.
#' @param tuning_matrix A matrix of tuning parameters.
#' @param p_threshold A vector of p-values to use for the tuning parameters.
#' @param p_threshold_params A vector of parameters to use for the p-value tuning parameters.
#' @param tau_factor A vector of factors to multiply the median value by to get the tuning parameters.
#' @param sub_tuning The number of tuning parameters to use for the second iteration.
#' @param lim_lambda The range of tuning parameters to use for the first iteration.
#' @param len_lambda The number of tuning parameters to use for the second iteration.
#' @param lambda_vec A vector of tuning parameters to use for the first iteration.
#' @param lambda_vec_limit_len The number of tuning parameters to use for the first iteration.
#' @param df_max The maximum degrees of freedom for the model.
#' @param sparse_beta Whether to use the sparse version of the algorithm.
#' @param debug_output Whether to output the tuning combinations that did not converge.
#' @param verbose Whether to output information through the evaluation of the algorithm.
#' @return A named list containing the following elements:
#' beta_matrix: A matrix of the estimated coefficients for each SNP and trait.
#' num_iter_vec: A vector of the number of iterations for each tuning combination.
#' all_tuning_matrix: A matrix of the tuning parameters used for each tuning combination.
#' @importFrom stats median qnorm quantile
#' @examples
#' # Load the library and data
#' library(PANPRSnext)
#' data("summaryZ")
#' data("Nvec")
#' data("plinkLD")
#' data("funcIndex")
#'
#' # Take random subset of the data
#' subset <- sample(nrow(summaryZ), 100)
#' subset_summary_z <- summaryZ[subset, ]
#' subset_func_index <- funcIndex[subset, ]
#'
#' # Run gsfPEN
#' output <- gsfPEN_R(
#'   summary_z = subset_summary_z,
#'   n_vec = Nvec,
#'   plinkLD = plinkLD,
#'   func_index = subset_func_index
#' )
#' @export
gsfPEN_R <- function(
    summary_z,
    n_vec,
    plinkLD,
    func_index,
    n_iter = 1000,
    upper_val = NULL,
    breaking = 1,
    z_scale = 1,
    tuning_matrix = NULL,
    p_threshold = NULL,
    p_threshold_params = c(0.5, 10^-4, 4),
    tau_factor = c(1 / 25, 1, 3),
    sub_tuning = 4,
    lim_lambda = c(0.5, 0.9),
    len_lambda = 4,
    lambda_vec = NULL,
    lambda_vec_limit_len = c(1.5, 3),
    df_max = NULL,
    sparse_beta = FALSE,
    debug_output = FALSE,
    verbose = FALSE) {
  if (z_scale != 1) {
    stop("Tuning values set-up for multiple traits analysis requires z_scale = 1") # nolint: object_usage_linter.
  }

  num_func <- ncol(func_index)

  P <- nrow(summary_z)
  Q <- length(n_vec)

  summary_betas <- matrix(0, nrow = P, ncol = Q)
  SD_vec <- matrix(0, nrow = P, ncol = Q)

  for (i in 1:Q) {
    summary_betas[, i] <- summary_z[, i] / sqrt(n_vec[i])
    SD_vec[, i] <- 1 / sqrt(n_vec[i])
  }

  rownames(summary_betas) <- rownames(summary_z)

  if (is.null(df_max)) {
    df_max <- ceiling(0.7 * P)
  }

  if (is.null(p_threshold)) {
    p_threshold <- seq(p_threshold_params[1], p_threshold_params[2], length.out = p_threshold_params[3])
  }

  if (any(c(is.null(tuning_matrix), is.null(lambda_vec)))) {
    median_val <- median(apply(abs(summary_betas), 1, sum), na.rm = TRUE)
    tau_vec <- sort(median_val * tau_factor)
    lim_lambda <- quantile(abs(summary_z[, 1]), lim_lambda)

    output <- Tuning_setup_group_func(
      lambda_vec,
      lambda_vec_limit_len,
      #   p_threshold,
      num_func,
      tau_vec,
      sub_tuning,
      lim_lambda,
      len_lambda,
      # len_lim_lambda,
      median_val,
      equal_lambda = TRUE,
      Q
    )

    func_lambda <- output$func_lambda
    lambda_vec <- output$lambda_vec
    tuning_matrix <- output$tuning_matrix

    rm(output)
  } else {
    func_lambda <- permutations( # nolint: object_usage_linter.
      length(lambda_vec),
      num_func,
      repeats.allowed = TRUE
    ) - 1
  }

  beta_index <- c(seq_len(nrow(summary_betas))) - 1
  SNP_names <- rownames(summary_betas)

  # Takes only the SNPs that are present in both the PLINK data set and the GWAs
  ld_J <- PlinkLD_transform(plinkLD, SNP_names)

  J_index_matrix <- matrix(nrow = nrow(ld_J), ncol = 2)

  mat1 <- match(ld_J[, 1], SNP_names)
  mat2 <- match(ld_J[, 2], SNP_names)

  J_index_matrix[, 1] <- beta_index[mat1]
  J_index_matrix[, 2] <- beta_index[mat2]

  ld_J[, 1] <- J_index_matrix[, 1]
  ld_J[, 2] <- J_index_matrix[, 2]

  od <- order(J_index_matrix[, 1], J_index_matrix[, 2], decreasing = FALSE)
  ld_J <- ld_J[od, ]


  wind <- which(!beta_index %in% ld_J[, 1])
  num_indices <- length(wind)

  index_J <- -1
  if (num_indices > 0) {
    index_J <- beta_index[wind]
  }

  counts <- table(ld_J[, 1])
  num_SNP <- length(counts)

  index_S <- c(0, cumsum(counts)[-num_SNP])
  index_E <- cumsum(counts) - 1

  index_matrix <- matrix(nrow = num_SNP, ncol = 3)
  index_matrix[, 1] <- as.numeric(names(counts))

  index_matrix[, 2] <- index_S
  index_matrix[, 3] <- index_E

  ld_vec <- ld_J[, 3]
  ld_J <- ld_J[, 2]


  if (is.null(upper_val)) {
    upper_val <- ceiling(max(abs(summary_betas), na.rm = TRUE) * 50)
  }

  if (nrow(func_index) != nrow(summary_betas)) {
    stop("nrow of summary_betas and row of func_index do not match")
  }


  z_matrix <- as.matrix(1 - func_index)

  sum_func_index <- apply(z_matrix, 1, sum)
  Ifunc_SNP <- rep(0, P)
  Ifunc_SNP[which(sum_func_index != 0)] <- 1

  lambda0_vec <- abs(-qnorm(p_threshold / 2))

  nrow_index_matrix <- nrow(index_matrix)
  ncol_index_matrix <- ncol(index_matrix)

  nrow_z_matrix <- nrow(z_matrix)
  ncol_z_matrix <- ncol(z_matrix)

  nrow_func_lambda <- nrow(func_lambda)
  ncol_func_lambda <- ncol(func_lambda)

  nrow_tuning_matrix <- nrow(tuning_matrix)
  ncol_tuning_matrix <- ncol(tuning_matrix)

  nrow_all_tuning_matrix <- nrow_tuning_matrix * length(p_threshold) * nrow_func_lambda
  ncol_all_tuning_matrix <- (num_func + 1) + 2

  nrow_beta_matrix <- nrow_all_tuning_matrix
  ncol_beta_matrix <- P * Q

  dims <- c(
    num_SNP, # 1
    P, # 2
    Q, # 3
    nrow_index_matrix, # 4
    ncol_index_matrix, # 5
    nrow_z_matrix, # 6
    ncol_z_matrix, # 7
    nrow_func_lambda, # 8
    ncol_func_lambda, # 9
    nrow_tuning_matrix, # 10
    ncol_tuning_matrix, # 11
    nrow_all_tuning_matrix, # 12 == num_tuning
    ncol_all_tuning_matrix, # 13
    nrow_beta_matrix, # 14
    ncol_beta_matrix # 15
  )

  params <- c(
    upper_val, # 1
    n_iter, # 2
    breaking, # 3
    z_scale, # 4
    df_max, # 5
    length(p_threshold), # 6
    num_indices # 7
  )

  args <- list(
    summary_betas,
    ld_J,
    index_matrix,
    index_J,
    ld_vec,
    SD_vec,
    tuning_matrix,
    lambda0_vec,
    z_matrix,
    lambda_vec,
    func_lambda,
    Ifunc_SNP,
    dims,
    params
  )

  if (verbose) {
    cat("Number of total tuning combinations =", nrow_all_tuning_matrix, "\n")
  }

  if (sparse_beta) {
    Z <- do.call(
      gsfPEN_sparse_cpp,
      args
    )
  } else {
    Z <- do.call(
      gsfPEN_cpp,
      args
    )
  }

  beta_matrix <- Z$beta_matrix
  colnames(beta_matrix) <- paste0(
    rep(SNP_names, times = Q),
    ".trait",
    rep(c(1:Q), each = P)
  )

  all_tuning_matrix <- Z$all_tuning_matrix
  colnames(all_tuning_matrix) <- c(
    "lambda0",
    paste0("lambdaf", c(1:num_func)),
    "lambda2",
    "tau2"
  )

  num_iter_vec <- Z$num_iter_vec

  # Remove the tuning combinations that did not converge (correspons to -2 in num_iter_vec)
  if (!debug_output) {
    converge_index <- which(num_iter_vec > 0)
    if (verbose) {
      cat("Removing (", length(num_iter_vec) - length(converge_index), ") tuning combinations that did not converge\n")
    }
    num_iter_vec <- num_iter_vec[converge_index]
    beta_matrix <- beta_matrix[converge_index, ]
    all_tuning_matrix <- all_tuning_matrix[converge_index, ]
  }

  output <- Clean_results(
    beta_matrix = beta_matrix,
    num_iter_vec = num_iter_vec,
    all_tuning_matrix = all_tuning_matrix
  )

  return(output)
}
