Tuning_setup_group_only <- function(
    tau_vec,
    sub_tuning,
    lim_lambda,
    len_lambda,
    median_val,
    Q) {
  print(paste0(Q))
  lambda_vec <- seq(min(lim_lambda), max(lim_lambda), len = len_lambda)

  if (Q == 1) {
    ratios <- 1
  } else {
    ratios <- seq(1, 0, len = sub_tuning)
  }
  ratios <- 1

  start <- TRUE
  for (i in seq_along(tau_vec)) {
    tau <- tau_vec[i]

    for (t in 2:length(lambda_vec)) {
      total <- lambda_vec[t]
      lower <- total * ratios
      upper <- total * (1 - ratios)
      for (idx in seq_along(lower)) {
        temp <- cbind(lower[idx], 0, upper[idx] * tau, tau)
        if (start) {
          tuning_matrix <- temp
          start <- FALSE
        } else {
          tuning_matrix <- rbind(tuning_matrix, temp)
        }
      }
    }
  }

  string_tuning <- apply(tuning_matrix, 1, paste0, collapse = ":")
  unique_tuning <- unique(string_tuning)
  tuning_matrix <- tuning_matrix[match(unique_tuning, string_tuning), ]

  return(tuning_matrix)
}

#' @importFrom gtools permutations
Tuning_setup_group_func <- function(
    lambda_vec,
    lambda_vec_limit_len,
    # p_threshold,
    num_func,
    tau_vec,
    sub_tuning,
    lim_lambda,
    len_lambda,
    # len_lim_lambda,
    median_val,
    equal_lambda = TRUE,
    Q) {
  if (is.null(lambda_vec)) {
    lambda_vec <- seq(0, lambda_vec_limit_len[1], length.out = lambda_vec_limit_len[2])
  }

  if (equal_lambda) {
    func_lambda_inputs <- matrix(rep(lambda_vec, each = num_func), length(lambda_vec), num_func, byrow = TRUE)
  } else {
    func_lambda0 <- permutations(length(lambda_vec), num_func, repeats.allowed = TRUE)
    func_lambda_inputs <- matrix(lambda_vec[c(func_lambda0)], nrow(func_lambda0), ncol(func_lambda0), byrow = FALSE)
  }

  if (Q == 1) {
    ratios <- 1
  } else {
    ratios <- seq(1, 0, len = sub_tuning)
  }

  lim_lambda_vec <- seq(min(lim_lambda), max(lim_lambda), len = len_lambda)

  start <- TRUE
  for (i in seq_along(tau_vec)) {
    tau <- tau_vec[i]

    for (t in seq_along(lim_lambda_vec)) {
      total <- lim_lambda_vec[t]
      lower <- total * ratios
      upper <- total * (1 - ratios)
      for (idx in seq_along(lower)) {
        temp <- cbind(lower[idx], func_lambda_inputs, upper[idx] * tau, tau)
        if (start) {
          tuning_matrix <- temp
          start <- FALSE
        } else {
          tuning_matrix <- rbind(tuning_matrix, temp)
        }
      }
    }
  }

  func_lambda <- permutations( # nolint: object_usage_linter.
    length(lambda_vec),
    num_func,
    repeats.allowed = TRUE
  ) - 1

  string_tuning <- apply(tuning_matrix, 1, paste0, collapse = ":")
  unique_tuning <- unique(string_tuning)
  tuning_matrix <- tuning_matrix[match(unique_tuning, string_tuning), ]

  output <- list(
    func_lambda = func_lambda,
    lambda_vec = lambda_vec,
    tuning_matrix = tuning_matrix
  )

  return(output)
}

Non_zero <- function(xx) {
  return(length(which(xx != 0)))
}

Clean_results <- function(
    beta_matrix,
    num_iter_vec,
    all_tuning_matrix) {
  # Selects only the unique rows, corresponding to unique tuning parameters
  tuning_vec <- apply(all_tuning_matrix, 1, paste0, collapse = ":")
  unique_vec <- unique(tuning_vec)
  mat <- match(unique_vec, tuning_vec)

  num_iter_vec <- num_iter_vec[mat]
  beta_matrix <- beta_matrix[mat, ]
  all_tuning_matrix <- all_tuning_matrix[mat, ]

  # Reorders the rows by the number of non-zero coefficients
  num_counts <- apply(beta_matrix, 1, Non_zero)
  ord <- order(num_counts)

  num_iter_vec <- num_iter_vec[ord]
  beta_matrix <- beta_matrix[ord, ]
  all_tuning_matrix <- all_tuning_matrix[ord, ]

  output <- list(
    beta_matrix = beta_matrix,
    num_iter_vec = num_iter_vec,
    all_tuning_matrix = all_tuning_matrix
  )
  return(output)
}
