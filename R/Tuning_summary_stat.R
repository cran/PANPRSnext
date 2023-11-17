Tuning_summary_stat <- function(
    beta_vec,
    family = "gaussian",
    penalty = "LOG",
    num_tau = 6,
    num_lambda = 100,
    lambda_min = NULL,
    pfactor = 0.1,
    min_to_max = FALSE
) {
    output <- list()
    output[["lambda"]] <- NULL
    output[["tau"]] <- NULL

    if (family == "gaussian") {
        if (penalty == "LOG") {
            l1_max <- max(abs(beta_vec))
            if (is.null(lambda_min)) {
                min_lambda <- l1_max * pfactor
            } else {
                min_lambda <- lambda_min
            }
            maxTau <- max(abs(beta_vec))
            if (min_to_max) {
                thresholds <- c(exp(seq(log(min_lambda), log(l1_max), len = num_lambda)))
            } else {
                thresholds <- c(exp(seq(log(l1_max), log(min_lambda), len = num_lambda)))
            }
            tauset <- c(exp(seq(log(1e-6), log(maxTau), len = num_tau)))

            tau <- lambda <- c()

            for (t in seq_along()(thresholds)) {
                thres <- thresholds[t]
                slambda <- tauset * thres

                tau <- c(tau, tauset)
                lambda <- c(lambda, slambda)
            }

            if (length(tau) != length(lambda)) {
                stop("error")
            }

            output[["lambda"]] <- lambda
            output[["tau"]] <- tau
        }
    }
    return(output)
}
