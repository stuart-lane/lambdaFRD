#' Lambda-Class Estimator for Fuzzy Regression Discontinuity Designs
#'
#' Implements the lambda-class of generalised estimators for fuzzy RD designs that
#' have finite moments of all orders in finite samples. The standard FRD estimator 
#' (lambda = 1) lacks finite moments, leading to instability. This function 
#' implements estimators (lambda < 1) that significantly improve finite-sample performance.
#' This function also includes the standard sharp RD estimator (lambda = 0).
#'
#' @param Y Numeric vector of the outcome
#' @param D Numeric vector of the treatment
#' @param X Numeric vector of the running variable
#' @param x0 Cutoff value (default = 0)
#' @param exog Optional matrix of exogenous covariates to include
#' @param bandwidth Bandwidth for local polynomial estimation
#' @param Lambda Logical. If TRUE, compute lambda using Lambda(psi) function. 
#'   If FALSE, use value provided in lambda argument (default = TRUE)
#' @param psi Tuning parameter for Lambda function. Recommended values: 
#'   psi = 4 for best performance in terms of median bias, median absolute deviation,
#'   root mean squared error, and coverage of confidence intervals. psi = 1 for conservative choice.
#'   psi = 0 gives the standard FRD estimator, which does not have finite sample moments. 
#'   (default = 4)
#' @param lambda Manual specification of lambda in [0,1]. Ignored if Lambda = TRUE.
#'   lambda = 1 gives standard FRD estimator (not recommended) In general, this parameter
#'   should be ignored, and the (Lambda = TRUE, psi = 4) combination should be used instead.
#' @param tau_0 Null hypothesis value for testing (default = 0)
#' @param p Order of local polynomial regression (default = 1, i.e., local linear)
#' @param kernel Kernel function: "uniform", "triangular", or "epanechnikov" 
#'   (default = "uniform")
#' @param robust Logical. If TRUE, compute heteroskedasticity-robust standard 
#'   errors (default = FALSE)
#' @param alpha Significance level for confidence intervals (default = 0.05)
#'
#' @return A list containing:
#' \describe{
#'   \item{tau_lambda}{Point estimate of treatment effect using lambda-class estimator}
#'   \item{tau_Y_lambda}{Numerator of the ratio estimator}
#'   \item{tau_D_lambda}{Denominator of the ratio estimator}
#'   \item{coefficients}{Full vector of estimated coefficients using lambda-class estimator}
#'   \item{t_stat}{t-statistic for testing H0: tau = tau_0}
#'   \item{var_est}{Variance estimate (homoskedastic)}
#'   \item{std_error}{Standard error (homoskedastic)}
#'   \item{ci_lower}{Lower bound of confidence interval}
#'   \item{ci_upper}{Upper bound of confidence interval}
#'   \item{reject_t}{Indicator for rejecting H0 at level alpha}
#'   \item{var_robust}{Variance estimate (heteroskedasticity-robust), if robust = TRUE}
#'   \item{std_robust}{Standard error (robust), if robust = TRUE}
#'   \item{t_robust}{t-statistic (robust), if robust = TRUE}
#'   \item{ci_lower_robust}{Lower bound of CI (robust), if robust = TRUE}
#'   \item{ci_upper_robust}{Upper bound of CI (robust), if robust = TRUE}
#'   \item{df}{Degrees of freedom}
#'   \item{M}{Effective sample size within bandwidth}
#' }
#'
#' @details
#' The lambda-class estimator takes the form:
#' \deqn{\hat{\tau}_{\lambda,p} = \frac{\tilde{D}' (I - \lambda M_{\tilde{Z}}) \tilde{Y}}
#'                                      {\tilde{D}' (I - \lambda M_{\tilde{Z}}) \tilde{D}}}
#' where tildes denote kernel-weighted vectors and matrices.
#' 
#' When lambda = 1, this reduces to the standard FRD estimator. For lambda < 1,
#' the estimator has finite moments of all orders and improved finite-sample
#' properties, particularly with small samples or a smaller discontinuity in
#' the probability of assignment to treatment at the cutoff.
#' 
#' The Lambda(psi) function is defined as:
#' \deqn{\Lambda(\psi) = 1 - \frac{\psi}{M - 2(p+1)}}
#' 
#' Recommended values based on extensive simulations:
#' \itemize{
#'   \item psi = 4: Best performance in simulations for median bias, median absolute deviation,
#'                  root mean squared error and confidence interval coverage
#'   \item psi = 1: More conservative, still substantial improvement over standard FRD estimator
#' }
#'
#' @references
#' Lane, S. (2025). "The moment is here: a generalised class of estimators for 
#' fuzzy regression discontinuity designs." Working paper.
#' 
#'
#' @examples
#' # Simulate data
#' set.seed(123)
#' n <- 500
#' X <- rnorm(n)
#' P <- ifelse(X < 0, p0, p1)
#' D <- rbinom(n, 1, P)
#' Y <- 0.5 * D + X + X^2 + rnorm(n, 0, 0.3)
#' 
#' # Estimate with recommended settings (psi = 4)
#' result <- lambda_class(Y, D, X, x0 = 0, bandwidth = 1, psi = 4)
#' print(result$tau_lambda)
#' print(c(result$ci_lower, result$ci_upper))
#' 
#' # Compare with standard estimator (lambda = 1, not recommended)
#' result_std <- lambda_class(Y, D, X, x0 = 0, bandwidth = 1, 
#'                           Lambda = FALSE, lambda = 1)
#' 
#' # With heteroskedasticity-robust standard errors
#' result_robust <- lambda_class(Y, D, X, x0 = 0, bandwidth = 1, 
#'                              psi = 4, robust = TRUE)
#' print(c(result_robust$ci_lower_robust, result_robust$ci_upper_robust))
#' 
#' NOTE ON BANDWIDTH SELECTION: ================================================
#' For bandwidth, I recommend using the coverage-optimal bandwidth from the 
#' `rdrobust` package, which with the data above would be
#' 
#' bw_rd <- rdbwselect(y = Y, x = X, c = 0, fuzzy = D, bwselect = "cerrd")
#' h_opt <- bw_rd$bws[1, 1]
#' 
#' The new estimator could then be computed as
#' #' result <- lambda_class(Y, D, X, x0 = 0, bandwidth = h_opt, psi = 4)
#'
#' @export

lambda_class <- function(
    Y, D, X, x0 = 0, exog = NULL, bandwidth, Lambda = TRUE, psi = 4,
    lambda = NULL, tau_0 = NULL, p = 1, kernel = "uniform", robust = FALSE,
    alpha = 0.05
) {
  
  ## ERROR MESSAGES AND WARNINGS ===============================================
  
  if (is.null(Y) || is.null(D) || is.null(X)) {
    stop("Please supply outcome Y, treatment D and running variable X.")
  }
  
  if (length(Y) != length(D) || length(Y) != length(X)) {
    stop("Y, D and X must be the same length")
  }
  
  if (is.null(bandwidth)) {
    stop("Please supply a bandwidth h.")
  }
  
  if (missing(x0)) {
    warning("Default value of x0 is 0. Please supply another value if the cut-off is different.")
  }
  
  if (p %% 1 != 0 || p < 0) {
    stop("Local polynomial order p must be a non-negative integer.")
  }
  
  if (!(alpha > 0 & alpha < 1)) {
    stop("Alpha must be strictly between 0 and 1.")
  }
  
  ## BASIC SETUP ===============================================================
  
  # Ensure inputs are numeric vectors
  Y <- as.numeric(Y)
  D <- as.numeric(D)
  X <- as.numeric(X)
  Z <- as.numeric(X >= x0)
  n <- length(Y)
  
  # Compute kernel weights
  if (kernel == "uniform") {
    weights <- 0.5 * diag(as.numeric(abs(X - x0) <= bandwidth))
  } else if (kernel == "triangular") {
    kernel_argument = 1 - abs(X - x0)/bandwidth
    weights <- diag(sqrt(pmax(0, kernel_argument))) 
  } else if (kernel == "epanechnikov") {
    kernel_argument = 0.75 * (1 - (abs(X - x0)/bandwidth)^2)
    weights <- diag(sqrt(pmax(0, kernel_argument))) 
  } else {
    stop("Kernel argument must take value 'uniform', 'triangular' or 'epanechnikov'.")
  }
  
  # Construct V matrix
  V <- matrix(1, nrow = n, ncol = 1)
  if (p >= 1) {
    for (j in 1:p) {
      V <- cbind(V, (1 - Z) * ((X - x0)^j), Z * ((X - x0)^j))
    } 
  }
  
  # Add exogenous regressors if provided
  if (!is.null(exog)) {
    exog <- as.matrix(exog)
    V <- cbind(V, exog)
  }
  
  # Apply triangular kernel weights
  Vw <- weights %*% V
  Zw <- weights %*% Z
  Yw <- weights %*% Y
  Dw <- weights %*% D
  
  # # Remove rows with zeros
  rowswithzeros <- (Vw[,1] == 0)
  keep_indices <- which(!rowswithzeros)
  
  Vw <- Vw[keep_indices, ]
  Yw <- Yw[keep_indices]
  Dw <- Dw[keep_indices]
  Zw <- Zw[keep_indices]
  
  # Construct Lambda function
  M <- nrow(Vw)
  
  if (Lambda) {
    lambda <- 1 - psi / (M - ncol(V) - 1)
  }
  
  ## CONSTRUCT ESTIMATOR =======================================================
  
  # Bind matrix terms for estimator computation
  VDw <- cbind(Vw, Dw)
  VZw <- cbind(Vw, Zw)
  izz <- solve(t(VZw) %*% VZw)
  PVZw <- VZw %*% izz %*% t(VZw)
  MVZw <- diag(M) - PVZw
  inner = diag(M) - lambda * MVZw
  
  coeffs_num <- t(VDw) %*% inner %*% Yw
  coeffs_den <- t(VDw) %*% inner %*% VDw
  coefficients <- solve(coeffs_den) %*% coeffs_num
  
  # =======================================================
  
  MV <- diag(M) - Vw %*% solve(t(Vw) %*% Vw) %*% t(Vw)
  Y_tilde <- MV %*% Yw
  D_tilde <- MV %*% Dw
  Z_tilde <- MV %*% Zw
  
  izz_tilde <- solve(t(Z_tilde) %*% Z_tilde)
  Pz_tilde <- Z_tilde %*% izz_tilde %*% t(Z_tilde)
  Mz_tilde <- diag(M) - Pz_tilde
  
  tau_Y_lambda = t(D_tilde) %*% (diag(M) - lambda*Mz_tilde) %*% Y_tilde
  tau_D_lambda = t(D_tilde) %*% (diag(M) - lambda*Mz_tilde) %*% D_tilde
  
  # Calculate estimator
  tau_lambda <- tau_Y_lambda / tau_D_lambda
  
  ## INFERENCE =================================================================
  
  # Calculate Fuller t-statistic
  u_hat <- Y_tilde - D_tilde %*% tau_lambda
  
  sigma_2_hat <- t(u_hat) %*% u_hat
  middle_term <- t(D_tilde) %*% Pz_tilde %*% D_tilde
  denom_term <- tau_D_lambda^2
  
  var_est <- as.numeric((sigma_2_hat * middle_term) / denom_term)
  
  var_robust <- NULL
  if (robust) {
    sigma_2_hat_robust <- t(D_tilde) %*% Pz_tilde %*% diag(as.vector(u_hat)^2) %*% Pz_tilde %*% D_tilde
    var_robust <- solve(tau_D_lambda) %*% sigma_2_hat_robust %*% solve(tau_D_lambda)
  }
  
  # Standard t-statistic
  t_stat <- NULL
  ci_lower <- NULL
  ci_upper <- NULL
  
  if (is.null(tau_0)) {
    tau_0 <- 0  # Default to testing against zero if not provided
  }
  
  # Calculate t-statistic and confidence interval
  df <- M - ncol(V) - 1 
  t_critical <- qt(1 - alpha/2, df)
  
  t_stat <- sqrt(df) * (tau_lambda - tau_0) / sqrt(var_est)
  
  # Confidence interval using homoskedastic standard errors
  std_error <- sqrt(var_est/df)
  ci_lower <- as.numeric(tau_lambda - t_critical * std_error)
  ci_upper <- as.numeric(tau_lambda + t_critical * std_error)
  reject_t <- as.numeric(t_stat > t_critical)
  
  # Heteroskedasticity-robust t-statistic and confidence interval
  t_robust <- NULL
  ci_lower_robust <- NULL
  ci_upper_robust <- NULL
  
  if (!is.null(var_robust)) {
    t_robust <- (tau_lambda - tau_0) / sqrt(var_robust)
    
    # Confidence interval using robust standard errors
    std_robust <- sqrt(var_robust)
    ci_lower_robust <- as.numeric(tau_lambda - t_critical * std_robust)
    ci_upper_robust <- as.numeric(tau_lambda + t_critical * std_robust)
  }
  
  ## COLLECT OUTPUTS ===========================================================
  
  # Return both the estimate and test statistics with clustering
  return(list(
    tau_lambda = as.numeric(tau_lambda),
    tau_Y_lambda = as.numeric(tau_Y_lambda),
    tau_D_lambda = as.numeric(tau_D_lambda),
    coefficients = as.numeric(coefficients),
    t_stat = as.numeric(t_stat),
    var_est = as.numeric(var_est),
    std_error = std_error,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    reject_t = reject_t,
    var_robust = var_robust,
    std_robust = if(!is.null(var_robust)) sqrt(var_robust) else NULL,
    t_robust = t_robust,
    ci_lower_robust = ci_lower_robust,
    ci_upper_robust = ci_upper_robust,
    df = df,
    M = as.numeric(M)
  ))
}