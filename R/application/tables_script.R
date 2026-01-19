### ============================================================================
### PACKAGE LOADING
### ============================================================================

packages <- c("dplyr", "haven", "rdrobust", "RDHonest")

for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

### ============================================================================
### USER CONFIGURATION/CUSTOMISATION
### ============================================================================


## Set working directory -------------------------------------------------------
setwd("~")
readRenviron(".Renviron")
workdir <- Sys.getenv("APPLICATION_WORKDIR")
setwd(workdir)

## Do you want to save the output? ---------------------------------------------
SAVE_OUTPUT_TO_CSV = FALSE
OUTPUT_DESTINATION <- "anglavy99_output.csv" 

## Parameter values for printing results ---------------------------------------
PRINT_OUTPUT = TRUE
DECIMAL_PLACES = 2
WIDTH = DECIMAL_PLACES + 4


### ============================================================================
### ANALYSIS CONFIGURATION
### ============================================================================


## Load dataset and functions --------------------------------------------------
data <- read_dta("./final4.dta")
source("../lambdaFRD.R")
source("./application_utils.R")

## Parameter configurations ----------------------------------------------------
tests <- c("verb", "math")
cutoffs <- c(40, 80, 120)
bandwidths <- c(6, 8, 10, 12, 14, 16, 18)

# Bias-aware Anderson-Rubin test grid parameters
AR_LOWER <- -5
AR_UPPER <- 5
GRID_POINTS <- 250

## Results dataframe initialisation --------------------------------------------
total_combinations <- length(tests) * length(cutoffs) * length(bandwidths)

df_results <- data.frame(
  test = character(total_combinations),
  cutoff = integer(total_combinations),
  bandwidth = integer(total_combinations),
  n = integer(total_combinations),
  iv_coeff = numeric(total_combinations),
  frd_coeff = numeric(total_combinations),
  bc_robust_coeff = numeric(total_combinations),
  lambda_1_coeff = numeric(total_combinations),
  lambda_4_coeff = numeric(total_combinations),
  iv_CI_lower = numeric(total_combinations),
  iv_CI_upper = numeric(total_combinations),
  frd_CI_lower = numeric(total_combinations),
  frd_CI_upper = numeric(total_combinations),
  bc_robust_CI_lower = numeric(total_combinations),
  bc_robust_CI_upper = numeric(total_combinations),
  lambda_1_CI_lower = numeric(total_combinations),
  lambda_1_CI_upper = numeric(total_combinations),
  lambda_4_CI_lower = numeric(total_combinations),
  lambda_4_CI_upper = numeric(total_combinations),
  AR_2_CI_lower = numeric(total_combinations),
  AR_2_CI_upper = numeric(total_combinations)
)

df_idx <- 1

tau_grid <- seq(AR_LOWER, AR_UPPER, length.out = GRID_POINTS)


### ============================================================================
### MAIN ANALYSIS LOOP
### ============================================================================


for (test in tests) {
  
  if (PRINT_OUTPUT) {
    cat("\n\n")
    cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", "\n")
    cat(toupper(test), "SCORES", "\n")
    cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", "\n\n") 
  }
  
  test_column = paste0("avg", test)
  
  for (cutoff in cutoffs) {
    
    if (PRINT_OUTPUT) {
      cat("\n")
      cat("=============================================", "\n")
      cat("CUTOFF:", cutoff, "\n")
      cat("=============================================", "\n\n")      
    }
    
    
    ## =======================================================================
    ## BANDWIDTH CALCULATIONS (REFERENCE ONLY)
    ## =======================================================================
    
    # Mean-squared error optimal bandwidth
    bw_mse <- rdbwselect(y = scale(data[[test_column]])[, 1], x = data$cohsize,
                         c = cutoff, covs = data$tipuach, fuzzy = data$classize,
                         bwselect = "mserd", masspoints = "off")
    mse_optimal_bandwidth <- bw_mse$bws[1, 1]
    
    # Coverage optimal bandwidth
    bw_cov <- rdbwselect(y = scale(data[[test_column]])[, 1], x = data$cohsize,
                         c = cutoff, covs = data$tipuach, fuzzy = data$classize,
                         bwselect = "cerrd", masspoints = "off")
    cov_optimal_bandwidth <- bw_cov$bws[1, 1]
    
    for (bandwidth in bandwidths) {
      
      ## =======================================================================
      ## FILTER DATA
      ## =======================================================================
      
      filtered_data <- data %>% 
        filter(abs(cohsize - cutoff) < bandwidth)
      
      Y <- scale(filtered_data[[test_column]])[,1] 
      D <- as.numeric(filtered_data$classize)
      X <- as.numeric(filtered_data$cohsize)
      W <- as.numeric(filtered_data$tipuach)
      C <- filtered_data$schlcode
      Z <- as.numeric(X >= cutoff)
      n <- nrow(filtered_data)
      
      ## =======================================================================
      ## STANDARD 2SLS ESTIMATOR
      ## =======================================================================
      
      iv_estimate <- lambdaFRD(
        Y = Y, D = D, X = X, x0 = cutoff, exog = W, bandwidth = bandwidth,
        Lambda = TRUE, psi = 0, lambda = NULL, tau_0 = 0, p = 1,
        kernel = "uniform", robust = TRUE, alpha = 0.05
      )
      
      iv_coeff <- iv_estimate$tau_lambda
      iv_CI_lower <- iv_estimate$ci_lower_robust
      iv_CI_upper <- iv_estimate$ci_upper_robust
      
      ## =======================================================================
      ## LAMBDA CLASS ESTIMATOR WITH λ = Λ(1)
      ## =======================================================================
      
      lambda_1_estimate <- lambdaFRD(
        Y = Y, D = D, X = X, x0 = cutoff, exog = W, bandwidth = bandwidth,
        Lambda = TRUE, psi = 1, lambda = NULL, tau_0 = 0, p = 1,
        kernel = "uniform", robust = TRUE, alpha = 0.05
      )
      
      lambda_1_coeff <- lambda_1_estimate$tau_lambda
      lambda_1_CI_lower <- lambda_1_estimate$ci_lower_robust
      lambda_1_CI_upper <- lambda_1_estimate$ci_upper_robust

      ## =======================================================================
      ## LAMBDA CLASS ESTIMATOR WITH λ = Λ(4)
      ## =======================================================================
      
      lambda_4_estimate <- lambdaFRD(
        Y = Y, D = D, X = X, x0 = cutoff, exog = W, bandwidth = bandwidth,
        Lambda = TRUE, psi = 4, lambda = NULL, tau_0 = 0, p = 1,
        kernel = "uniform", robust = TRUE, alpha = 0.05
      )
      
      lambda_4_coeff <- lambda_4_estimate$tau_lambda
      lambda_4_CI_lower <- lambda_4_estimate$ci_lower_robust
      lambda_4_CI_upper <- lambda_4_estimate$ci_upper_robust
      
      ## =======================================================================
      ## STANDARD FRD AND BIAS-CORRECTED ESTIMATORS
      ## =======================================================================
      
      BC_CI_results <- rdrobust(y = Y, x = X, covs = W, c = cutoff, fuzzy = D, 
                                h = bandwidth, kernel = "triangular",
                                masspoints = "off")
      
      frd_coeff <- BC_CI_results$coef[1]
      frd_CI_lower <- BC_CI_results$ci[1,1]
      frd_CI_upper <- BC_CI_results$ci[1,2]
      
      bc_robust_coeff <- BC_CI_results$coef[3]
      bc_robust_CI_lower <- BC_CI_results$ci[3,1]
      bc_robust_CI_upper <- BC_CI_results$ci[3,2]
      
      ## =======================================================================
      ## BIAS AWARE AR CONFIDENCE INTERVALS
      ## =======================================================================
      
      d <- list()
      d$Y <- Y
      d$D <- D
      d$X <- X
      d$ind.X <- (d$X >= cutoff)
      
      M_rot2 <- tryCatch({
        calculate_ROT2(d, cutoff)
      }, error = function(e) {
        cat("ROT2 calculation failed: ", e$message, "\n")
        return(c(NA, NA))
      })
      
      df <- data.frame(Y = d$Y, D = d$D, X = d$X)
      
      # Create grid of values for AR test
      accepted_AR2 <- logical(length(tau_grid))
      
      # Loop over grid
      for (j in seq_along(tau_grid)) {
        tau0 <- tau_grid[j]
        
        df_transformed <- data.frame(
          Y = d$Y - tau0 * d$D,
          X = d$X
        )
        
        AR2 <- suppressMessages(try({
          RDHonest(
            Y ~ X,
            data = df_transformed,
            M = M_rot2[1],
            cutoff = cutoff,
            kern = "triangular",
            h = bandwidth,
            sclass = "H",
            opt.criterion = "FLCI"
          )
        }, silent = TRUE))
        
        if (!inherits(AR2, "try-error")) {
          accepted_AR2[j] <- (0 >= AR2$coef$conf.low) &
            (0 <= AR2$coef$conf.high)
        } else {
          accepted_AR2[j] <- FALSE
        }
      }
      
      # Assess shape of the AR confidence interval
      if (!any(accepted_AR2)) {
        # Checks for empty
        AR2_CI_type  <- "empty"
        AR2_CI_lower <- NA
        AR2_CI_upper <- NA
      } else {
        # Checks for union of disconnected intervals
        idx_true <- which(accepted_AR2)
        gaps     <- diff(idx_true)
        if (any(gaps > 1)) {
          AR2_CI_type  <- "disconnected"
          AR2_CI_lower <- max(tau_grid[idx_true[idx_true < mean(idx_true)]])
          AR2_CI_upper <- min(tau_grid[idx_true[idx_true > mean(idx_true)]])
        } else {
          # Checks for unbounded
          if (min(tau_grid[idx_true]) == AR_LOWER && max(tau_grid[idx_true]) == AR_UPPER) {
            AR2_CI_type <- "unbounded"
            AR2_CI_lower <- -Inf
            AR2_CI_upper <- Inf
          } else {
            # Standard closed interval
            AR2_CI_type  <- "bounded"
            AR2_CI_lower <- min(tau_grid[idx_true])
            AR2_CI_upper <- max(tau_grid[idx_true]) 
          }
        }
      }
      
      # Store output in the results dataframe
      df_results[df_idx, ] <- c(
        test,
        cutoff,
        bandwidth,
        n,
        iv_coeff,
        frd_coeff,
        bc_robust_coeff,
        lambda_1_coeff,
        lambda_4_coeff,
        iv_CI_lower,
        iv_CI_upper,
        frd_CI_lower,
        frd_CI_upper,
        bc_robust_CI_lower,
        bc_robust_CI_upper,
        lambda_1_CI_lower,
        lambda_1_CI_upper,
        lambda_4_CI_lower,
        lambda_4_CI_upper,
        AR_2_CI_lower,
        AR_2_CI_upper
      )
      
      df_idx <- df_idx + 1
      
      ### ======================================================================
      ### PRINT OUTPUT
      ### ======================================================================

      if (PRINT_OUTPUT) {
        # Bandwidth used in estimation -----------------------------------------
        cat("Bandwith:", bandwidth, "N =", n, "\n")
        
        # Reference bandwidths -------------------------------------------------
        cat("---------------------------------------------", "\n")
        cat("MSE optimal bandwidth:",
            fmt_align(mse_optimal_bandwidth, DECIMAL_PLACES, WIDTH), "\n")
        cat("CER optimal bandwidth:",
            fmt_align(cov_optimal_bandwidth, DECIMAL_PLACES, WIDTH), "\n")
        cat("---------------------------------------------", "\n")
        
        # Estimator and confidence interval results ----------------------------
        cat("Estimator       Coeff   Confidence interval\n")
        cat("---------       -----   -------------------\n")      
        cat("IV:           ", fmt_align(iv_coeff, DECIMAL_PLACES, WIDTH), "  [",
            fmt_align(iv_CI_lower, DECIMAL_PLACES, WIDTH), ",",
            fmt_align(iv_CI_upper, DECIMAL_PLACES, WIDTH), "]\n")
        cat("FRD:          ", fmt_align(frd_coeff, DECIMAL_PLACES, WIDTH), "  [",
            fmt_align(frd_CI_lower, DECIMAL_PLACES, WIDTH), ",",
            fmt_align(frd_CI_upper, DECIMAL_PLACES, WIDTH), "]\n")
        cat("BC robust:    ", fmt_align(bc_robust_coeff, DECIMAL_PLACES, WIDTH), "  [",
            fmt_align(bc_robust_CI_lower, DECIMAL_PLACES, WIDTH), ",",
            fmt_align(bc_robust_CI_upper, DECIMAL_PLACES, WIDTH), "]\n")
        cat("Lambda 1:     ", fmt_align(lambda_1_coeff, DECIMAL_PLACES, WIDTH), "  [",
            fmt_align(lambda_1_CI_lower, DECIMAL_PLACES, WIDTH), ",",
            fmt_align(lambda_1_CI_upper, DECIMAL_PLACES, WIDTH), "]\n")
        cat("Lambda 4:     ", fmt_align(lambda_4_coeff, DECIMAL_PLACES, WIDTH), "  [",
            fmt_align(lambda_4_CI_lower, DECIMAL_PLACES, WIDTH), ",",
            fmt_align(lambda_4_CI_upper, DECIMAL_PLACES, WIDTH), "]\n")
        cat("Bias-Aware AR        ",
            fmt_ar_ci(AR2_CI_lower, AR2_CI_upper, AR2_CI_type),
            "\n")
        cat("=============================================", "\n") 
      }

    }
    
  }
  
}


if (SAVE_OUTPUT_TO_CSV) {
  write.csv(df_results, OUTPUT_DESTINATION, row.names = FALSE) 
  if (PRINT_OUTPUT) {
    print(paste0("Output saved to: '", OUTPUT_DESTINATION, "'."))
  }
}

if (PRINT_OUTPUT) {
  cat("\n\n Analysis completed successfully!\n")  
} 

### ============================================================================
### END OF SCRIPT
### ============================================================================