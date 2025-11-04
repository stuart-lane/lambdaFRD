### ============================================================================
### SETUP
### ============================================================================

packages <- c("parallel", "doParallel", "foreach", "progress", "MASS",
              "triangle", "expm",  "dplyr", "haven", "rdrobust", "RDHonest")

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
workdir <- Sys.getenv(SIMULATIONS_WORKDIR)
setwd(workdir)

csv_folder <- "output/"

source("../lambda_class_function.R")
source("./simulation_utils.R")

### ============================================================================
### PARAMETER CONFIGURATION
### ============================================================================

result_index <- 1

x0 = 0
nsim = 10000
n = 500

coeffscont_lee <- matrix(c(0.48, 1.27, 7.18, 20.21, 21.54, 7.33), ncol = 1)
coeffstreat_lee <- matrix(c(0.52, 0.84, -3, 7.99, -9.01, 3.56), ncol = 1)
coeffscont_lm <- matrix(c(3.7, 2.99, 3.28, 1.45, 0.22, 0.03), ncol = 1)
coeffstreat_lm <- matrix(c(0.26, 18.49, -54.8, 74.3, -45.02, 9.83), ncol = 1)

### ============================================================================
### DATA FRAME CONFIGURATION
### ============================================================================

dgps = c(1, 2)
running_variables = c(1, 2, 3)
setups = c(1, 2, 3)
p1_values = c(0.6, 0.7, 0.8, 0.9)

total_combinations <- length(running_variables) * length(setups) * length(p1_values) * length(dgps)

results_df <- data.frame(
  
  ## Paramater configurations
  dgp = integer(total_combinations),
  running_variable = integer(total_combinations),
  setup = integer(total_combinations),
  p1 = numeric(total_combinations),
  n = numeric(total_combinations),
  true_treatment_effect = numeric(total_combinations),
  
  ## Bias
  bias_iv_ccf = numeric(total_combinations),
  bias_iv_ik = numeric(total_combinations),
  bias_lambda_1_ccf = numeric(total_combinations),
  bias_lambda_1_ik = numeric(total_combinations),
  bias_lambda_4_ccf = numeric(total_combinations),
  bias_lambda_4_ik = numeric(total_combinations),
  
  ## Median bias
  median_bias_iv_ccf = numeric(total_combinations),
  median_bias_iv_ik = numeric(total_combinations),
  median_bias_lambda_1_ccf = numeric(total_combinations),
  median_bias_lambda_1_ik = numeric(total_combinations),
  median_bias_lambda_4_ccf = numeric(total_combinations),
  median_bias_lambda_4_ik = numeric(total_combinations),
  
  ## Median asbolute deviation
  mad_iv_ccf = numeric(total_combinations),
  mad_iv_ik = numeric(total_combinations),
  mad_lambda_1_ccf = numeric(total_combinations),
  mad_lambda_1_ik = numeric(total_combinations),
  mad_lambda_4_ccf = numeric(total_combinations),
  mad_lambda_4_ik = numeric(total_combinations),
  
  ## Root mean squared error
  rmse_iv_ccf = numeric(total_combinations),
  rmse_iv_ik = numeric(total_combinations),
  rmse_lambda_1_ccf = numeric(total_combinations),
  rmse_lambda_1_ik = numeric(total_combinations),
  rmse_lambda_4_ccf = numeric(total_combinations),
  rmse_lambda_4_ik = numeric(total_combinations),
  
  ## Confidence interval coverage
  cov_iv_ccf = numeric(total_combinations),
  cov_iv_ik = numeric(total_combinations),
  cov_lambda_1_ccf = numeric(total_combinations),
  cov_lambda_1_ik = numeric(total_combinations),
  cov_lambda_4_ccf = numeric(total_combinations),
  cov_lambda_4_ik = numeric(total_combinations),
  cov_ar_ba_1 = numeric(total_combinations),
  cov_ar_ba_2 = numeric(total_combinations),
  cov_honest_1 = numeric(total_combinations),  
  cov_honest_2 = numeric(total_combinations),  
  
  ## Confidence interval length
  len_iv_ccf = numeric(total_combinations),
  len_iv_ik = numeric(total_combinations),
  len_lambda_1_ccf = numeric(total_combinations),
  len_lambda_1_ik = numeric(total_combinations),
  len_lambda_4_ccf = numeric(total_combinations),
  len_lambda_4_ik = numeric(total_combinations),
  len_ar_ba_1 = numeric(total_combinations),
  len_ar_ba_2 = numeric(total_combinations),
  len_honest_1 = numeric(total_combinations),  
  len_honest_2 = numeric(total_combinations)  
)

### ============================================================================
### PARALLELISED SIMULATIONS FUNCTION
### ============================================================================

parallelise_simulation <- function(
    dgps = c(2, 1),
    running_variables = c(1, 2, 3),
    setups = c(1, 2, 3),
    p1_values = c(0.6, 0.7, 0.8, 0.9),
    ns = c(300, 600)
) {
  dir.create(csv_folder, showWarnings = FALSE, recursive = TRUE)
  
  ## Setup parallel processing -------------------------------------------------
  cat("\n----- New Simulation Started -----\n")
  cat("Date and Time:", as.character(round(Sys.time(), units = "secs")), "\n")
  cat("================================== \n")
  total_cores = detectCores()
  num_cores <- max(1, total_cores - 1)
  cat("Using", num_cores, "/", total_cores, "cores\n")
  
  ## Loop over DGP and running variable combinations ---------------------------
  for (dgp in dgps) {
    for (running_variable in running_variables) {
      
      cat("\n")
      cat("===========================================\n")
      cat("Processing DGP =", dgp, ", Running Variable =", running_variable, "\n")
      cat("===========================================\n")
      
      start_time_combo <- Sys.time()
      
      ## Create parameter grid -------------------------------------------------
      combinations <- expand.grid(
        dgp = dgp,                  
        running_variable = running_variable,
        setup = setups,               
        p1 = p1_values,               
        n = ns                        
      )
      
      cat("Parameter combinations for this file:", nrow(combinations), "\n")
      cat("-------------------------------------------\n\n")
      
      ## Create cluster and export ---------------------------------------------
      cl <- makeCluster(num_cores)
      registerDoParallel(cl)
      
      clusterExport(cl, c(
        "generate_rdd_data", "calculate_true_parameters",
        "calculate_ROT1", "calculate_ROT2", "lambda_class",
        "coeffscont_lee", "coeffstreat_lee",
        "coeffscont_lm", "coeffstreat_lm",
        "x0", "nsim", "csv_folder"
      ))
      
      clusterEvalQ(cl, {
        library(rdrobust)
        library(RDHonest)
        library(triangle)
        library(dplyr)
      })
      
      ## Run parallel simulation for this combination --------------------------
      results_list <- foreach(
        idx = 1:nrow(combinations),
        .packages = c("rdrobust", "RDHonest", "triangle", "dplyr"),
        .combine = rbind,
        .errorhandling = "pass"
      ) %dopar% {
        
        # Extract parameters for this iteration
        p1 <- combinations$p1[idx]
        setup <- combinations$setup[idx]
        running_variable <- combinations$running_variable[idx]
        dgp <- combinations$dgp[idx]
        n <- combinations$n[idx]
        
        true_params <- calculate_true_parameters(dgp, setup)
        true_treatment_effect <- true_params$true_treatment_effect
        M_true <- true_params$M_true
        
        true_treatment_effect <- ifelse(dgp == 1, -3.44, 0.04)
        
        running_variable_label <- switch(
          as.character(running_variable),
          "1" = "Normal",
          "2" = "Uniform", 
          "3" = "Beta"
        )
        
        # Initialise vectors to store results
        results_iv_ccf <- numeric(nsim)
        results_iv_ik <- numeric(nsim)
        results_lambda_1_ik <- numeric(nsim)
        results_lambda_1_ccf <- numeric(nsim)
        results_lambda_4_ik <- numeric(nsim)
        results_lambda_4_ccf <- numeric(nsim)
        coverage_iv_ccf <- numeric(nsim)
        coverage_iv_ik <- numeric(nsim)
        coverage_lambda_1_ccf <- numeric(nsim)
        coverage_lambda_1_ik <- numeric(nsim)
        coverage_lambda_4_ik <- numeric(nsim) 
        coverage_lambda_4_ccf <- numeric(nsim)
        coverage_ar_ba_1 <- numeric(nsim)
        coverage_ar_ba_2 <- numeric(nsim)
        coverage_honest_1 <- numeric(nsim)
        coverage_honest_2 <- numeric(nsim)
        length_iv_ccf <- numeric(nsim)
        length_iv_ik <- numeric(nsim)
        length_lambda_1_ccf <- numeric(nsim)
        length_lambda_1_ik <- numeric(nsim)
        length_lambda_4_ik <- numeric(nsim)
        length_lambda_4_ccf <- numeric(nsim)
        length_ar_ba_1 <- numeric(nsim)
        length_ar_ba_2 <- numeric(nsim)
        length_honest_1 <- numeric(nsim)
        length_honest_2 <- numeric(nsim)
        
        iter_seeds <- numeric(nsim)
        
        ## =====================================================================
        ## MAIN SIMULATION LOOP
        ## =====================================================================
        
        for (i in 1:nsim) {
          param_id <- dgp * 10000 + running_variable * 1000 + setup * 100 + 
            round(p1 * 10) + (n %/% 100)
          iter_seed <- 1234 + param_id * nsim + i
          set.seed(iter_seed)
          iter_seeds[i] <- iter_seed
          
          tryCatch({
            sim_data <- generate_rdd_data(n, running_variable, setup, p1, dgp)
            
            df <- data.frame(y = sim_data$y, X = sim_data$x, D = sim_data$D)
            
            sim_data$y <- as.matrix(sim_data$y, ncol = 1)
            sim_data$x <- as.matrix(sim_data$x, ncol = 1)
            sim_data$D <- as.matrix(sim_data$D, ncol = 1)
            
            ## =================================================================
            ## BANDWIDTH SELECTION
            ## =================================================================
            
            ## Coverage optimal bandwidth --------------------------------------
            bw_ccf <- rdbwselect(y = sim_data$y, x = sim_data$x, c = x0,
                                 fuzzy = sim_data$D, bwselect = "cerrd")
            ccf_bw <- bw_ccf$bws[1, 1]
            
            ## MSE optimal bandwidth -------------------------------------------
            bw_ik <- rdbwselect(y = sim_data$y, x = sim_data$x, c = x0,
                                fuzzy = sim_data$D)
            ik_bw <- bw_ik$bws[1, 1]
            
            ## =================================================================
            ## LAMBDA CLASS ESTIMATOR WITH λ = Λ(1)
            ## =================================================================
            
            ## Lambda 1 estimator with mse optimal bandwidth -------------------
            result_lambda_1_ik <- lambda_class(
              Y = sim_data$y, D = sim_data$D, X = sim_data$x, x0 = x0, exog = NULL,
              bandwidth = ik_bw, Lambda = TRUE, psi = 1,  lambda = NULL, 
              tau_0 = true_treatment_effect, p = 1, kernel = "triangular", cluster = NULL,
              robust = FALSE, alpha = 0.05
            )
            
            results_lambda_1_ik[i] <- result_lambda_1_ik$tau_est
            coverage_lambda_1_ik[i] <- as.numeric(1 - result_lambda_1_ik$reject_t)
            length_lambda_1_ik[i] <- as.numeric(result_lambda_1_ik$ci_upper -
                                                  result_lambda_1_ik$ci_lower)
            
            ## Lambda 1 estimator with coverage optimal bandwidth --------------
            result_lambda_1_ccf <- lambda_class(
              Y = sim_data$y, D = sim_data$D, X = sim_data$x, x0 = 0, exog = NULL,
              bandwidth = ccf_bw, Lambda = TRUE, psi = 1, lambda = NULL, 
              tau_0 = true_treatment_effect, p = 1, kernel = "triangular", cluster = NULL,
              robust = FALSE, alpha = 0.05
            )
            
            results_lambda_1_ccf[i] <- result_lambda_1_ccf$tau_est
            coverage_lambda_1_ccf[i] <- as.numeric(1 - result_lambda_1_ccf$reject_t)
            length_lambda_1_ccf[i] <-  as.numeric(result_lambda_1_ccf$ci_upper -
                                                    result_lambda_1_ccf$ci_lower) 
            
            ## =================================================================
            ## LAMBDA CLASS ESTIMATOR WITH λ = Λ(4)
            ## =================================================================
            
            ## Lambda 4 estimator with mse optimal bandwidth -------------------
            result_lambda_4_ik <- lambda_class(
              Y = sim_data$y, D = sim_data$D, X = sim_data$x, x0 = 0, exog = NULL,
              bandwidth = ik_bw, Lambda = TRUE, psi = 4, lambda = NULL, 
              tau_0 = true_treatment_effect, p = 1, kernel = "triangular", cluster = NULL,
              robust = FALSE, alpha = 0.05
            )
            
            results_lambda_4_ik[i] <- result_lambda_4_ik$tau_est
            coverage_lambda_4_ik[i] <- as.numeric(1 -result_lambda_4_ik$reject_t)
            length_lambda_4_ik[i] <-  as.numeric(result_lambda_4_ik$ci_upper -
                                                   result_lambda_4_ik$ci_lower) 
            
            ## Lambda 1 estimator with cov optimal bandwidth -------------------
            result_lambda_4_ccf <- lambda_class(
              Y = sim_data$y, D = sim_data$D, X = sim_data$x, x0 = 0, exog = NULL,
              bandwidth = ccf_bw, Lambda = TRUE, psi = 4,  lambda = NULL, 
              tau_0 = true_treatment_effect, p = 1, kernel = "triangular", cluster = NULL,
              robust = FALSE, alpha = 0.05
            )
            
            results_lambda_4_ccf[i] <- result_lambda_4_ccf$tau_est
            coverage_lambda_4_ccf[i] <- as.numeric(1 - result_lambda_4_ccf$reject_t)
            length_lambda_4_ccf[i] <-  as.numeric(result_lambda_4_ccf$ci_upper -
                                                    result_lambda_4_ccf$ci_lower) 
            
            ## =================================================================
            ## RD ROBUST CONFIDENCE INTERVALS
            ## =================================================================

            ## FRD estimator with MSE optimal bandwidth ------------------------
            rd_result_ccf <- rdrobust(y = sim_data$y, x = sim_data$x, c = x0,
                                      fuzzy = sim_data$D, h = ccf_bw)
            results_iv_ccf[i] <- rd_result_ccf$coef[1]

            ## Coverage for FRD estimator with MSE optimal bandwidth -----------
            ci_lower_iv_ccf <- rd_result_ccf$ci[3,1]
            ci_upper_iv_ccf <- rd_result_ccf$ci[3,2]
            coverage_iv_ccf[i] <- as.numeric(true_treatment_effect >= ci_lower_iv_ccf &
                                               true_treatment_effect <= ci_upper_iv_ccf)
            length_iv_ccf[i] <- as.numeric(ci_upper_iv_ccf - ci_lower_iv_ccf)

            ## FRD estimator with Imbens-Kalyanaraman bandwidth ----------------
            rd_result_ik <- rdrobust(y = sim_data$y, x = sim_data$x,
                                     c = x0, fuzzy = sim_data$D, h = ik_bw)
            results_iv_ik[i] <- rd_result_ik$coef[1]

            ## Coverage for FRD estimator with Imbens-Kalyanaraman bandwidth -
            ci_lower_iv_ik <- rd_result_ik$ci[3,1]
            ci_upper_iv_ik <- rd_result_ik$ci[3,2]
            coverage_iv_ik[i] <- as.numeric(true_treatment_effect >= ci_lower_iv_ik &
                                              true_treatment_effect <= ci_upper_iv_ik)
            length_iv_ik[i] <- as.numeric(ci_upper_iv_ik - ci_lower_iv_ik)

            ## =================================================================
            ## BIAS AWARE CONFIDENCE INTERVALS
            ## =================================================================

            d <- list()
            d$Y <- sim_data$y
            d$D <- sim_data$D
            d$X <- sim_data$x
            d$ind.X <- (d$X >= 0)

            M_rot1 <- calculate_ROT1(d, x0)
            M_rot2 <- calculate_ROT2(d, x0)

            df_original <- data.frame(Y = d$Y, D = d$D, X = d$X)

            df_transformed <- data.frame(
              Y = d$Y - true_treatment_effect * d$D,
              X = d$X
            )

            ## AR1 -------------------------------------------------------------
            reg_ar_1 <- suppressMessages(try({
              RDHonest(Y ~ X,
                       data = df_transformed,
                       M = M_rot1[1] + abs(true_treatment_effect) * M_rot1[2],
                       cutoff = 0,
                       kern = "triangular",
                       sclass = "H",
                       opt.criterion = "FLCI")
            }, silent = TRUE))

            if (!inherits(reg_ar_1, "try-error")) {
              ar_ccf_coverage <- as.numeric((0 >= reg_ar_1$coef$conf.low) &
                                              (0 <= reg_ar_1$coef$conf.high))
              ar_ccf_bandwidth <- reg_ar_1$coef$bandwidth
            } else {
              ar_ccf_coverage <- NA
              ar_ccf_bandwidth <- NA
            }

            coverage_ar_ba_1[i] <- ar_ccf_coverage
            length_ar_ba_1[i] <- as.numeric(reg_ar_1$coef$conf.high - reg_ar_1$coef$conf.low)

            ## AR2 -------------------------------------------------------------
            reg_ar_2 <- suppressMessages(try({
              RDHonest(Y ~ X,
                       data = df_transformed,
                       M = M_rot2[1] + abs(true_treatment_effect) * M_rot2[2],
                       cutoff = 0,
                       kern = "triangular",
                       sclass = "H",
                       opt.criterion = "FLCI")
            }, silent = TRUE))

            if (!inherits(reg_ar_2, "try-error")) {
              ar_2_coverage <- as.numeric((0 >= reg_ar_2$coef$conf.low) &
                                            (0 <= reg_ar_2$coef$conf.high))
              ar_2_bandwidth <- reg_ar_2$coef$bandwidth
            } else {
              ar_2_coverage <- NA
              ar_2_bandwidth <- NA
            }

            coverage_ar_ba_2[i] <- ar_2_coverage
            length_ar_ba_2[i] <- as.numeric(reg_ar_2$coef$conf.high - reg_ar_2$coef$conf.low)

            ## =================================================================
            ## HONEST CONFIDENCE INTERVALS
            ## =================================================================

            ## Honest confidence interval 1 ------------------------------------
            frd_result <- suppressMessages(try({RDHonest(
              y | D ~ X,
              data = df,
              kern = "triangular",
              M = abs(M_rot1),
              sclass = "H",
              opt.criterion = "FLCI")
            }, silent = TRUE))

            summary <- frd_result$coefficients
            ci_lower_honest_1 <- summary$conf.low
            ci_upper_honest_1 <- summary$conf.high
            coverage_honest_1[i] <- as.numeric(true_treatment_effect >= ci_lower_honest_1 &
                                                 true_treatment_effect <= ci_upper_honest_1)
            length_honest_1[i] <- as.numeric(ci_upper_honest_1 - ci_lower_honest_1)

            ## Honest confidence interval 2 ------------------------------------
            frd_result <- suppressMessages(try({RDHonest(
              y | D ~ X,
              data = df,
              kern = "triangular",
              M = abs(M_rot2),
              sclass = "H",
              opt.criterion = "FLCI")
            }, silent = TRUE))

            summary <- frd_result$coefficients
            ci_lower_honest_2 <- summary$conf.low
            ci_upper_honest_2 <- summary$conf.high
            coverage_honest_2[i] <- as.numeric(true_treatment_effect >= ci_lower_honest_2 &
                                                 true_treatment_effect <= ci_upper_honest_2)
            length_honest_2[i] <- as.numeric(ci_upper_honest_2 - ci_lower_honest_2)
            
          }, error = function(e) {
            warning(paste("Error in simulation iteration:", e$message))
          })
        }
        
        ## =====================================================================
        ## COMPUTE SUMMARY STATISTICS FOR PARAMETER CONFIGURATION 
        ## =====================================================================
        
        # Calculate Bias
        bias_iv_ik <- mean(results_iv_ik - true_treatment_effect)
        bias_iv_ccf <- mean(results_iv_ccf - true_treatment_effect)
        bias_lambda_1_ik <- mean(results_lambda_1_ik - true_treatment_effect)
        bias_lambda_1_ccf <- mean(results_lambda_1_ccf - true_treatment_effect)
        bias_lambda_4_ik <- mean(results_lambda_4_ik - true_treatment_effect)
        bias_lambda_4_ccf <- mean(results_lambda_4_ccf - true_treatment_effect)
        
        median_bias_iv_ik <- median(results_iv_ik - true_treatment_effect)
        median_bias_iv_ccf <- median(results_iv_ccf - true_treatment_effect)
        median_bias_lambda_1_ik <- median(results_lambda_1_ik - true_treatment_effect)
        median_bias_lambda_1_ccf <- median(results_lambda_1_ccf - true_treatment_effect)
        median_bias_lambda_4_ik <- median(results_lambda_4_ik - true_treatment_effect)
        median_bias_lambda_4_ccf <- median(results_lambda_4_ccf - true_treatment_effect)
        
        # Calculate Median Absolute Deviation
        mad_iv_ik <- median(abs(results_iv_ik - true_treatment_effect))
        mad_iv_ccf <- median(abs(results_iv_ccf - true_treatment_effect))
        mad_lambda_1_ik <- median(abs(results_lambda_1_ik - true_treatment_effect))
        mad_lambda_1_ccf <- median(abs(results_lambda_1_ccf - true_treatment_effect))
        mad_lambda_4_ik <- median(abs(results_lambda_4_ik - true_treatment_effect))
        mad_lambda_4_ccf <- median(abs(results_lambda_4_ccf - true_treatment_effect))
        
        # Calculate Mean Squared Errors
        rmse_iv_ik <- sqrt(mean((results_iv_ik - true_treatment_effect)^2))
        rmse_iv_ccf <- sqrt(mean((results_iv_ccf - true_treatment_effect)^2))
        rmse_lambda_1_ik <- sqrt(mean((results_lambda_1_ik - true_treatment_effect)^2))
        rmse_lambda_1_ccf <- sqrt(mean((results_lambda_1_ccf - true_treatment_effect)^2))
        rmse_lambda_4_ik <- sqrt(mean((results_lambda_4_ik - true_treatment_effect)^2))
        rmse_lambda_4_ccf <- sqrt(mean((results_lambda_4_ccf - true_treatment_effect)^2))
        
        # Calculate Coverage Probabilities
        cov_iv_ik <- mean(coverage_iv_ik)
        cov_iv_ccf <- mean(coverage_iv_ccf)
        cov_lambda_1_ik <- mean(coverage_lambda_1_ik)
        cov_lambda_1_ccf <- mean(coverage_lambda_1_ccf)
        cov_lambda_4_ik <- mean(coverage_lambda_4_ik)
        cov_lambda_4_ccf <- mean(coverage_lambda_4_ccf)
        cov_ar_ba_1 <- mean(coverage_ar_ba_1)
        cov_ar_ba_2 <- mean(coverage_ar_ba_2)
        cov_honest_1 <- mean(coverage_honest_1)
        cov_honest_2 <- mean(coverage_honest_2)
        
        # Calculate Average CI Lengths
        len_iv_ik <- mean(length_iv_ik)
        len_iv_ccf <- mean(length_iv_ccf)
        len_lambda_1_ik <- mean(length_lambda_1_ik)
        len_lambda_1_ccf <- mean(length_lambda_1_ccf)
        len_lambda_4_ik <- mean(length_lambda_4_ik)
        len_lambda_4_ccf <- mean(length_lambda_4_ccf)
        len_ar_ba_1 <- mean(length_ar_ba_1)
        len_ar_ba_2 <- mean(length_ar_ba_2)
        len_honest_1 <- mean(length_honest_1)
        len_honest_2 <- mean(length_honest_2)
        
        # Add metrics to the results dataframe
        data.frame(
          # Store parameter configuration
          dgp = dgp,
          running_variable = running_variable,
          setup = setup,
          p1 = p1,
          n = n,
          true_treatment_effect = true_treatment_effect,
          
          # Store bias
          bias_iv_ccf = bias_iv_ccf,
          bias_iv_ik = bias_iv_ik,
          bias_lambda_1_ccf = bias_lambda_1_ccf,
          bias_lambda_1_ik = bias_lambda_1_ik,
          bias_lambda_4_ccf = bias_lambda_4_ccf,
          bias_lambda_4_ik = bias_lambda_4_ik,
          
          # Store median bias
          median_bias_iv_ccf = median_bias_iv_ccf,
          median_bias_iv_ik = median_bias_iv_ik,
          median_bias_lambda_1_ccf = median_bias_lambda_1_ccf,
          median_bias_lambda_1_ik = median_bias_lambda_1_ik,
          median_bias_lambda_4_ccf = median_bias_lambda_4_ccf,
          median_bias_lambda_4_ik = median_bias_lambda_4_ik,
          
          # Store median absolute deviation
          mad_iv_ccf = mad_iv_ccf,
          mad_iv_ik = mad_iv_ik,
          mad_lambda_1_ccf = mad_lambda_1_ccf,
          mad_lambda_1_ik = mad_lambda_1_ik,
          mad_lambda_4_ccf = mad_lambda_4_ccf,
          mad_lambda_4_ik = mad_lambda_4_ik,
          
          # Store root mean squared error
          rmse_iv_ccf = rmse_iv_ccf,
          rmse_iv_ik = rmse_iv_ik,
          rmse_lambda_1_ccf = rmse_lambda_1_ccf,
          rmse_lambda_1_ik = rmse_lambda_1_ik,
          rmse_lambda_4_ccf = rmse_lambda_4_ccf,
          rmse_lambda_4_ik = rmse_lambda_4_ik,
          
          # Store empirical coverage
          cov_iv_ccf = cov_iv_ccf,
          cov_iv_ik = cov_iv_ik,
          cov_lambda_1_ccf = cov_lambda_1_ccf,
          cov_lambda_1_ik = cov_lambda_1_ik,
          cov_lambda_4_ccf = cov_lambda_4_ccf,
          cov_lambda_4_ik = cov_lambda_4_ik,
          cov_ar_ba_1 = cov_ar_ba_1,
          cov_ar_ba_2 = cov_ar_ba_2,
          cov_honest_1 = cov_honest_1,
          cov_honest_2 = cov_honest_2,
          
          # Store empirical length
          len_iv_ccf = len_iv_ccf,
          len_iv_ik = len_iv_ik,
          len_lambda_1_ccf = len_lambda_1_ccf,
          len_lambda_1_ik = len_lambda_1_ik,
          len_lambda_4_ccf = len_lambda_4_ccf,
          len_lambda_4_ik = len_lambda_4_ik,
          len_ar_ba_1 = len_ar_ba_1,
          len_ar_ba_2 = len_ar_ba_2,
          len_honest_1 = len_honest_1,
          len_honest_2 = len_honest_2
        )
      }
      
      ## Stop the cluster ------------------------------------------------------
      stopCluster(cl)
      
      ## Save results for this DGP × running variable combination --------------
      final_results_df <- as.data.frame(results_list)
      
      # Create descriptive filename
      running_var_label <- switch(
        as.character(running_variable),
        "1" = "Normal",
        "2" = "Uniform",
        "3" = "Beta"
      )
      
      filename <- paste0("simulation_results_dgp", dgp, 
                         "_rv", running_var_label, ".csv")
      filepath <- file.path(csv_folder, filename)
      
      write.csv(final_results_df, filepath, row.names = FALSE)
      
      end_time_combo <- Sys.time()
      time_taken <- difftime(end_time_combo, start_time_combo, units = "mins")
      
      cat("\n✓ Completed DGP =", dgp, ", Running Variable =", running_variable, "\n")
      cat("  Time taken:", round(time_taken, 2), "minutes\n")
      cat("  Results saved to:", filename, "\n\n")
    }
  }
  
  cat("\n")
  cat("========================================\n")
  cat("ALL SIMULATIONS COMPLETE!\n")
  cat("========================================\n")
}

## Run the simulation ----------------------------------------------------------
start_time_overall <- Sys.time()
parallelise_simulation()
end_time_overall <- Sys.time()

cat("\nTotal simulation time:", difftime(end_time_overall, start_time_overall, units = "hours"), "hours\n\n")