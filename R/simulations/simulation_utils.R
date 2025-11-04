generate_rdd_data <- function(n, running_variable, setup, p1, dgp) {
  p0 <- 1 - p1
  
  # Generate the running variable x based on the specified type
  if (running_variable == 1) {
    x <- rnorm(n, mean = 0, sd = 1)
  } else if (running_variable == 2) {
    x <- runif(n, min = -1, max = 1)
  } else {
    x <- 2 * rbeta(n, shape1 = 2, shape2 = 4) - 1
  }
  
  # Generate error terms
  u <- rnorm(n, mean = 0, sd = 0.3)
  
  # Generate treatments
  if (setup == 1) {
    P <- ifelse(x < 0, p0, p1)
    T <- rbinom(n, 1, P)
  } else if (setup == 2) {
    P <- ifelse(x < -1, 0, ifelse(x <= 0, (1 - p1) * x + (1 - p1),
                                  ifelse(x < 1, (1 - p1) * x + p1, 1)))
    T <- rbinom(n, 1, P)
  } else {
    P <- (x < 0) * p0 * exp(0.2*x) + (x >= 0) * (p1 + p0 * (1 - exp(-0.2*x)))
    T <- rbinom(n, 1, P)
  }
  
  # Generate structural function
  if (dgp == 1) {
    Xmlee <- cbind(1, x, x^2, x^3, x^4, x^5)
    y <- (Xmlee %*% coeffscont_lm) * (x < 0) + (Xmlee %*% coeffstreat_lm) * (x > 0) + u 
  } else {
    Xmlee <- cbind(1, x, x^2, x^3, x^4, x^5)
    y <- (Xmlee %*% coeffscont_lee) * (x < 0) + (Xmlee %*% coeffstreat_lee) * (x > 0) + u 
  }
  
  return(list(y = y, D = T, x = x, exog = NULL))
}

calculate_ROT1 <- function(d, cutoff) {
  # Fit fourth-order polynomials on either side of the cutoff
  fit_y_below <- lm(d$Y ~ d$X + I(d$X^2) + I(d$X^3) + I(d$X^4), subset = which(d$X < cutoff))
  fit_y_above <- lm(d$Y ~ d$X + I(d$X^2) + I(d$X^3) + I(d$X^4), subset = which(d$X >= cutoff))
  fit_d_below <- lm(d$D ~ d$X + I(d$X^2) + I(d$X^3) + I(d$X^4), subset = which(d$X < cutoff))
  fit_d_above <- lm(d$D ~ d$X + I(d$X^2) + I(d$X^3) + I(d$X^4), subset = which(d$X >= cutoff))
  
  # Create a grid of x values to evaluate the second derivatives
  x_grid_below <- seq(min(d$X[d$X < cutoff]), 0, length.out = 100)
  x_grid_above <- seq(0, max(d$X[d$X >= cutoff]), length.out = 100)
  
  # Calculate second derivatives for outcome function
  y_second_deriv_below <- 2*coef(fit_y_below)[3] +
    6*coef(fit_y_below)[4]*x_grid_below +
    12*coef(fit_y_below)[5]*x_grid_below^2
  
  y_second_deriv_above <- 2*coef(fit_y_above)[3] +
    6*coef(fit_y_above)[4]*x_grid_above +
    12*coef(fit_y_above)[5]*x_grid_above^2
  
  # Calculate second derivatives for treatment function
  d_second_deriv_below <- 2*coef(fit_d_below)[3] +
    6*coef(fit_d_below)[4]*x_grid_below +
    12*coef(fit_d_below)[5]*x_grid_below^2
  
  d_second_deriv_above <- 2*coef(fit_d_above)[3] +
    6*coef(fit_d_above)[4]*x_grid_above +
    12*coef(fit_d_above)[5]*x_grid_above^2
  
  # Find maximum absolute second derivatives
  B_Y <- max(max(abs(y_second_deriv_below)), max(abs(y_second_deriv_above)))
  B_T <- max(max(abs(d_second_deriv_below)), max(abs(d_second_deriv_above)))
  
  return(c(B_Y, B_T))
}

calculate_ROT2 <- function(d, cutoff) {
  # Fit quadratic polynomials on either side of the cutoff
  fit_y_below <- lm(d$Y ~ d$X + I(d$X^2), subset = which(d$X < cutoff))
  fit_y_above <- lm(d$Y ~ d$X + I(d$X^2), subset = which(d$X >= cutoff))
  fit_d_below <- lm(d$D ~ d$X + I(d$X^2), subset = which(d$X < cutoff))
  fit_d_above <- lm(d$D ~ d$X + I(d$X^2), subset = which(d$X >= cutoff))
  
  # Get the second derivatives (twice the coefficient of x²)
  y_second_deriv_below <- 2 * coef(fit_y_below)[3]
  y_second_deriv_above <- 2 * coef(fit_y_above)[3]
  d_second_deriv_below <- 2 * coef(fit_d_below)[3]
  d_second_deriv_above <- 2 * coef(fit_d_above)[3]
  
  # Find maximum absolute second derivatives and multiply by 2
  B_Y <- 2 * max(abs(y_second_deriv_below), abs(y_second_deriv_above))
  B_T <- 2 * max(abs(d_second_deriv_below), abs(d_second_deriv_above))
  
  return(c(B_Y, B_T))
}

calculate_true_parameters <- function(dgp, setup) {
  if (dgp == 1) {
    # Get the true treatment effect from original matrices
    true_treatment_effect <- coeffstreat_lm[1] - coeffscont_lm[1]
    
    # Define second derivative of f(x) = x^2 (with a jump at x=0)
    # The second derivative of x^2 is constant 2
    second_deriv_f <- function(x) {
      rep(2, length(x))  # Return 2 for all values of x
    }
    
    # Compute B_Y_true and B_T_true for different setups
    x_grid <- seq(-3, 3, by=0.01)  # Grid for evaluating second derivatives
    
    if (setup == 1) {
      B_Y_true <- 2  # Constant second derivative of x^2
      B_T_true <- 0  # Step function still has zero second derivative
    } else if (setup == 2) {
      B_Y_true <- 2  # Same constant second derivative
      B_T_true <- 0  # Piecewise linear treatment remains unchanged
    } else if (setup == 3) {
      # Setup 3 involves pnorm, so treatment function second derivative remains the same
      d2pnorm <- function(x) dnorm(x) * (-x)
      B_T_true <- max(abs(d2pnorm(x_grid)))
      
      # Outcome function is now x^2, so second derivative is constant 2
      B_Y_true <- 2
    }
  } else {
    true_treatment_effect <- coeffstreat_lee[1] - coeffscont_lee[1]
    
    # Calculate true smoothness bounds for DGP 2 (polynomial model)
    x_grid <- seq(-1, 1, by=0.001)
    
    # Second derivatives of outcome polynomial
    second_deriv_y_below <- function(x) {
      2*coeffscont_lee[3] + 6*coeffscont_lee[4]*x + 12*coeffscont_lee[5]*x^2 + 20*coeffscont_lee[6]*x^3
    }
    
    second_deriv_y_above <- function(x) {
      2*coeffstreat_lee[3] + 6*coeffstreat_lee[4]*x + 12*coeffstreat_lee[5]*x^2 + 20*coeffstreat_lee[6]*x^3
    }
    
    # For treatment function, depends on the setup
    if (setup == 1) {
      B_T_true <- 0  # Step function has zero second derivative
    } else if (setup == 2) {
      B_T_true <- 0  # Piecewise linear function has zero second derivative
    } else if (setup == 3) {
      # For setup 3 with pnorm
      d2pnorm <- function(x) dnorm(x) * (-x)
      B_T_true <- max(abs(d2pnorm(x_grid)))
    }
    
    # Calculate overall bounds for outcome
    B_Y_true <- max(
      max(abs(second_deriv_y_below(x_grid[x_grid < 0]))),
      max(abs(second_deriv_y_above(x_grid[x_grid > 0])))
    )
  }
  
  # Return the parameters
  return(list(
    true_treatment_effect = true_treatment_effect,
    M_true = c(B_Y_true, B_T_true)
  ))
}