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

fmt_align <- function(x, dp, width = 9) {
  sprintf(paste0("%", width, ".", dp, "f"), x)
}