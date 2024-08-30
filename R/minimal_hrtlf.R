#' @export
minimal_hrtlf<- function(Number_of_Factors){
  factor_levels <- rep(2, (Number_of_Factors - 1))
  result <- minimal.factorial(factor_levels)
  run_sequence <- result$run.sequence
  new_column <- apply(run_sequence, 1, prod)
  new_column <- matrix(new_column, ncol = 1)
  run_sequence <- cbind(run_sequence, new_column)
  new_perFactorChange <- sum(result$perFactorChange)
  updated_perFactorChange <- c(result$perFactorChange, new_perFactorChange)
  updated_totalChange <- sum(updated_perFactorChange)
  mean_column <- matrix(1, nrow = nrow(run_sequence), ncol = 1)
  design_matrix <- cbind(mean_column, run_sequence)
  A1 <- t(design_matrix) %*% design_matrix
  D <- det(A1)
  Orthogonal_polynomial <- poly(1:nrow(run_sequence), degree = 1, simple = TRUE)
  A2 <- t(design_matrix) %*% Orthogonal_polynomial
  A3 <- t(A2)
  A4 <- t(Orthogonal_polynomial) %*% Orthogonal_polynomial
  Dt <- round(det(rbind(cbind(A1, A2), cbind(A3, A4))), 2)
  Trend_factor <- round((Dt / D)^(1 / ncol(design_matrix)), 2)
  return(list("Run_Sequence" = run_sequence,
              "Factor_Wise_Change" = updated_perFactorChange,
              "Total_Change" = updated_totalChange,
              "Trend_Factor" = Trend_factor))
}
