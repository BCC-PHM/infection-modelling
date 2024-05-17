# Define the function to minimize
test_function <- function(variance, mean, lowerCI, upperCI) {
  rate <- mean / variance
  shape <- mean ^ 2 / variance
  
  val1 <- abs(upperCI - qgamma(0.975, shape, rate = rate))
  val2 <- abs(lowerCI - qgamma(0.025, shape, rate = rate))
  
  test_val <- val1 + val2
  return(test_val)
}

gamma_params <- function(mean, lowerCI, upperCI, Distribution) {
  # Optimize the function to find the variance that minimizes test_val
  # Provide a reasonable range for variance based on the context
  if (Distribution == "Gamma") {
    result <- optimize(test_function, 
                       interval = c(1e-4, 100), 
                       tol = 1e-7, 
                       mean = mean, lowerCI = lowerCI, upperCI = upperCI)
    
    # Extract the optimal variance
    optimal_variance <- result$minimum
    
    rate <- mean / optimal_variance
    shape <- mean ^ 2 / optimal_variance
    
  } else {
    rate <- NA
    shape <- NA
  }
  return(list("rate" = rate, "shape" = shape))
}

# Example usage
# params = gamma_params(5.5, 3, 8)
