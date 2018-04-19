simulateDeSeq2 <- function(num_transcripts, num_replicates, percent_with_effect, effect_size, num_conditions = 2) {
  dispersion_sigma = abs(rnorm(1,0,1));
  #dispersion_sigma = 1;
  #asymptotic_dispersion = abs(rcauchy(1,0,0.2)) + 0.01;
  asymptotic_dispersion = 0.1
  
  #dispersion_from_mean_count = abs(rcauchy(1,0,10));
  dispersion_from_mean_count = 3
  
  
  #intercept_sigma = abs(rnorm(1,0,1));
  #intercept_mean = abs(rnorm(1,0,1)) + 2;
  #coefficients_sigma = abs(rcauchy(num_conditions - 1, 0, 0.2)) + 0.1;

  coefficients = array(0, c(num_transcripts, num_conditions));
  for(f in 1:(num_conditions - 1)) {
    signs <- (-0.5 + rbinom(num_transcripts, 1, 0.5)) * 2
    coefficients[,f] = rbinom(num_transcripts, 1, percent_with_effect) * signs * effect_size; 
  }
  #The intercept
  coefficients[,num_conditions] = runif(num_transcripts, 3, 5) #abs(rnorm(num_transcripts, 0 , intercept_sigma)) + intercept_mean;

  num_samples = num_replicates * num_conditions;
  
  design_matrix = array(0, c(num_conditions, num_conditions));
  for(c in 1:(num_conditions - 1)) {
    design_matrix[c,c] = 1;
  }
  #The intercept coefficients
  design_matrix[num_conditions, ] = 1;
  
  sample_conditions = integer(num_samples);
  for(c in 1:num_conditions) {
    startIndex = (c - 1) * num_replicates + 1;
    sample_conditions[startIndex:(startIndex + num_replicates - 1)] = c;
  }
  
  
  normalization = rep_len(1, num_samples);
  
  unnormalized_count_means = 2^(coefficients %*% design_matrix);
  
  mean_normalized_read_counts = array(0, c(num_transcripts, num_samples));
  for(s in 1:num_samples) {
    mean_normalized_read_counts[,s] = unnormalized_count_means[,sample_conditions[s]] * normalization[s];
  }
  
    
  dispersion_trend = (dispersion_from_mean_count / rowMeans(mean_normalized_read_counts)) + asymptotic_dispersion;  
  dispersion = rlnorm(num_transcripts,  dispersion_trend, dispersion_sigma)
  #dispersion = 1;

  counts = array(as.integer(0), c(num_samples, num_transcripts))
  for(s in 1:num_samples) {
    counts[s,] = as.integer(rnbinom(num_transcripts, mu = mean_normalized_read_counts[,s], size = dispersion));
  }

  return (
    list(
      observed = list(
        counts = t(counts),
        normalization = normalization,
        design_matrix = array(t(design_matrix[1:num_conditions - 1,]), c(num_conditions, num_conditions - 1)),
        sample_info = data.frame(id = 1:num_samples, condition = as.factor(sample_conditions))
      ),
      true = list(
        mean_normalized_read_counts = mean_normalized_read_counts,
        coefficients = coefficients,
        dispersion = dispersion,
        asymptotic_dispersion = asymptotic_dispersion,
        dispersion_from_mean_count = dispersion_from_mean_count,
        dispersion_sigma = dispersion_sigma
        #intercept_mean = intercept_mean,
        #intercept_sigma = intercept_sigma,
        #coefficients_sigma = array(coefficients_sigma, num_conditions - 1)
      )
    )
  )
}