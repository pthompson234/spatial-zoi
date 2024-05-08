#include "TMB.hpp"
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_SCALAR( t_weights_log_sum ); // sum of seasonality values
  DATA_VECTOR( t_covariates_sum ); // "traditional" covariate values at time of detection. first entry is # of data points (for intercept)
  DATA_MATRIX( t_distances ); // distances from each site to nearby sites
  DATA_ARRAY( t_spcov_values ); // special covariate values at all focal locations for each detection (can be a 3-D array if >1 sp. cov.)
  DATA_MATRIX( t_alcov_values ); // covariates affecting reach of special covariates
  DATA_VECTOR( site_weights ); // seasonality index for each site history
  DATA_MATRIX( site_covariates ); // "traditional" covariates corresponding to each site interval. first column should be all 1's (for intercept)
  DATA_MATRIX( site_distances ); // distances from each site to nearby sites
  DATA_ARRAY( site_spcov_values ); // special covariate values at all focal locations for each site (can be a 3-D array if >1 sp. cov.)
  DATA_MATRIX( site_alcov_values ); // covariates affecting reach of spatial covariates for sites
  
  // Parameters
  PARAMETER_VECTOR( lambda_i ); // coefficients for intensity covariates, including intercept
  PARAMETER_MATRIX( beta_i ); // coefficients for strength of "special" covariates; includes same covariates as alpha
  PARAMETER_MATRIX( logalpha_i ); // coefficients for distance decay of "special" covariates; log-transformed so always >0; first column is intercept
  //      A note: should we allow there to be different alcovs for each spcov? In other words, t_alcov_values would become a 3-D array, and then we'd have 
  // to account for potentially having different #'s of parameters in logalpha_i matrix. Could do this with a DATA_IMATRIX of 1's and 0's?
  
  int n_site_changes = site_weights.size();
  int n_betas = beta_i.rows();
  int n_detections = t_distances.rows();
  int n_site_connections = t_distances.cols(); // This includes the site itself!! First row of distances will be 0
  
  Type val = -t_weights_log_sum - (lambda_i * t_covariates_sum).sum(); // negative log likelihood
  
  vector<Type> site_total_likelihoods(n_site_changes);
  
  for (int i = 0; i < n_detections; i++) {
    for (int r = 0; r < n_betas; r++) {
      vector<Type> this_alpha = logalpha_i.row(r);
      vector<Type> this_alpha_covs = t_alcov_values.row(i);
      Type this_alpha_i = exp((this_alpha * this_alpha_covs).sum()); // take exp so it's always greater than 0 -> decreasing spatial effect
      
      vector<Type> this_beta = beta_i.row(r);
      Type this_beta_i = (this_beta * this_alpha_covs).sum();
      
      Type h_s = 0.0;
      for (int s = 0; s < n_site_connections; s++) {
        h_s += (pow(t_distances(i, s) + 1.0, -this_alpha_i) * t_spcov_values(i, s, r));
        // h_s += (exp(-this_alpha_i * t_distances(i, s)) * t_spcov_values(i, s, r));
      }
      val -= (h_s * this_beta_i);
      // val -= (log(h_s) * this_beta_i);
    }
  }
  
  for (int j = 0; j < n_site_changes; j++) {
    Type lik_tot = 0;
    for (int r = 0; r < n_betas; r++) {
      vector<Type> this_alpha = logalpha_i.row(r);
      vector<Type> this_alpha_covs = site_alcov_values.row(j);
      Type this_alpha_i = exp((this_alpha * this_alpha_covs).sum()); // take exp so it's always greater than 0 -> decreasing spatial effect
      
      vector<Type> this_beta = beta_i.row(r);
      Type this_beta_i = (this_beta * this_alpha_covs).sum();
      
      Type h_s = 0.0;
      for (int s = 0; s < n_site_connections; s++) {
        h_s += (pow(site_distances(j, s) + 1.0, -this_alpha_i) * site_spcov_values(j, s, r));
        // h_s += (exp(-this_alpha_i * site_distances(j, s)) * site_spcov_values(j, s, r));
      }
      lik_tot += (h_s * this_beta_i);
      // lik_tot += (log(h_s) * this_beta_i);
    }
    
    vector<Type> this_row_of_covs = site_covariates.row(j);
    site_total_likelihoods(j) = lik_tot + (this_row_of_covs * lambda_i).sum();
  }
  
  val += (exp(site_total_likelihoods) * site_weights).sum();
  
  // Reporting
  REPORT( val );
  REPORT( lambda_i );
  REPORT( beta_i );
  REPORT( logalpha_i );
  
  // for debugging purposes
  // REPORT( val_aftercases );
  
  ADREPORT( lambda_i );
  ADREPORT( beta_i );
  ADREPORT( logalpha_i );
  
  return val;
}
