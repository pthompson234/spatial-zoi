#include "TMB.hpp"
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_MATRIX( t_covariates ); // "traditional" covariate values at time of detection. first entry is # of data points (for intercept)
  DATA_MATRIX( t_distances ); // distances from each site to nearby sites
  DATA_MATRIX( t_spcov_values ); // special covariate values at all focal locations for each detection
  DATA_VECTOR( site_weights ); // length of time for each "integral"
  DATA_MATRIX( site_covariates ); // "traditional" covariates corresponding to each site interval. first column should be all 1's (for intercept)
  DATA_MATRIX( site_distances ); // distances from each site to nearby sites
  DATA_MATRIX( site_spcov_values ); // special covariate values at all focal locations for each site
  DATA_VECTOR( connection_ones ); // a vector of length n_site_connections with all entries being 1; for easy summing; better to make it ahead of time
  
  // Parameters
  PARAMETER_VECTOR( lambda_i ); // coefficients for intensity covariates, including intercept
  PARAMETER( logalpha_i ); // coefficients for distance decay of "special" covariates; log-transformed so always >0; first column is intercept
  PARAMETER( lognu_i ); // coefficients for how the human and environmental portions of the model are weighted
  
  Type alpha_i = exp(logalpha_i);
  Type nu_i = exp(lognu_i);
  
  matrix<Type> t_cor_weights = pow(t_distances.array(), -alpha_i);
  vector<Type> t_weighted_cor_values_sum = (t_cor_weights * t_spcov_values.transpose()).diagonal();
  vector<Type> t_cor_weights_sum = t_cor_weights * connection_ones;
  vector<Type> t_cor_weights_01 = Type(1.0) - exp(-nu_i * t_cor_weights_sum);

  // negative log likelihood - first for the detections
  vector<Type> val1 = (t_cor_weights_01 - Type(1.0)) * (t_covariates * lambda_i);
  vector<Type> val2 = t_cor_weights_01 * log(t_weighted_cor_values_sum / t_cor_weights_sum);
  Type val = (val1 - val2).sum();

  matrix<Type> site_cor_weights = pow(site_distances.array(), -alpha_i);
  vector<Type> site_weighted_cor_values_sum = (site_cor_weights * site_spcov_values.transpose()).diagonal();
  vector<Type> site_cor_weights_sum = site_cor_weights * connection_ones;
  vector<Type> site_cor_weights_01 = Type(1.0) - exp(-nu_i * site_cor_weights_sum);

  val += (site_weights * exp((Type(1.0) - site_cor_weights_01) * (site_covariates * lambda_i) + site_cor_weights_01 * log(site_weighted_cor_values_sum / site_cor_weights_sum))).sum();
  // likelihood portion from "non-detections"
  
  // Reporting
  REPORT( val );
  REPORT( lambda_i );
  REPORT( logalpha_i );
  REPORT( lognu_i );
  
  // for debugging purposes
  // REPORT( val_aftercases );
  
  ADREPORT( lambda_i );
  ADREPORT( logalpha_i );
  ADREPORT( lognu_i );
  
  return val;
}
