#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::mat opt_mu_Lambda_nll (const arma::mat& Lambda,
			     const arma::mat& y,
			     Rcpp::List& M) {
  
  int J = M.size();
  int T = y.n_rows;

  arma::colvec opt_mu = arma::zeros(T);
  
  for (int j = 0; j < J; ++j) {
    arma::mat M_j = M[j];
    arma::colvec y_j = y.col(j);
    arma::mat  inv_Lam_Mj = inv_sympd(Lambda + M_j);
 
    opt_mu = opt_mu + inv_Lam_Mj * y_j;
  }

  arma::mat opt_mu_mat(T, J, arma::fill::zeros);
  opt_mu_mat.each_col() = opt_mu / J;
    
  return opt_mu_mat;
}


//[[Rcpp::export]]
arma::mat opt_mu_Lambda_ol (const arma::mat& thetas,
			    const arma::mat& Lambda,
			    const arma::mat& y,
			    Rcpp::List& M) {
  
  int J = M.size();
  int T = y.n_rows;

  arma::mat denom = arma::zeros(T, T);
  arma::colvec numer = arma::zeros(T);
  
  for (int j = 0; j < J; ++j) {
    arma::mat M_j = M[j];
    arma::colvec y_j = y.col(j);
    arma::colvec theta_j = thetas.col(j);
    arma::mat inv_Lam_Mj = inv_sympd(Lambda + M_j);
    arma::mat shr_mat = Lambda * inv_Lam_Mj;
    arma::mat inv_Lam_Mj_Mj = (arma::eye(T,T) - shr_mat).t();

    denom = denom +  inv_Lam_Mj_Mj *  inv_Lam_Mj_Mj.t();
    numer = numer + inv_Lam_Mj_Mj * (theta_j - shr_mat * y_j);
  }

  arma::colvec opt_mu = inv_sympd(denom) * numer; 
  arma::mat opt_mu_mat(T, J, arma::fill::zeros);
 
  opt_mu_mat.each_col() = opt_mu;
    
  return opt_mu_mat;
}
