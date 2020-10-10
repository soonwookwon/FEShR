#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
double nll (const arma::mat& mu,
	    const arma::mat& Lambda,
	    const arma::mat& y,
	    Rcpp::List& M) {

  int J = M.size();
  
  double nll = 0;
  
  for (int j = 0; j < J; ++j) {
    arma::mat M_j = M[j];
    arma::colvec y_j = y.col(j);
    arma::colvec mu_j = mu.col(j);
    nll = nll + log(det(Lambda + M_j)) +
      as_scalar((y_j - mu_j).t() * inv(Lambda + M_j) * (y_j - mu_j));
  }
  return(nll / J);
}
