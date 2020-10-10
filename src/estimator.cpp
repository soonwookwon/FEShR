#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::mat get_thetahat (const arma::mat& mu,
			const arma::mat& Lambda,
			const arma::mat& y,
			Rcpp::List& M) {

  int J = M.size();
  int T = y.n_rows;
  
  arma::mat thetahat = arma::zeros(T, J);
  
  for (int j = 0; j < J; ++j) {
    arma::mat M_j = M[j];
    arma::colvec y_j = y.col(j);
    arma::colvec mu_j = mu.col(j);
    arma::mat shr_mat = Lambda * inv_sympd(Lambda + M_j);
    thetahat.col(j) = (arma::eye(T,T) - shr_mat) * mu_j + shr_mat * y_j;
  }
  
  return(thetahat);
}
