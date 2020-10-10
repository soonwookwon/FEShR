#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::mat nll_deriv (const arma::mat& mu,
		     const arma::mat& Lambda,
		     const arma::mat& y,
		     Rcpp::List& M) {

  int J = M.size();
  int T = y.n_rows;
    
  arma::mat nll_deriv = arma::zeros(T,T);
  
  for (int j = 0; j < J; ++j) {
    arma::mat M_j = M[j];
    arma::colvec y_j = y.col(j);
    arma::colvec mu_j = mu.col(j);
    arma::mat inv_Lam_Mj = inv_sympd(Lambda + M_j); 
    // double URE_j = -2 * trace(inv_Lam_Mj * M_j * M_j);
    // URE_j = URE_j + accu(square(M_j * inv_Lam_Mj * (y_j - mu_j)));
    nll_deriv = nll_deriv +
      (arma::eye(T, T)  - inv_Lam_Mj * (y_j - mu_j) * (y_j - mu_j).t()) *
      inv_Lam_Mj;
  }
  return(nll_deriv / J);
}
