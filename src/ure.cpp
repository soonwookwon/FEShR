#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
double URE (const arma::mat& mu,
	    const arma::mat& Lambda,
	    const arma::mat& y,
	    Rcpp::List& M,
	    const arma::mat& W) {

  int J = M.size();
  
  double URE = 0;

  
  for (int j = 0; j < J; ++j) {
    arma::mat M_j = M[j];
    arma::colvec y_j = y.col(j);
    arma::colvec mu_j = mu.col(j);
    arma::mat inv_Lam_Mj = inv_sympd(Lambda + M_j);
    arma::colvec  Mj_inv_Lam_Mj_yj = M_j * inv_Lam_Mj * (y_j - mu_j);
    // double URE_j = -2 * trace(inv_Lam_Mj * M_j * M_j);
    // URE_j = URE_j + accu(square(M_j * inv_Lam_Mj * (y_j - mu_j)));
    URE = URE - 2 * trace(inv_Lam_Mj * M_j * W * M_j) +
      arma::as_scalar(Mj_inv_Lam_Mj_yj.t() * W * Mj_inv_Lam_Mj_yj);
  }
  return(URE / J);
}
