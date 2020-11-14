#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
double URE (const arma::mat& mu,
	    const arma::mat& Lambda,
	    const arma::mat& y,
	    Rcpp::List& M,
	    const arma::mat& W,
	    bool missing_cells,
	    Rcpp::List& O) {

  int J = M.size();
  
  double URE = 0;

  if (!missing_cells) {
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
  } else {
    
    for (int j = 0; j < J; ++j) {
      arma::uvec O_j = O[j];
      O_j = O_j - 1;
      double o_j = O_j.n_elem;
      arma::mat M_j = M[j];
      M_j = M_j.submat(O_j, O_j);
      arma::colvec y_j = y.col(j);
      y_j = y_j.elem(O_j);
      arma::colvec mu_j =  mu.col(j);
      mu_j = mu_j.elem(O_j);
      arma::mat Lam_j = Lambda.submat(O_j, O_j);
      arma::mat W_j = W.submat(O_j, O_j);
      arma::mat inv_Lam_Mj = inv_sympd(Lam_j + M_j);
      arma::colvec  Mj_inv_Lam_Mj_yj = M_j * inv_Lam_Mj * (y_j - mu_j);
      // double URE_j = -2 * trace(inv_Lam_Mj * M_j * M_j);
      // URE_j = URE_j + accu(square(M_j * inv_Lam_Mj * (y_j - mu_j)));
      URE = URE - (2 / o_j) * trace(inv_Lam_Mj * M_j * W_j * M_j) +
	(1 / o_j) * arma::as_scalar(Mj_inv_Lam_Mj_yj.t() * W_j * Mj_inv_Lam_Mj_yj);
    }
  }    
 
  return(URE / J);
}
