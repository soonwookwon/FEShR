#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::mat nll_deriv (const arma::mat& mu,
		     const arma::mat& Lambda,
		     const arma::mat& y,
		     Rcpp::List& M,
		     bool missing_cells,
		     Rcpp::List& O) {

  int J = M.size();
  int T = y.n_rows;
    
  arma::mat nll_deriv = arma::zeros(T,T);

  if (!missing_cells) {
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
  } else {
    arma::mat I = arma::eye(T,T);
    for (int j = 0; j < J; ++j) {
      arma::uvec O_j = O[j];
      O_j = O_j - 1;
      arma::mat O_j_mat = I.rows(O_j);
      double o_j = O_j.n_elem;
      arma::mat M_j = M[j];
      M_j = M_j.submat(O_j, O_j);
      arma::colvec y_j = y.col(j);
      y_j = y_j.elem(O_j);
      arma::colvec mu_j =  mu.col(j);
      mu_j = mu_j.elem(O_j);
      arma::mat Lam_j = Lambda.submat(O_j, O_j);
      arma::mat inv_Lam_Mj = inv_sympd(Lam_j + M_j); 
      // double URE_j = -2 * trace(inv_Lam_Mj * M_j * M_j);
      // URE_j = URE_j + accu(square(M_j * inv_Lam_Mj * (y_j - mu_j)));
      nll_deriv = nll_deriv +
	O_j_mat.t() * (arma::eye(o_j,o_j) -
		       inv_Lam_Mj * (y_j - mu_j) * (y_j - mu_j).t()) *
	inv_Lam_Mj * O_j_mat;
    }
  }
    
  
  return(nll_deriv / J);
}
