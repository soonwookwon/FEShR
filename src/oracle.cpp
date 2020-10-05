#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
double oracle_loss_cpp (const arma::mat& thetas,
			const arma::mat& mu,
			const arma::mat& Lambda,
			const arma::mat& y,
			Rcpp::List& M) {

  int J = M.size();
  int T = y.n_rows;
  
  double loss = 0;
  
  for (int j = 0; j < J; ++j) {
    arma::mat M_j = M[j];
    arma::colvec y_j = y.col(j);
    arma::colvec theta_j = thetas.col(j);
    arma::colvec mu_j = mu.col(j);
    arma::mat shr_mat = Lambda * inv_sympd(Lambda + M_j);
    arma::colvec thetahat_j = (arma::eye(T,T) - shr_mat) * mu_j + shr_mat * y_j;
    // double URE_j = -2 * trace(inv_Lam_Mj * M_j * M_j);
    // URE_j = URE_j + accu(square(M_j * inv_Lam_Mj * (y_j - mu_j)));
    loss = loss + accu(square(thetahat_j - theta_j));
  }
  
  return(loss / J);
}


//[[Rcpp::export]]
arma::mat oracle_cpp_deriv (const arma::mat& thetas,
			    const arma::mat& mu,
			    const arma::mat& Lambda,
			    const arma::mat& y,
			    Rcpp::List& M) {

  int J = M.size();
  int T = y.n_rows;
  
  arma::mat deriv = arma::zeros(T,T);
  
  for (int j = 0; j < J; ++j) {
    arma::mat M_j = M[j];
    arma::colvec y_j = y.col(j);
    arma::colvec theta_j = thetas.col(j);
    arma::colvec mu_j = mu.col(j);
    arma::mat inv_Lam_Mj = inv_sympd(Lambda + M_j);
    arma::mat Mj_inv_Lam_Mj = M_j * inv_Lam_Mj;
    deriv = deriv + 2 * Mj_inv_Lam_Mj.t() * (y_j - theta_j) * (y_j - mu_j).t() *
      inv_Lam_Mj - 2 * inv_Lam_Mj * (y_j - mu_j) * (y_j - mu_j).t() *
      Mj_inv_Lam_Mj.t() * Mj_inv_Lam_Mj;
  }
  
  return(deriv / J);
}
