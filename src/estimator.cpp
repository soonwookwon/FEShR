#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::mat get_thetahat (const arma::mat& mu,
			const arma::mat& Lambda,
			const arma::mat& y,
			Rcpp::List& M,
			bool missing_cells,
			Rcpp::List& O) {

  int J = M.size();
  int T = y.n_rows;
  
  arma::mat thetahat = arma::zeros(T, J);

  if (! missing_cells) {
    for (int j = 0; j < J; ++j) {
      arma::mat M_j = M[j];
      arma::colvec y_j = y.col(j);
      arma::colvec mu_j = mu.col(j);
      arma::mat shr_mat = Lambda * inv_sympd(Lambda + M_j);
      thetahat.col(j) = (arma::eye(T,T) - shr_mat) * mu_j + shr_mat * y_j;
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
      arma::mat shr_mat = Lam_j * inv_sympd(Lam_j + M_j);
      
      thetahat.col(j) = O_j_mat.t() *
	((arma::eye(o_j, o_j) - shr_mat) * mu_j + shr_mat * y_j);
    }
  }
    
  
  return thetahat;
}
