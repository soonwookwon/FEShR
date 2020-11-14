#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
double nll (const arma::mat& mu,
	    const arma::mat& Lambda,
	    const arma::mat& y,
	    Rcpp::List& M,
	    bool missing_cells,
	    Rcpp::List& O) {

  int J = M.size();
  
  double nll = 0;

  if (!missing_cells) {
    for (int j = 0; j < J; ++j) {
      arma::mat M_j = M[j];
      arma::colvec y_j = y.col(j);
      arma::colvec mu_j = mu.col(j);
      nll = nll + log(det(Lambda + M_j)) +
	as_scalar((y_j - mu_j).t() * inv(Lambda + M_j) * (y_j - mu_j));
    }
  } else {

    for (int j = 0; j < J; ++j) {
      arma::uvec O_j = O[j];
      O_j = O_j - 1;
      arma::mat M_j = M[j];
      M_j = M_j.submat(O_j, O_j);
      arma::colvec y_j = y.col(j);
      y_j = y_j.elem(O_j);
      arma::colvec mu_j =  mu.col(j);
      mu_j = mu_j.elem(O_j);
      arma::mat Lam_j = Lambda.submat(O_j, O_j);
      nll = nll + log(det(Lam_j + M_j)) +
	as_scalar((y_j - mu_j).t() * inv(Lam_j + M_j) * (y_j - mu_j));
    }
    
  }
   
  return(nll / J);
}
