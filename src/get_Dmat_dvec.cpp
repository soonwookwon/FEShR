#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
Rcpp::List get_Dmat_dvec (const arma::mat& Lambda,
			  const arma::mat& y,
			  Rcpp::List& M) {
  
  int J = M.size();
  int T = y.n_rows;

  arma::mat Dmat = arma::zeros(T, T);
  arma::colvec dvec = arma::zeros(T);
  
  for (int j = 0; j < J; ++j) {
    arma::mat M_j = M[j];
    arma::colvec y_j = y.col(j);
    arma::mat inv_Lam_Mj = inv_sympd(Lambda + M_j);
    arma::mat inv_Lam_Mj_Mj = inv_Lam_Mj * M_j;
    arma::mat common_part = inv_Lam_Mj_Mj * inv_Lam_Mj_Mj.t();
    Dmat = Dmat + common_part;
    dvec = dvec + common_part * y_j;
  }

  return Rcpp::List::create(Rcpp::Named("Dmat")=Dmat/J,
			    Rcpp::Named("dvec")=dvec/J);
}
