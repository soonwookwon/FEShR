#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::colvec opt_LamTmT_Lambda(const arma::mat& Lambda,
			       const arma::mat& y,
			       Rcpp::List& M,
			       bool missing_cells,
			       Rcpp::List& O) {

  int J = M.size();
  int T = y.n_rows;
  
  arma::mat Lam_TmT_denom = arma::zeros(T-1, T-1);
  arma::colvec Lam_TmT_numer = arma::zeros(T-1);
  
  if (!missing_cells) {

    for (int j = 0; j < J; ++j) {
      arma::mat M_j = M[j];
      M_j = M_j.submat(0, 0, T-2, T-2);
      arma::colvec y_j = y.col(j);
      y_j = y_j.head(T-1);
      double y_jt = arma::as_scalar(y_j.tail(1));
      
      arma::mat inv_Lam_Mj = inv_sympd(Lambda + M_j);
      arma::colvec y_til_j = inv_Lam_Mj * y_j;
      Lam_TmT_denom = Lam_TmT_denom + y_til_j * y_til_j.t();
      Lam_TmT_numer = Lam_TmT_numer + y_jt * y_til_j;
    }
  }
  arma::colvec Lam_TmT = arma::solve(Lam_TmT_denom, Lam_TmT_numer);

  return Lam_TmT;
}


//[[Rcpp::export]]
double UPE (const arma::mat& Lambda,
	    const arma::mat& y,
	    Rcpp::List& M,
	    bool missing_cells,
	    Rcpp::List& O,
	    double bounds) {

  int J = M.size();
  int T = y.n_rows;
  
  double UPE = 0;
  
  // arma::mat Lam_TmT_denom = arma::zeros(T-1, T-1);
  // arma::colvec Lam_TmT_numer = arma::zeros(T-1);
  
  if (!missing_cells) {

    // for (int j = 0; j < J; ++j) {
    //   arma::mat M_j = M[j];
    //   M_j = M_j.submat(0, 0, T-2, T-2);
    //   arma::colvec y_j = y.col(j);
    //   y_j = y_j.head(T-1);
    //   double y_jt = arma::as_scalar(y_j.tail(1));

    //   arma::mat inv_Lam_Mj = inv_sympd(Lambda + M_j);
    //   arma::colvec y_til_j = inv_Lam_Mj * y_j;
    //   Lam_TmT_denom = Lam_TmT_denom + y_til_j * y_til_j.t();
    //   Lam_TmT_numer = Lam_TmT_numer + y_jt * y_til_j;
    // }

    // arma::colvec Lam_TmT = arma::solve(Lam_TmT_denom, Lam_TmT_numer);

    arma::colvec Lam_TmT = opt_LamTmT_Lambda(Lambda, y, M, missing_cells, O);
      
    for (int j = 0; j < J; ++j) {
      arma::mat M_j = M[j];
      M_j = M_j.submat(0, 0, T-2, T-2);
      arma::colvec y_j = y.col(j);
      y_j = y_j.head(T-1);
      double y_jt = arma::as_scalar(y_j.tail(1));
       
      arma::mat inv_Lam_Mj = inv_sympd(Lambda + M_j);
            
      UPE = UPE +
	arma::as_scalar(arma::square(Lam_TmT.t() * inv_Lam_Mj * y_j - y_jt));
    }
  } // else {
    
  //   for (int j = 0; j < J; ++j) {
  //     arma::uvec O_j = O[j];
  //     O_j = O_j - 1;
  //     double o_j = O_j.n_elem;
  //     arma::mat M_j = M[j];
  //     M_j = M_j.submat(O_j, O_j);
  //     arma::colvec y_j = y.col(j);
  //     y_j = y_j.elem(O_j);
  //     arma::colvec mu_j =  mu.col(j);
  //     mu_j = mu_j.elem(O_j);
  //     arma::mat Lam_j = Lambda.submat(O_j, O_j);
  //     arma::mat W_j = W.submat(O_j, O_j);
  //     arma::mat inv_Lam_Mj = inv_sympd(Lam_j + M_j);
  //     arma::colvec  Mj_inv_Lam_Mj_yj = M_j * inv_Lam_Mj * (y_j - mu_j);
  //     // double URE_j = -2 * trace(inv_Lam_Mj * M_j * M_j);
  //     // URE_j = URE_j + accu(square(M_j * inv_Lam_Mj * (y_j - mu_j)));
  //     URE = URE - (2 / o_j) * trace(inv_Lam_Mj * M_j * W_j * M_j) +
  // 	(1 / o_j) * arma::as_scalar(Mj_inv_Lam_Mj_yj.t() * W_j * Mj_inv_Lam_Mj_yj);
  //   }
  // }    
 
  return(UPE / J);
}


//[[Rcpp::export]]
arma::mat get_thetahat_fc(const arma::colvec& LamTmT,
			  const arma::mat& Lambda,
			  const arma::mat& y,
			  Rcpp::List& M,
			  bool missing_cells,
			  Rcpp::List& O) {

  int J = M.size();
  int T = y.n_rows;
  
  arma::mat thetahat = arma::zeros(1, J);

  if (! missing_cells) {
    for (int j = 0; j < J; ++j) {
      arma::mat M_j = M[j];
      M_j = M_j.submat(0, 0, T-2, T-2);
      arma::colvec y_j = y.col(j);
      y_j = y_j.tail(T-1);
      arma::mat inv_Lam_Mj = inv_sympd(Lambda + M_j);
      thetahat.col(j) = LamTmT.t() * inv_Lam_Mj * y_j;
    }
  } 
  
  return thetahat;
}
