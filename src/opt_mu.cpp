#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


//[[Rcpp::export]]
arma::mat gam_to_mu (const arma::colvec& gam_opt,
		     Rcpp::List& Z,
		     int T) {
  
  int J = Z.size();
  
  arma::mat mu = arma::zeros(T, J);
  
  for (int j = 0; j < J; ++j) {
    arma::mat Z_j = Z[j];
    mu.col(j) = Z_j * gam_opt; 
  }

  return mu;
}


//[[Rcpp::export]]
Rcpp::List get_Dmat_dvec (const arma::mat& Lambda,
			  const arma::mat& y,
			  Rcpp::List& M,
			  const arma::mat& W,
			  bool missing_cells,
			  Rcpp::List& O,
			  Rcpp::List& Z) {
  
  int J = M.size();
  int T = y.n_rows;
  int Z_length = Z.size();

  arma::mat Dmat = arma::zeros(T, T);
  arma::colvec dvec = arma::zeros(T);

  if (Z_length != 0) {
     arma::mat Z_1 = Z[1];
     int k = Z_1.n_cols;
     Dmat = Dmat.zeros(k, k);
     dvec = dvec.zeros(k);
  }

  if (!missing_cells) {
    for (int j = 0; j < J; ++j) {
      arma::mat M_j = M[j];
      arma::colvec y_j = y.col(j);
      arma::mat inv_Lam_Mj = inv_sympd(Lambda + M_j);
      arma::mat inv_Lam_Mj_Mj = inv_Lam_Mj * M_j;
      arma::mat common_part = inv_Lam_Mj_Mj * W * inv_Lam_Mj_Mj.t();
      if (Z_length != 0) {
	arma::mat Z_j = Z[j];
        Dmat = Dmat + Z_j.t() * common_part * Z_j;
	dvec = dvec + Z_j.t() * common_part * y_j;
      } else {
	Dmat = Dmat + common_part;
	dvec = dvec + common_part * y_j;
      }
    }
  } else {
    arma::mat I = arma::eye(T, T);
    for (int j = 0; j < J; ++j) {
      arma::uvec O_j = O[j];
      O_j = O_j - 1;
      double o_j = O_j.n_elem;
      arma::mat M_j = M[j];
      M_j = M_j.submat(O_j, O_j);
      arma::colvec y_j = y.col(j);
      y_j = y_j.elem(O_j);
      arma::mat Lam_j = Lambda.submat(O_j, O_j);
      arma::mat W_j = W.submat(O_j, O_j);
      arma::mat inv_Lam_Mj = inv_sympd(Lam_j + M_j);
      arma::mat inv_Lam_Mj_Mj = inv_Lam_Mj * M_j;
      arma::mat common_part = inv_Lam_Mj_Mj * W_j * inv_Lam_Mj_Mj.t();
      if (Z_length != 0) {
	arma::mat Z_j = Z[j];
	Z_j = Z_j.rows(O_j);
        Dmat = Dmat + (1 / o_j) * Z_j.t() * common_part * Z_j;
	dvec = dvec + (1 / o_j) * Z_j.t() * common_part * y_j;
      } else {
	arma::mat O_j_mat = I.rows(O_j);
	Dmat = Dmat + (1 / o_j) * O_j_mat.t() * common_part * O_j_mat;
	dvec = dvec + (1 / o_j) * O_j_mat.t() * common_part * y_j;
      }      
    }
  }

  return Rcpp::List::create(Rcpp::Named("Dmat")=Dmat/J,
			    Rcpp::Named("dvec")=dvec/J);
}


//[[Rcpp::export]]
arma::mat opt_mu_Lambda_nll (const arma::mat& Lambda,
			     const arma::mat& y,
			     Rcpp::List& M,
			     bool missing_cells,
			     Rcpp::List& O,
			     Rcpp::List& Z) {
  int J = M.size();
  int T = y.n_rows;
  int Z_length = Z.size();
 
  arma::mat denom = arma::zeros(T,T);
  arma::colvec numer = arma::zeros(T);

  if (Z_length != 0) {
    arma::mat Z_1 = Z[1];
    int k = Z_1.n_cols;  
    denom = denom.zeros(k, k);
    numer = numer.zeros(k);
  }
  
  if (!missing_cells) {
    
    for (int j = 0; j < J; ++j) {
      arma::mat M_j = M[j];
      arma::colvec y_j = y.col(j);
      arma::mat  inv_Lam_Mj = inv_sympd(Lambda + M_j);

      if (Z_length != 0) {
	arma::mat Z_j = Z[j];
	denom = denom + Z_j.t() * inv_Lam_Mj * Z_j;
	numer = numer + Z_j.t() * inv_Lam_Mj * y_j;
      } else {
	denom = denom + inv_Lam_Mj;
	numer = numer + inv_Lam_Mj * y_j;
      }
    }
  } else {
    arma::mat I = arma::eye(T,T);
    for (int j = 0; j < J; ++j) {
      arma::uvec O_j = O[j];
      O_j = O_j - 1;
      arma::mat O_j_mat = I.rows(O_j);
      arma::mat M_j = M[j];
      M_j = M_j.submat(O_j, O_j);
      arma::colvec y_j = y.col(j);
      y_j = y_j.elem(O_j);
      arma::mat Lam_j = Lambda.submat(O_j, O_j);
      arma::mat inv_Lam_Mj = inv_sympd(Lam_j + M_j);

      if (Z_length != 0) {
	arma::mat Z_j = Z[j];
	Z_j = Z_j.rows(O_j);
	denom = denom + Z_j.t() * inv_Lam_Mj * Z_j;
	numer = numer + Z_j.t() * inv_Lam_Mj * y_j;
      } else {
	denom = denom + O_j_mat.t() * inv_Lam_Mj * O_j_mat;
	numer = numer + O_j_mat.t() * inv_Lam_Mj * y_j;
      }
    }

  }


  arma::colvec opt_mu = arma::solve(denom, numer);
  arma::mat opt_mu_mat(T, J, arma::fill::zeros);
    

  if (Z_length == 0) {
    opt_mu_mat.each_col() = opt_mu;
  } else {
    opt_mu_mat = gam_to_mu(opt_mu, Z, T);
  }

  return opt_mu_mat;
}


//[[Rcpp::export]]
arma::mat opt_mu_Lambda_ol (const arma::mat& thetas,
			    const arma::mat& Lambda,
			    const arma::mat& y,
			    Rcpp::List& M,  
			    const arma::mat& W,
			    bool missing_cells,
			    Rcpp::List& O,
			    Rcpp::List& Z) {
  
  int J = M.size();
  int T = y.n_rows;
  int Z_length = Z.size();

  arma::mat denom = arma::zeros(T,T);
  arma::colvec numer = arma::zeros(T);

  if (Z_length != 0) {
    arma::mat Z_1 = Z[1];
    int k = Z_1.n_cols;
    denom = denom.zeros(k, k);
    numer = numer.zeros(k);
  }
  
  if (!missing_cells) {
    
    for (int j = 0; j < J; ++j) {
      arma::mat M_j = M[j];
      arma::colvec y_j = y.col(j);
      arma::colvec theta_j = thetas.col(j);
      arma::mat inv_Lam_Mj = inv_sympd(Lambda + M_j);
      arma::mat shr_mat = Lambda * inv_Lam_Mj;
      arma::mat inv_Lam_Mj_Mj = (arma::eye(T,T) - shr_mat).t();

      if (Z_length != 0) {
	arma::mat Z_j = Z[j];
	denom = denom + Z_j.t() * inv_Lam_Mj_Mj * inv_Lam_Mj_Mj.t() * Z_j;
	numer = numer + Z_j.t() * inv_Lam_Mj_Mj * (theta_j - shr_mat * y_j);
      } else {
	denom = denom + inv_Lam_Mj_Mj * inv_Lam_Mj_Mj.t();
	numer = numer + inv_Lam_Mj_Mj * (theta_j - shr_mat * y_j);
      }
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
      arma::colvec theta_j =  thetas.col(j);
      theta_j = theta_j.elem(O_j);
    
      arma::mat Lam_j = Lambda.submat(O_j, O_j); 
      arma::mat inv_Lam_Mj = inv_sympd(Lam_j + M_j);
      arma::mat shr_mat = Lam_j * inv_Lam_Mj;
      arma::mat inv_Lam_Mj_Mj = (arma::eye(o_j,o_j) - shr_mat).t();

      
      if (Z_length != 0) {
	arma::mat Z_j = Z[j];
	Z_j = Z_j.rows(O_j);
	denom = denom +
	  (1/o_j) * Z_j.t() * inv_Lam_Mj_Mj *  inv_Lam_Mj_Mj.t() * Z_j;
	numer = numer +
	  (1/o_j) * Z_j.t() * inv_Lam_Mj_Mj * (theta_j - shr_mat * y_j);
      } else {
	denom = denom +
	  (1/o_j) * O_j_mat.t() * inv_Lam_Mj_Mj *  inv_Lam_Mj_Mj.t() * O_j_mat;
	numer = numer +
	  (1/o_j) * O_j_mat.t() * inv_Lam_Mj_Mj * (theta_j - shr_mat * y_j);
      }
      
    }
  }
  

  arma::colvec opt_mu = arma::solve(denom, numer);
  arma::mat opt_mu_mat(T, J, arma::fill::zeros);

  if (Z_length == 0) {
    opt_mu_mat.each_col() = opt_mu;
  } else {
    opt_mu_mat = gam_to_mu(opt_mu, Z, T);
  }

  return opt_mu_mat;
}  
