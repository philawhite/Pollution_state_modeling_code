#include <stdlib.h>            // malloc
#include <stdio.h>             // printf
#include <math.h>              // fabs, sqrt, etc.
#include <Rmath.h>              // fabs, sqrt, etc.
#include <time.h>              // time
#include <cmath>
#include <unistd.h>            // getpid
#include <string>
#include <RcppArmadillo.h>
//#include <Rcpp.h>
using namespace Rcpp;
using namespace R;
using namespace arma;  // use the Armadillo library for matrix computations

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec my_mvrnorm(vec mu, mat sigma) {
  int ncols = sigma.n_cols;
  vec Z = vectorise(randn(1,ncols) * chol(sigma));
  return mu + Z;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double sig21_update(List Z1,List X,mat bet1, List L1, mat gam1, double a_11, 
                    vec V1,double a_s, double b_s,int ns,int nt){
  
  double a_new = a_s + 0.5 * (double (ns))*(double (nt)) ;
  double b_new = b_s;
  mat bet_temp = bet1.t();
  mat gam_temp = gam1.t();
  
  for(int i=0; i < ns ; i++){
    mat X_temp = X[i]; 
    mat L_temp = L1[i]; 
    vec aV =  a_11 * ones(nt)*V1[i];
    vec Z_temp = Z1[i];
    b_new += sum(pow( Z_temp- X_temp * bet_temp.col(i) - L_temp * gam_temp.col(i) - aV,2.0)) / 2.0 ; 
  }
  double sig2 = 1/rgamma(1,a_new,1/b_new)[0];
  
  return sig2;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double sig22_update(List Z2,List X,mat bet2, List L2, mat gam2, double a_12, 
                    vec V1,vec V2,double a_s, double b_s,int ns,int nt){
  
  double a_new = a_s + 0.5 * (double (ns))*(double (nt)) ;
  double b_new = b_s;
  mat bet_temp = bet2.t();
  mat gam_temp = gam2.t();
  
  for(int i=0; i < ns ; i++){
    mat X_temp = X[i]; 
    mat L_temp = L2[i]; 
    vec aV =  (a_12*V1[i] + V2[i]) * ones(nt);
    vec Z_temp = Z2[i];
    b_new += sum(pow( Z_temp - X_temp * bet_temp.col(i) - L_temp * gam_temp.col(i) - aV,2.0)) / 2.0 ; 
  }
  double sig2 = 1/rgamma(1,a_new,1/b_new)[0];
  
  return sig2;
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat bet1_update(List Z1,List X, List L1,mat Sig_b1_inv,vec bet10,mat gam1, double a_11, 
                    vec V1,double sig21, List XtX, int p, int ns,int nt){
  mat bet_new(ns,p);
  mat gam_temp = gam1.t();
  
  for(int i = 0; i < ns ; i++){
    mat XtX_temp = XtX[i];
    mat V_new_inv = Sig_b1_inv + (1/sig21) * XtX_temp ;
    mat V_new = V_new_inv.i();
    
    mat X_temp = X[i]; 
    mat L_temp = L1[i]; 
    vec aV =  a_11 * ones(nt)*V1[i];
    vec Z_temp = Z1[i];
    vec m_new = Sig_b1_inv * bet10 + (1/sig21) * X_temp.t() * (Z_temp - L_temp * gam_temp.col(i) - aV);
    bet_new.row(i) = my_mvrnorm(V_new * m_new, V_new).t();
  }
  
  return bet_new;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat bet2_update(List Z2,List X, List L2,mat Sig_b2_inv,vec bet20,mat gam2, double a_12, 
                vec V1,vec V2,double sig22, List XtX, int p, int ns,int nt){
  mat bet_new(ns,p);
  mat gam_temp = gam2.t();
  
  for(int i = 0; i < ns ; i++){
    mat XtX_temp = XtX[i];
    mat V_new_inv = Sig_b2_inv + (1/sig22) * XtX_temp ;
    mat V_new = V_new_inv.i();
    
    mat X_temp = X[i]; 
    mat L_temp = L2[i]; 
    vec aV =  (a_12*V1[i] + V2[i]) * ones(nt);
    vec Z_temp = Z2[i];
    vec m_new = Sig_b2_inv * bet20 + (1/sig22) * X_temp.t() * (Z_temp - L_temp * gam_temp.col(i) - aV);
    bet_new.row(i) = my_mvrnorm(V_new * m_new, V_new).t();
  }
  
  return bet_new;
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec bet10_update(mat Sig1_inv,mat bet1,mat S1_inv, vec Sm1, int ns){

  mat V_new_inv = S1_inv + double(ns)*Sig1_inv ;
  mat V_new = V_new_inv.i();
  vec m_new = Sm1 + Sig1_inv * sum(bet1,0).t() ;
  
  vec out = my_mvrnorm(V_new * m_new, V_new);

  return out;
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec bet20_update(mat Sig2_inv,mat bet2,mat S2_inv, vec Sm2, int ns){
  
  mat V_new_inv = S2_inv + double(ns)*Sig2_inv ;
  mat V_new = V_new_inv.i();
  vec m_new = Sm2 + Sig2_inv * sum(bet2,0).t() ;
  vec out = my_mvrnorm(V_new * m_new, V_new);
  
  return out;
  
}

// mat Sig10b_update(mat bet1,vec bet10,mat M_b1,double nu_b1 ,int ns){
//  return 0; 
// }

// mat Sig20b_update(mat bet2,vec bet20,mat M_b2,double nu_b2 ,int ns){
//   return 0; 
// }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat gam1_update(List Z1,List X, List L1,mat Sig_g1_inv,vec gam10,mat bet1, double a_11, 
                vec V1,double sig21, List LtL, int n_lags, int ns,int nt){
  mat gam_new(ns,n_lags);
  mat bet_temp = bet1.t();
  
  for(int i = 0; i < ns ; i++){
    mat LtL_temp = LtL[i];
    mat V_new_inv = Sig_g1_inv + (1/sig21) * LtL_temp ;
    mat V_new = V_new_inv.i();
    
    mat X_temp = X[i]; 
    mat L_temp = L1[i]; 
    vec aV =  a_11 * ones(nt)*V1[i];
    vec Z_temp = Z1[i];
    vec m_new = Sig_g1_inv * gam10 + (1/sig21) * L_temp.t() * (Z_temp - X_temp * bet_temp.col(i) - aV);
    gam_new.row(i) = my_mvrnorm(V_new * m_new, V_new).t();
  }
  
  return gam_new;
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat gam2_update(List Z2,List X, List L2,mat Sig_g2_inv,vec gam20,mat bet2, double a_12, 
                vec V1,vec V2,double sig22, List LtL, int n_lags, int ns,int nt){
  mat gam_new(ns,n_lags);
  mat bet_temp = bet2.t();
  
  for(int i = 0; i < ns ; i++){
    mat LtL_temp = LtL[i];
    mat V_new_inv = Sig_g2_inv + (1/sig22) * LtL_temp ;
    mat V_new = V_new_inv.i();
    
    mat X_temp = X[i]; 
    mat L_temp = L2[i]; 
    vec aV =  (a_12*V1[i] + V2[i]) * ones(nt);
    vec Z_temp = Z2[i];

    vec m_new = Sig_g2_inv * gam20 + (1/sig22) * L_temp.t() * (Z_temp - X_temp * bet_temp.col(i) - aV);
    gam_new.row(i) = my_mvrnorm(V_new * m_new, V_new).t();
  }
  
  return gam_new;
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec gam10_update(mat Sig1_inv,mat gam1,mat S1_inv, vec Sm1, int ns){
  
  mat V_new_inv = S1_inv + double(ns)*Sig1_inv ;
  mat V_new = V_new_inv.i();
  vec m_new = Sm1 + Sig1_inv * sum(gam1,0).t() ;
  vec out = my_mvrnorm(V_new * m_new, V_new);
  
  return out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec gam20_update(mat Sig2_inv,mat gam2,mat S2_inv, vec Sm2, int ns){
  
  mat V_new_inv = S2_inv + double(ns)*Sig2_inv ;
  mat V_new = V_new_inv.i();
  vec m_new = Sm2 + Sig2_inv * sum(gam2,0).t() ;
  vec out = my_mvrnorm(V_new * m_new, V_new);
  
  return out;
}


// mat Sig10g_update(mat gam1,vec gam10,mat M_b1,double nu_b1 ,int ns){
//   return 0; 
// }

// mat Sig20g_update(mat gam2,vec gam20,mat M_b2,double nu_b2 ,int ns){
//   return 0; 
//   
// }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec V1_update(List Z1,List Z2,List X, List L1, List L2,double a_11, double a_12, 
              vec V2,double sig21,double sig22,double tau21,mat bet1,mat bet2, 
              mat gam1,mat gam2,mat Q,int nt, int ns){
  
  mat V_new_inv = Q/tau21 + (double(nt) * pow(a_11,2.0) / sig21 + double(nt) * pow(a_12,2.0)/sig22) *eye<mat>(ns,ns)  ;
  mat V_new = V_new_inv.i();
  mat bet1_temp = bet1.t();
  mat bet2_temp = bet2.t();
  mat gam1_temp = gam1.t();
  mat gam2_temp = gam2.t();
  
  vec m_new(ns,fill::zeros);
  
  for(int i = 0 ; i < nt ; i++){
    mat X_temp = X[i]; 
    mat L1_temp = L1[i]; 
    mat L2_temp = L2[i]; 
    vec Z1_temp = Z1[i];
    vec Z2_temp = Z2[i];
    vec xb1(ns,fill::zeros);
    vec xb2(ns,fill::zeros);
    vec LG1(ns,fill::zeros);
    vec LG2(ns,fill::zeros);
    
    for(int j = 0; j < ns ; j ++){
      xb1[j] = sum(X_temp.row(j) * bet1_temp.col(j)) ;
      xb2[j] = sum(X_temp.row(j) * bet2_temp.col(j)) ;
      LG1[j] = sum(L1_temp.row(j) * gam1_temp.col(j)) ;
      LG2[j] = sum(L2_temp.row(j) * gam2_temp.col(j)) ;
    }
    
    m_new += (a_11/sig21) * (Z1_temp - xb1 - LG1) + (a_12/sig22)*(Z2_temp - xb2 - LG2 - V2);
  }
  
  vec out = my_mvrnorm(V_new * m_new, V_new);
  return out;
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec V2_update(List Z2,List X, List L2, double a_12, vec V1, double sig22,double tau22,
              mat bet2,mat gam2,mat Q,int nt, int ns){
  
  mat V_new_inv = Q/tau22  + ( double(nt) / sig22) * eye<mat>(ns,ns) ;
  mat V_new = V_new_inv.i();
  mat bet2_temp = bet2.t();
  mat gam2_temp = gam2.t();
  
  vec m_new(ns,fill::zeros);
  
  for(int i = 0 ; i < nt ; i++){
    mat X_temp = X[i]; 
    mat L2_temp = L2[i]; 
    vec Z2_temp = Z2[i];
    vec xb1(ns,fill::zeros);
    vec xb2(ns,fill::zeros);
    vec LG1(ns,fill::zeros);
    vec LG2(ns,fill::zeros);
    
    for(int j = 0; j < ns ; j ++){
      xb2[j] = sum(X_temp.row(j) * bet2_temp.col(j)) ;
      LG2[j] = sum(L2_temp.row(j) * gam2_temp.col(j)) ;
    }
    
    m_new +=  (1/sig22)*(Z2_temp - xb2 - LG2 - a_12 * V1);
  }
  
  vec out = my_mvrnorm(V_new * m_new, V_new);
  return out;
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double tau21_update(mat Q, vec V1,double a_t1,double b_t1,int ns){
  double a_new = a_t1 + (double(ns))/2 ;
  mat quad = 0.5 * V1.t() * (Q * V1);
  double b_new = b_t1 + quad(0,0);
  double out = 1/rgamma(1, a_new,1/b_new)[0];
  return out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double tau22_update(mat Q, vec V2,double a_t2,double b_t2,int ns){
  double a_new = a_t2 + (double(ns))/2 ;
  mat quad = 0.5 * V2.t() * (Q * V2);
  double b_new = b_t2 + quad(0,0);
  double out = 1/rgamma(1, a_new,1/b_new)[0];
  return out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double a12_update(List Z2,List X, List L2,mat bet2,mat gam2, vec V1,vec V2,
                  double sig22, double m, double s2, int ns,int nt){
  
  double m_new = m/s2;
  vec summed_V = sum(pow(V1,2.0) );
  
  double V_new = 1/(1/s2 + (double(nt))*summed_V[0]/sig22 );
  
  mat gam_temp = gam2.t();
  mat bet_temp = bet2.t();
  
  for(int i = 0; i < ns ; i++){
    mat X_temp = X[i]; 
    mat L_temp = L2[i]; 
    vec Z_temp = Z2[i];
    vec aV = ones(nt)*V2[i];
    
    m_new += sum(V1[i] * ( Z_temp - X_temp * bet_temp.col(i) - L_temp * gam_temp.col(i) - aV ) ) / sig22 ;
    
  }
  
  double out = rnorm(1,V_new * m_new , sqrt(V_new) )[0];
  
  return out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double Z1_update_alt( int t_ind, vec Z_temp,mat X_temp, mat L_temp,rowvec b1,rowvec g1, double a_11, 
                      double V_temp,double sig21, int n_lag,IntegerVector lags){
  
  int t_idx = t_ind - 1;
  double m_new = (sum(b1 * X_temp.row(t_idx).t() + g1 * L_temp.row(t_idx).t()) + V_temp) / sig21;
  double v_inv_new = 1 / sig21 ;
  
  if(n_lag != 0){
    for(int i = 0; i < n_lag ; i++){
      int t_idx_new = t_idx + lags[i] ;
      vec L_row = L_temp.row(t_idx_new).t();
      vec X_row = X_temp.row(t_idx_new).t();
      m_new += g1[i] * (Z_temp[t_idx_new] - sum(b1 * X_row + g1 * L_row - g1[i] * L_row[i] + V_temp)) / sig21;
      v_inv_new += pow(g1[i],2.0) / sig21 ;
    }
  }
  
  double v_new = 1 / v_inv_new;
  double z_out = rnorm(1,v_new * m_new, sqrt(v_new) )[0];
  
  return z_out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double Z2_update_alt(int t_ind, vec Z_temp,mat X_temp, mat L_temp,rowvec b1,rowvec g1, double a_12, 
                     double V1,double V2,double sig22, int n_lag,IntegerVector lags){
  
  int t_idx = t_ind - 1;
  double V_temp = V1 * a_12 + V2;
  double m_new = (sum(b1 * X_temp.row(t_idx).t() + g1 * L_temp.row(t_idx).t()) + V_temp) / sig22;
  double v_inv_new = 1 / sig22 ;
  
  if( n_lag != 0){
    for(int i = 0; i < n_lag ; i++){
      int t_idx_new = t_idx + lags[i] ;
      vec L_row = L_temp.row(t_idx_new).t();
      vec X_row = X_temp.row(t_idx_new).t();
      m_new += g1[i] * (Z_temp[t_idx_new] - sum(b1 * X_row + g1 * L_row - g1[i] * L_row[i] + V_temp)) / sig22;
      v_inv_new += pow(g1[i],2.0) / sig22 ;
    }  
  }
  
  
  double v_new = 1 / v_inv_new;
  double z_out = rnorm(1,v_new * m_new, sqrt(v_new) )[0];
  
  return z_out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double Z1_update_alt_fut( int t_ind, vec Z_temp,mat X_temp, mat L_temp,rowvec b1,rowvec g1, double a_11, 
                      double V_temp,double sig21, int n_lag,IntegerVector lags){
  
  int t_idx = t_ind - 1;
  double m_new = (sum(b1 * X_temp.row(t_idx).t() + g1 * L_temp.row(t_idx).t()) + V_temp) ;
  double z_out = rnorm(1,m_new, sqrt(sig21) )[0];
  
  return z_out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double Z2_update_alt_fut(int t_ind, vec Z_temp,mat X_temp, mat L_temp,rowvec b1,rowvec g1, double a_12, 
                     double V1,double V2,double sig22, int n_lag,IntegerVector lags){
  
  int t_idx = t_ind - 1;
  double V_temp = V1 * a_12 + V2;
  double m_new = (sum(b1 * X_temp.row(t_idx).t() + g1 * L_temp.row(t_idx).t()) + V_temp);
  double z_out = rnorm(1,m_new, sqrt(sig22) )[0];
  
  return z_out;
}

