#include <Rcpp.h>
using namespace Rcpp;
#include <stdlib.h>            // malloc
#include <stdio.h>             // printf
#include <math.h>              // fabs, sqrt, etc.
#include <time.h>              // time
#include <unistd.h>            // getpid
//#include <gsl/gsl_rng.h>    

NumericMatrix euc_dist(NumericMatrix x){
  int n = x.nrow();
  double d;
  NumericMatrix out(n,n);
  
  for (int i = 0; i < n - 1; i++){
    NumericVector v1 = x.row(i);
    for (int j = i + 1; j < n ; j ++){
      d = sqrt(sum(pow(v1-x.row(j), 2.0)));
      out(j,i)=d;
      out(i,j)=d;
    }
  }
  
  return out;
}

// [[Rcpp::export]]
double pairmean(NumericVector x) {
  int n = x.size();
  double total = 0;

  for(int i = 0; i < n ; i++) {
    for(int j = 0; j < n; j++) {
      total += fabs(x[i] - x[j]);
    }
  }
  double numb = double(n)*double(n);
  double avg = double(total) / numb;
  return avg;
}

// [[Rcpp::export]]
double pairnorm(NumericMatrix x) {
  return mean(euc_dist(x)) / 2.0 ; 
}
