#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// numerically stable version of log(inv_logit(x))
// [[Rcpp::export]]
NumericVector log_inv_logit_fast(const NumericVector & l) {
    NumericVector vec(l.size());

    for(std::size_t i=0; i < l.size(); ++i) {
      if(l[i] < 0.0) {
        vec[i] = l[i] - log1p(exp(l[i]));
      } else {
        vec[i] = -log1p(exp(-l[i]));
      }
    }
    
    return vec;
}

// numerically stable version of log1m_exp(x) = log(1-exp(x))
// [[Rcpp::export]]
NumericVector log1m_exp_max0_fast(const NumericVector & l) {
    NumericVector vec(l.size());

    for(std::size_t i=0; i < l.size(); ++i) {
      if(l[i] > 0.0) {
        vec[i] = NAN;
      } else if(l[i] > -0.693147) {
        vec[i] = log(-expm1(l[i]));
      } else {
        vec[i] = log1p(-exp(l[i]));
      }
    }
   
    return vec;
}
