#include <Rcpp.h>
using namespace Rcpp;

//quantile implementation for default R method (type 7)
//' @export
// [[Rcpp::export]]
NumericVector fquantile(NumericVector x,
                        NumericVector probs,
                        bool na_rm = true) {
  if (all(is_na(x)))
  {
    NumericVector out(probs.size(), NA_REAL);
    return(out);
  }
  if (na_rm) x = na_omit(x);
  const int n = x.size();
  NumericVector out(probs.size());
  IntegerVector ii(probs.size());
  NumericVector h(probs.size());
  NumericVector index = 1 + (n - 1) * probs;
  NumericVector lo = floor(index); //floor
  //ceiling
  NumericVector hi = ceiling(index);
  for(int i = 0; i < probs.size(); i++)
  {//catch corner case when index element is int and ceiling = floor
    h[i] = index[i] - lo[i];
  }
  std::sort(x.begin(), x.end());
  out = x[as<IntegerVector>(lo) - 1];
  if (all(is_na(ii)))
  {
    return(out);
  } else
  {
    x = x[as<IntegerVector>(hi) - 1];
    for(int i = 0; i < probs.size(); i++)
    {
      out[i] = (1 - h[i]) * out[i] + h[i] * x[i];
    }
    return(out);
  }
}
