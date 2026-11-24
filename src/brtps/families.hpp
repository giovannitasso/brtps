// Observation models (likelihoods) shared by all brtps models.
//
// family_loglik() returns the log-likelihood of the whole response vector
// given the linear predictor eta and the log of the family's dispersion
// parameters (log_phi). The integer codes must match R/families.R.
//
// To add a family: add a code to the enum, a case to the switch and the
// matching entry in R/families.R.

#ifndef BRTPS_FAMILIES_HPP
#define BRTPS_FAMILIES_HPP

namespace brtps {

enum family_code {
  gaussian_identity = 0
  // poisson_log = 1, binomial_logit = 2, gamma_log = 3, lognormal = 4
};

template<class Type>
Type family_loglik(vector<Type> y, vector<Type> eta, vector<Type> log_phi,
                   vector<Type> size, int family)
{
  Type ll = Type(0);
  switch (family) {
  case gaussian_identity: {
    Type sigma = exp(log_phi(0));
    ll = dnorm(y, eta, sigma, true).sum();
    break;
  }
  default:
    Rf_error("brtps: family code %d is not implemented", family);
  }
  return ll;
}

} // namespace brtps

#endif
