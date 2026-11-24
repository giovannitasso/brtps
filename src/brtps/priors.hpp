// Prior densities shared by all brtps models.
//
// Every prior is identified by an integer code and up to three numeric
// parameters (p1, p2, p3). The codes must match the table in R/priors.R.

#ifndef BRTPS_PRIORS_HPP
#define BRTPS_PRIORS_HPP

namespace brtps {

// Log density of a real-valued parameter x (regression coefficients).
//   0 flat (improper, contributes nothing)
//   1 normal(mean = p1, sd = p2)
//   2 student_t(df = p1, location = p2, scale = p3)
template<class Type>
Type log_prior_real(Type x, int code, Type p1, Type p2, Type p3)
{
  switch (code) {
  case 0:
    return Type(0);
  case 1:
    return dnorm(x, p1, p2, true);
  case 2:
    return dt((x - p2) / p3, p1, true) - log(p3);
  default:
    Rf_error("brtps: unknown prior code %d for a real parameter", code);
  }
  return Type(0);
}

// Log density of theta = log(x) for a positive parameter x.
//
// The samplers work on theta, so the Jacobian |dx/dtheta| = x is included
// (the "+ theta" terms) whenever the prior is stated on x.
//   0 flat on theta (improper)
//   1 exponential(rate = p1); a PC prior for a standard deviation
//   2 half-normal(scale = p1)
//   3 half-Cauchy(scale = p1)
//   4 half-t(df = p1, scale = p2)
//   5 lognormal(meanlog = p1, sdlog = p2), i.e. theta ~ normal(p1, p2)
//   6 gamma(shape = p1, rate = p2)
template<class Type>
Type log_prior_logpos(Type theta, int code, Type p1, Type p2, Type p3)
{
  Type x = exp(theta);
  switch (code) {
  case 0:
    return Type(0);
  case 1:
    return log(p1) - p1 * x + theta;
  case 2:
    return log(Type(2)) + dnorm(x, Type(0), p1, true) + theta;
  case 3:
    return log(Type(2)) - log(Type(M_PI) * p1 * (Type(1) + (x / p1) * (x / p1))) + theta;
  case 4:
    return log(Type(2)) + dt(x / p2, p1, true) - log(p2) + theta;
  case 5:
    return dnorm(theta, p1, p2, true);
  case 6:
    return dgamma(x, p1, Type(1) / p2, true) + theta;
  default:
    Rf_error("brtps: unknown prior code %d for a positive parameter", code);
  }
  return Type(0);
}

} // namespace brtps

#endif
