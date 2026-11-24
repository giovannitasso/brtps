// Spatial regTPS model (Cavieres & Krumscheid, 2026, Sect. 3.5 and App. C.2).
//
//   eta = X beta + offset + Phi z
//   z_k = z_tilde_k / sqrt(1 + alpha * v_k)      (non-centred, eq. A20)
//   z_tilde_k ~ N(0, 1)
//
// Phi (N x M) holds the truncated discrete eigenfunctions evaluated at the
// observed locations and v (length M) the matching eigenvalues of the TPS
// penalty matrix S, with v_k = 0 for the polynomial null space. Both are
// computed in R (regtps_basis()).

#ifndef BRTPS_REGTPS_S_HPP
#define BRTPS_REGTPS_S_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type regtps_s(objective_function<Type>* obj)
{
  // Data ---------------------------------------------------------------------
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_VECTOR(offset);
  DATA_VECTOR(size);          // binomial trials; unused by other families
  DATA_MATRIX(Phi);
  DATA_VECTOR(v);
  DATA_INTEGER(family);

  DATA_IVECTOR(prior_beta_code);
  DATA_MATRIX(prior_beta_par);    // one row per coefficient: p1, p2, p3
  DATA_INTEGER(prior_alpha_code);
  DATA_VECTOR(prior_alpha_par);   // p1, p2, p3
  DATA_IVECTOR(prior_phi_code);
  DATA_MATRIX(prior_phi_par);     // one row per dispersion parameter

  // Parameters ---------------------------------------------------------------
  PARAMETER_VECTOR(beta);
  PARAMETER(log_alpha);
  PARAMETER_VECTOR(log_phi);
  PARAMETER_VECTOR(z_tilde);

  // regTPS field -------------------------------------------------------------
  Type alpha = exp(log_alpha);
  vector<Type> z = z_tilde / sqrt(Type(1) + alpha * v);
  vector<Type> eta = X * beta + offset + Phi * z;

  // Negative log posterior ---------------------------------------------------
  Type nll = Type(0);
  nll -= brtps::family_loglik(y, eta, log_phi, size, family);
  nll -= dnorm(z_tilde, Type(0), Type(1), true).sum();

  for (int j = 0; j < beta.size(); j++) {
    nll -= brtps::log_prior_real(beta(j), prior_beta_code(j),
                                 prior_beta_par(j, 0), prior_beta_par(j, 1),
                                 prior_beta_par(j, 2));
  }
  nll -= brtps::log_prior_logpos(log_alpha, prior_alpha_code,
                                 prior_alpha_par(0), prior_alpha_par(1),
                                 prior_alpha_par(2));
  for (int j = 0; j < log_phi.size(); j++) {
    nll -= brtps::log_prior_logpos(log_phi(j), prior_phi_code(j),
                                   prior_phi_par(j, 0), prior_phi_par(j, 1),
                                   prior_phi_par(j, 2));
  }

  REPORT(z);
  REPORT(eta);
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
