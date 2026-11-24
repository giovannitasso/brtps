// brtps: one compiled TMB library that holds every brtps model.
//
// Each model lives in its own header under src/brtps/ and is selected at run
// time through the DATA_STRING `model`, which the R side fills in. To add a
// model (for example the spatio-temporal regtps_st):
//   1. write src/brtps/<model>.hpp following regtps_s.hpp,
//   2. include it below,
//   3. add one branch to the dispatcher.
// Likelihoods and priors are shared by all models (families.hpp, priors.hpp).

#define TMB_LIB_INIT R_init_brtps
#include <TMB.hpp>

#include "brtps/priors.hpp"
#include "brtps/families.hpp"
#include "brtps/regtps_s.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(model);

  if (model == "regtps_s") {
    return regtps_s(this);
  }

  Rf_error("brtps: unknown model '%s'", model.c_str());
  return Type(0);
}
