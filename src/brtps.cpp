// brtps: one compiled TMB library that holds every brtps model.
//
// Each model will live in its own header under src/brtps/ and is selected at
// run time through the DATA_STRING `model`, which the R side fills in.

#define TMB_LIB_INIT R_init_brtps
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(model);
  Rf_error("brtps: unknown model '%s'", model.c_str());
  return Type(0);
}
