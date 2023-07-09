#pragma once

#include "types.h"

namespace ew_cmaes {
struct parameters {
  //runtime parameters
  unsigned max_iter_;
  unsigned max_fevals_;
  double xtol_;
  double ftol_;
  double ftarget_;

  double lower_bound_;
  double upper_bound_;

  // es parameters
  unsigned dim_;
  unsigned lambda_;
  unsigned mu_;
  double sigma_;

  // cma parameters
  double mueff_;
  double cc_;
  double cmu_;
  double ccov_;
  vec weights_;
};
}  // namespace ew_cmaes
