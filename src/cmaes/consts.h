#pragma once

#include "types.h"


namespace ew_cmaes::consts {
  static constexpr unsigned MAX_ITER_MUL{10000};
  static constexpr unsigned MAX_FEVALS_MUL{1000};
  static constexpr unsigned LAMBDA_MUL{4};
  static constexpr number_t DEFAULT_TOL{10e-6};
  static constexpr number_t DEFAULT_LOWER_BOUND{-100.0};
  static constexpr number_t DEFAULT_UPPER_BOUND{100.0};
  static constexpr number_t DEFAULT_STEP_SIZE{1};
}
