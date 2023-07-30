#pragma once

#include <blaze/math/dense/DynamicVector.h>
#include <blaze/math/expressions/DVecGenExpr.h>
#include <blaze/math/expressions/DVecMapExpr.h>
#include <blaze/math/simd/Log.h>
#include <blaze/math/simd/Pow.h>

#include <range/v3/all.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/iota.hpp>

#include "consts.h"
#include "types.h"
#include "util.h"

namespace ew_cmaes {

template <int Dimension> struct parameters {
  static constexpr unsigned dim_ = Dimension;

  parameters() {
    weights_ = blaze::log(mu_ + 1) -
               blaze::log(blaze::linspace(mu_, as<unsigned>(1), mu_));

    weights_ = weights_ / blaze::sum(weights_);

    mueff_ = blaze::pow(blaze::sum(weights_), 2) /
             blaze::sum(blaze::pow(weights_, 2));
    cmu_ = mueff_;

    ccov_ = (1.0 / cmu_) * 2.0 / blaze::pow(dim_ + 1.4, 2) +
            (1.0 - 1.0 / cmu_) *
                ((2.0 * cmu_ - 1.0) / (blaze::pow((dim_ + 2), 2.0) + 2.0 * cmu_));

    csigma_ = (mueff_ + 2) / (dim_ + mueff_ + 3);

   chin_ = blaze::sqrt(dim_) * (1.0 - (1.0/(4.0 * dim_ )) + (1.0 / ( 21.0 * dim_ * dim_)));

  }

  // runtime parameters
  unsigned max_fevals_{consts::MAX_FEVALS_MUL * dim_};
  unsigned lambda_{consts::LAMBDA_MUL * dim_};
  unsigned mu_{as<unsigned>(blaze::floor(as<number_t>(lambda_) / 2))};
  double xtol_{consts::DEFAULT_TOL};
  double lower_bound_{consts::DEFAULT_LOWER_BOUND};
  double upper_bound_{consts::DEFAULT_UPPER_BOUND};
  // es parameters
  double sigma_{consts::DEFAULT_STEP_SIZE};
  // cma parameters
  double cc_{4.0 / (dim_ + 4.0 + 0.0)};
  double csigma_{};

  double ccov_{};
  double chin_{};
  double cmu_{};
  double mueff_{};
  blaze::DynamicVector<number_t> weights_;
};
}  // namespace ew_cmaes
