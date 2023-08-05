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

struct parameters {
  parameters() = default;

  parameters(u32_t dimension) : dim_{dimension} {
    weights_ = blaze::log(mu_ + 1) - blaze::log(blaze::linspace(mu_, as<u32_t>(1), mu_));
    weights_ = weights_ / blaze::sum(weights_);
    weights_diag_ = math::diag(weights_);

    mueff_ = blaze::pow(blaze::sum(weights_), 2) / blaze::sum(blaze::pow(weights_, 2));
    cmu_ = mueff_;
    ccov_ = (1.0 / cmu_) * 2.0 / blaze::pow(dim_ + 1.4, 2) +
            (1.0 - 1.0 / cmu_) * ((2.0 * cmu_ - 1.0) / (blaze::pow((dim_ + 2), 2.0) + 2.0 * cmu_));

    csigma_ = (mueff_ + 2) / (dim_ + mueff_ + 3);
    chin_ = blaze::sqrt(dim_) * (1.0 - (1.0 / (4.0 * dim_)) + (1.0 / (21.0 * dim_ * dim_)));
    psigma_coeff_ = blaze::sqrt(csigma_ * (2.0 - csigma_) * mueff_);
    psigma_decay_factor_ = 1.0 - csigma_;
    pcov_decay_factor_ = 1.0 - cc_;
    pcov_coeff_ = blaze::sqrt(cc_ * (2.0 - cc_) * mueff_);
  }

  u32_t dim_;
  number_t hsig_coeff_{(1.4 + 2.0 / (dim_ + 1.0))};

  u64_t max_fevals_{consts::MAX_FEVALS_MUL * dim_};
  u32_t lambda_{consts::LAMBDA_MUL * dim_};
  u32_t mu_{as<u32_t>(blaze::floor(as<number_t>(lambda_) / 2))};

  number_t xtol_{consts::DEFAULT_TOL};
  number_t lower_bound_{consts::DEFAULT_LOWER_BOUND};
  number_t upper_bound_{consts::DEFAULT_UPPER_BOUND};

  number_t sigma_{consts::DEFAULT_STEP_SIZE};
  number_t cc_{4.0 / (dim_ + 4.0 + 0.0)};
  number_t csigma_{};

  number_t psigma_decay_factor_{};
  number_t psigma_coeff_{};

  number_t pcov_decay_factor_{};
  number_t pcov_coeff_{};

  number_t ccov_{};
  number_t chin_{};
  number_t cmu_{};
  number_t mueff_{};
  dvec_t weights_{};
  dmat_t weights_diag_{};
};
}  // namespace ew_cmaes
