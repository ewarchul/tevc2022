#pragma once

#include <blaze/math/simd/Exp.h>

#include "parameter.h"
namespace ew_cmaes {

struct dummy_updater {
  auto operator()(double step_size) -> double { return step_size * param_; }

  double param_{1};
};

struct csa_updater {
  csa_updater(const parameters& params) : chin_{params.chin_}, csigma_{params.csigma_} {
    damping_factor_ = 1 + 2 * blaze::max(0, blaze::sqrt((params.mueff_ - 1) / (params.dim_ + 1)) - 1) + csigma_;
  }

  auto operator()(double sigma, const auto& solutions) -> double {
    return sigma * blaze::exp(blaze::norm(solutions.psigma_) / chin_ - 1) * csigma_ / damping_factor_;
  }

  double chin_;
  double csigma_;
  double damping_factor_{};
};

struct ppmf_updater {};

struct msr_updater {};

struct tpa_updater {};

}  // namespace ew_cmaes
