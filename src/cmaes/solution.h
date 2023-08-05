#pragma once

#include <blaze/math/adaptors/symmetricmatrix/BaseTemplate.h>
#include <blaze/math/dense/DynamicMatrix.h>
#include <blaze/math/dense/DynamicVector.h>
#include <blaze/math/dense/StaticMatrix.h>
#include <blaze/math/dense/StaticVector.h>
#include <blaze/math/expressions/DMatGenExpr.h>
#include <blaze/math/sparse/ZeroVector.h>

#include <limits>

#include "types.h"
#include "util.h"
namespace ew_cmaes {

struct solutions {
  solutions() = default;

  solutions(u32_t dimension) : dim_{dimension} {
    psigma_ = blaze::zero<number_t>(dim_);
    pcov_ = blaze::zero<number_t>(dim_);

    eigen_vecs = math::diag(dim_, 1);
    BD = eigen_vecs * math::diag(dim_, 1);
    cov_mat = BD * blaze::trans(BD);
  }

  auto inc_feval(const u32_t lambda) { fevals_ += lambda; }
  auto inc_iter() { ++iter_; }

  auto update_best(const auto& current_best, const auto fitness) {
    if (fitness < best_so_far_fitness_) {
      best_so_far_ = current_best;
      best_so_far_fitness_ = fitness;
    }
    best_log_.push_back(current_best);
    best_log_fitness.push_back(fitness);
  }

  u32_t dim_{};
  unsigned iter_{};
  unsigned fevals_{};

  double sigma_{1};

  dvec_t mean_{};
  dmat_t population_{};

  dvec_t best_so_far_{};
  number_t best_so_far_fitness_{std::numeric_limits<number_t>{}.max()};
  std::vector<dvec_t> best_log_{};
  std::vector<number_t> best_log_fitness{};

  bool hsig_{true};
  dvec_t pcov_;
  dvec_t psigma_;
  dmat_t eigen_vecs;
  dmat_t BD;
  blaze::SymmetricMatrix<dmat_t> cov_mat{};
};

}  // namespace ew_cmaes
