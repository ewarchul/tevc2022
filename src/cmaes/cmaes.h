#pragma once

#include <blaze/math/dense/DynamicMatrix.h>
#include <blaze/math/dense/Eigen.h>
#include <blaze/math/expressions/DMatMapExpr.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/Matrix.h>
#include <blaze/math/simd/Sqrt.h>

#include <exception>
#include <utility>

#include "parameter.h"
#include "solution.h"
#include "step_size_update.h"
#include "termination.h"
#include "types.h"
#include "util.h"

namespace ew_cmaes {

template <int Dimension, typename StepSizeUpdater = dummy_step_size_updater> class cmaes {
 public:
  cmaes(dvec_t x0, fitness_function_t fit_fn, parameters<Dimension> params,
        termination_criteria<Dimension> stops = predefined::predefined_termination_criteria<Dimension>)
      : fn_{std::move(fit_fn)}, params_{std::move(params)}, stops_{std::move(stops)} {
    solutions_.mean_ = x0;
  }

  [[nodiscard]] auto run() -> solutions<Dimension> {
    while (!terminate()) {
      const auto [pop, diffs] = ask();
      auto fn_vals = evaluate(pop, fn_);
      solutions_.fevals_ += params_.lambda_;
      auto selected_indices = selection_sort(fn_vals, params_.mu_);
      auto selected_pop = blaze::columns(pop, selected_indices);
      auto selected_diffs = blaze::columns(diffs, selected_indices);
      tell(std::move(selected_pop), std::move(selected_diffs));

      std::cout << fn_vals << std::endl;
      fmt::print("[fevals = {}] Running cmaes solver...\n", solutions_.fevals_);
    }

    return solutions_;
  }

 private:
  [[nodiscard]] auto ask() {
    // ~ N(0, C)
    const dmat_t diffs = math::random(Dimension, params_.lambda_);
    // ~N(m, s^2C)
    dmat_t population = blaze::expand(solutions_.mean_, params_.lambda_) + solutions_.sigma_ * solutions_.BD * diffs;
    return std::make_pair(population, diffs);
  }

  auto tell(dmat_t&& selected_pop, dmat_t&& selected_diffs) -> void {
    // update mean-point
    solutions_.mean_ = selected_pop * params_.weights_;
    // update evolution paths
    const auto zmean = selected_diffs * params_.weights_;
    solutions_.psigma_ =
        params_.psigma_decay_factor_ * solutions_.psigma_ + params_.psigma_coeff_ * (solutions_.eigen_vecs * zmean);
    // calculate stall flag for pcov
    const auto hsig_exp = 2.0 * solutions_.fevals_ / params_.lambda_;
    solutions_.hsig_ = (blaze::norm(solutions_.psigma_) / blaze::sqrt(1.0 - blaze::pow(params_.psigma_decay_factor_, hsig_exp)) /
                        params_.chin_) < params_.hsig_coeff_;

    solutions_.pcov_ =
        params_.pcov_decay_factor_ * solutions_.pcov_ + solutions_.hsig_ * params_.pcov_coeff_ * (solutions_.BD * zmean);

    update_cma(std::move(selected_diffs));
    update_step_size();

    // eigendecomposition
    dvec_t eigen_vals;
    blaze::eigen(solutions_.cov_mat, eigen_vals, solutions_.eigen_vecs);
    const dmat_t eigen_diag = math::diag(eigen_vals, [](const auto& e) { return blaze::sqrt(e); });
    solutions_.BD = solutions_.eigen_vecs * eigen_diag;
  }

  auto update_step_size() -> void { solutions_.sigma_ = sigma_updater_(solutions_.sigma_); }

  auto update_cma(auto&& selz) -> void {
    const dmat_t BDz = solutions_.BD * selz;

    const smat_t<Dimension, Dimension> empirical_cov = (1.0 - params_.ccov_) * solutions_.cov_mat;
    const smat_t<Dimension, Dimension> rank_mu_cov = params_.ccov_ * (1.0 / params_.cmu_) *
                         (blaze::outer(solutions_.pcov_, solutions_.pcov_) +
                          (1.0 - solutions_.hsig_) * params_.cc_ * (2.0 - params_.cc_) * solutions_.cov_mat);
    const dmat_t rank_one_cov = params_.ccov_ * (1 - 1 / params_.cmu_) * BDz * params_.weights_diag_ * blaze::trans(BDz);

    solutions_.cov_mat = blaze::declsym(empirical_cov + rank_mu_cov + rank_one_cov);
  }

  [[nodiscard]] auto terminate() -> bool {
    const auto terminate_status = stops_.check(params_, solutions_);

    if (terminate_status.terminate_) {
      fmt::print("{}\n", terminate_status.msg_);
    }
    return terminate_status.terminate_;
  }

  fitness_function_t fn_;
  parameters<Dimension> params_;
  solutions<Dimension> solutions_{};
  termination_criteria<Dimension> stops_;

  StepSizeUpdater sigma_updater_;
};

}  // namespace ew_cmaes
