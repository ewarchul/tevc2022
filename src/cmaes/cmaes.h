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

auto print(auto&& x) -> void {
  std::cout << x << std::endl;
}

namespace ew_cmaes {

template <int Dimension, typename StepSizeUpdater = dummy_step_size_updater>
class cmaes {
 public:
  cmaes(blaze::DynamicVector<double> x0, fitness_function_t fit_fn,
        parameters<Dimension> params,
        termination_criteria<Dimension> stops =
            predefined::predefined_termination_criteria<Dimension>)
      : fn_{std::move(fit_fn)},
        params_{std::move(params)},
        stops_{std::move(stops)} {
    solutions_.mean_ = x0;
  }

  [[nodiscard]] auto run() -> solutions<Dimension> {
    while (!terminate()) {
      const auto [arx, arz] = ask();
      auto fn_vals = evaluate(arx, fn_);
      solutions_.fevals_ += params_.lambda_;
      auto selected_indices = selection_sort(fn_vals, params_.mu_);
      auto selx = blaze::columns(arx, selected_indices);
      auto selz = blaze::columns(arz, selected_indices);
      tell(selx, selz);

      fmt::print("[fevals = {}] Running cmaes solver...\n", solutions_.fevals_);
    }

    return solutions_;
  }

 private:
  [[nodiscard]] auto ask() {
    const blaze::DynamicMatrix<number_t> arz =
        math::random(Dimension, params_.lambda_);
    blaze::DynamicMatrix<number_t> arx =
        blaze::expand(solutions_.mean_, params_.lambda_) +
        solutions_.sigma_ * solutions_.BD * arz;
    return std::make_pair(arx, arz);
  }

  auto tell(blaze::DynamicMatrix<number_t>&& selx,
            blaze::DynamicMatrix<number_t>&& selz) -> void {
    solutions_.mean_ = selx * params_.weights_;
    auto zmean = selz * params_.weights_;

    solutions_.psigma_ =
        (1.0 - params_.csigma_) * solutions_.psigma_ +
        blaze::sqrt(params_.csigma_ * (2.0 - params_.csigma_) * params_.mueff_) *
            (solutions_.B * zmean);
    

    solutions_.hsig_ = (blaze::norm(solutions_.psigma_) /
                           blaze::sqrt(1.0 - blaze::pow((1.0 - params_.csigma_),
                                                      2.0 * solutions_.fevals_ /
                                                          params_.lambda_)) /
                                               params_.chin_) <
                       (1.4 + 2.0 / (params_.dim_ + 1.0));

    solutions_.pcov_ =
        (1.0 - params_.cc_) * solutions_.pcov_ +
        solutions_.hsig_ *
            blaze::sqrt(params_.cc_ * (2.0 - params_.cc_) * params_.mueff_) *
            (solutions_.BD * zmean);

    update_cma(selz);
    update_step_size();
  }

  auto update_step_size() -> void {
    solutions_.sigma_ = sigma_updater_(solutions_.sigma_);
  }

  auto update_cma(auto&& selz) -> void {
    blaze::DynamicMatrix<double> BDz = solutions_.BD * selz;
    const auto BDz_t = blaze::trans(BDz);
    const auto C_1 = (1.0 - params_.ccov_) * solutions_.C;
    const auto C_2 = params_.ccov_ * (1.0 / params_.cmu_) *
                     (blaze::outer(solutions_.pcov_, solutions_.pcov_) +
                      (1.0 - solutions_.hsig_) * params_.cc_ * (2.0 - params_.cc_) *
                          solutions_.C);
    const blaze::DynamicMatrix<double> C_3 =
        params_.ccov_ * (1 - 1 / params_.cmu_) * BDz *
        math::diag(params_.weights_) * BDz_t;
    
    solutions_.C = blaze::declsym(C_1 + C_2 + C_3);

    blaze::DynamicVector<double> w;  // The vector for the real eigenvalues
    blaze::eigen(solutions_.C, w, solutions_.B);
    const blaze::DynamicMatrix<double> w_mat =
        math::diag(w, [](const auto& e) { return blaze::sqrt(e); });
    solutions_.BD = solutions_.B * w_mat;
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
