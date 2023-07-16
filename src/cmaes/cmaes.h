#pragma once

#include <blaze/math/dense/DynamicMatrix.h>
#include <blaze/math/expressions/Matrix.h>
#include <utility>

#include "parameter.h"
#include "solution.h"
#include "termination.h"
#include "types.h"
#include "util.h"

namespace ew_cmaes {

template <int Dimension> class cmaes {
 public:
  cmaes(parameters<Dimension> params,
        termination_criteria<Dimension> stops =
            predefined::predefined_termination_criteria<Dimension>)
      : params_{std::move(params)}, stops_{std::move(stops)} {
    solutions_.iter_ = 0;
  }

  [[nodiscard]] auto run() -> solutions<Dimension> {
    while (!terminate()) {
      
      const auto population = ask();
      const auto fn_vals = evaluate(population);
      tell();

      fmt::print("[fevals = {}] Running cmaes solver...\n", solutions_.fevals_);
      solutions_.fevals_ += params_.lambda_;
    }

    return solutions_;
  }

 private:
  [[nodiscard]] auto ask() -> blaze::DynamicMatrix<number_t> {
    const auto arz = math::random(Dimension, params_.lambda_);
    return solutions_.mean_ + solutions_.sigma_ * solutions_.B_mat * solutions_.D * arz;
  }

  auto tell() -> void;
  [[nodiscard]] auto terminate() -> bool {
    const auto terminate_status = stops_.check(params_, solutions_);

    if (terminate_status.terminate_) {
      fmt::print("{}\n", terminate_status.msg_);
    }
    return terminate_status.terminate_;
  }

  //  fitness_function fn_;
  parameters<Dimension> params_;
  solutions<Dimension> solutions_{};
  termination_criteria<Dimension> stops_;
};

}  // namespace ew_cmaes
