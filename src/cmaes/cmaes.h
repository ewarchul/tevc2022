#pragma once

#include <utility>

#include "parameter.h"
#include "solution.h"
#include "termination.h"
#include "types.h"

namespace ew_cmaes {

class cmaes {
 public:
  cmaes(parameters params,
        termination_criteria stops = predefined::predefined_termination_criteria)
      : params_{std::move(params)}, stops_{std::move(stops)} {
    solutions_.iter_ = 0;
  }

  [[nodiscard]] auto run() -> solutions {
    while (!terminate()) {
      fmt::print("[iter = {}] Running cmaes solver...\n", solutions_.iter_);
      ++solutions_.iter_;
    }

    return solutions_;
  }

 private:
  [[nodiscard]] auto ask() -> matrix;
  auto tell() -> void;
  [[nodiscard]] auto terminate() -> bool {
    const auto terminate_status = stops_.check(params_, solutions_);

    if (terminate_status.terminate_) {
      fmt::print("{}\n", terminate_status.msg_);
    }
    return terminate_status.terminate_;
  }

  //  fitness_function fn_;
  parameters params_;
  solutions solutions_{};
  termination_criteria stops_;
};

}  // namespace ew_cmaes
