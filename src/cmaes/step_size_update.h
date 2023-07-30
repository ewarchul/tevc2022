#pragma once

namespace ew_cmaes {

struct dummy_step_size_updater {
  auto operator()(double step_size) -> double { return step_size * param_; }

  double param_{1};
};

}  // namespace ew_cmaes
