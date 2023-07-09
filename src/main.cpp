#include <fmt/core.h>

#include "cmaes/cmaes.h"
#include "cmaes/parameter.h"

auto main() -> int {
  auto params = ew_cmaes::parameters{};
  params.max_iter_ = 100;
  auto cma = ew_cmaes::cmaes(params);

  const auto ret = cma.run();

  return 0;
}
