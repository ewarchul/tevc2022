#include <blaze/math/TransposeFlag.h>
#include <blaze/math/dense/DynamicVector.h>
#include <blaze/math/dense/Forward.h>
#include <blaze/util/Random.h>
#include <fmt/core.h>

#include <array>
#include <iostream>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/iterator/operations.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/map.hpp>

#include "cmaes/cmaes.h"
#include "cmaes/parameter.h"
#include "cmaes/util.h"

auto sphere_fn(const auto& x) -> double {
  double result{0.0};
  for (int i{0}; i < x.size(); ++i) {
    result += x[i] * x[i];
  }
  return result;
}

auto main() -> int {
  constexpr int dim{5};
  const auto xx = ew_cmaes::math::diag(std::vector{2, 5, 1, 3, 4});

  std::cout << xx << std::endl;

  auto fit_vals =
      ew_cmaes::evaluate(xx, [](const auto& x) { return sphere_fn(x); });

  std::cout << fit_vals << std::endl;

  auto selected_indices = ew_cmaes::selection_sort(std::move(fit_vals), 2);

  std::cout << ranges::views::all(selected_indices) << std::endl;

  auto params = ew_cmaes::parameters<dim>{};

  auto fn = [](const auto& x) { return sphere_fn(x); };

  auto cma = ew_cmaes::cmaes<dim>(fn, params);

  const auto ret = cma.run();

  return 0;
}
