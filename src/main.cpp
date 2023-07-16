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

auto sphere_fn(const auto& x, auto size) -> double {
  double result{0.0};
  for (int i{0}; i < size; ++i) {
    result += x[i] * x[i];
  }
  return result;
}

auto main() -> int {
  constexpr int dim{5};
  const auto xx = ew_cmaes::math::diag(std::vector{2, 5, 1, 3, 4});

  std::cout << xx << std::endl;

  blaze::DynamicVector<double, blaze::rowVector> r(dim);
  std::cout << xx.columns() << std::endl;
  for (auto col = 0; col < xx.columns(); ++col) {
    auto col_view = blaze::column(xx, col);
    r[col] = sphere_fn(col_view, col_view.size());
  }

  std::cout << r << std::endl;

  auto enumrated = r | ranges::views::enumerate | ranges::to_vector;

  ranges::sort(enumrated, [](const auto& r1, const auto& r2) {
    return r1.second > r2.second;
  });

  auto indicies = enumrated | ranges::views::take(2) | ranges::views::keys | ranges::to_vector;

  auto sorted = blaze::columns(xx, indicies);

  std::cout << sorted << std::endl;

  /*  auto params = ew_cmaes::parameters<dim>{};*/
  /*auto cma = ew_cmaes::cmaes<dim>(params);*/
  /*const auto ret = cma.run();*/

  return 0;
}
