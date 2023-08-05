#include <bits/chrono.h>
#include <blaze/math/StorageOrder.h>
#include <blaze/math/TransposeFlag.h>
#include <blaze/math/dense/DynamicVector.h>
#include <blaze/math/dense/Forward.h>
#include <blaze/math/expressions/Forward.h>
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
#include "cmaes/random.h"
#include "cmaes/step_size_update.h"
#include "cmaes/util.h"

#include <chrono>

auto sphere_fn(const auto& x) -> double {
  double result{0.0};
  for (int i{0}; i < x.size(); ++i) {
    result += x[i] * x[i];
  }
  return result;
}

auto main() -> int {
  constexpr int dim{10};
  blaze::StaticVector<double, dim> vec(95);

  auto params = ew_cmaes::parameters{dim};
  auto csa = ew_cmaes::csa_updater{params};
  auto fn = [](const auto& x) { return sphere_fn(x); };

  auto cma = ew_cmaes::cmaes(vec, fn, csa, params);
  
  auto start = std::chrono::system_clock::now();
  const auto ret = cma.run();
  auto end = std::chrono::system_clock::now();

  fmt::print("elapsed time -> {}ms\n", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    
  std::cout << ret.best_so_far_fitness_ << std::endl;

  return 0;
}
