#pragma once
#include <blaze/math/expressions/DMatGenExpr.h>
#include <blaze/util/Exception.h>

#include <iostream>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/iterator/operations.hpp>
#include <range/v3/range/concepts.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/take.hpp>
#include <span>
#include <type_traits>

#include "random.h"
#include "types.h"

namespace ew_cmaes {

constexpr auto to_underlying(auto input)
    -> std::underlying_type_t<decltype(input)> {
  return static_cast<std::underlying_type_t<decltype(input)>>(input);
}

template <typename OutputType> constexpr auto as(auto input) -> OutputType {
  return static_cast<OutputType>(input);
}

namespace math {

inline auto diag(auto size, auto scalar) {
  return blaze::generate(size, size, [scalar](auto row, auto col) {
    return row == col ? scalar : 0;
  });
}

inline auto diag(const ranges::range auto& vector) {
  return blaze::generate(
      vector.size(), vector.size(),
      [vector](auto row, auto col) { return row == col ? vector[row] : 0; });
}

inline auto random(auto rows, auto cols) {
  return blaze::generate(rows, cols,
                         [](auto r, auto c) { return random::rnorm(); });
}

}  // namespace math

auto selection_sort(auto&& fitness_vals, const auto selected_num)
    -> std::vector<int> {
  auto zipped = fitness_vals | ranges::views::enumerate | ranges::to_vector;
  ranges::sort(
      zipped,
      [](const auto& r1, const auto& r2) { return r1.second > r2.second; });
  return zipped | ranges::views::take(selected_num) | ranges::views::keys | ranges::to<std::vector<int>>;
}

auto evaluate(const auto& population, auto&& fitness_fn) {
  blaze::DynamicVector<double> fitness_values(population.columns());
  for (auto column :
       ranges::views::iota(0) | ranges::views::take(population.columns())) {
    const auto& col_view = blaze::column(population, column);
    fitness_values[column] = fitness_fn(col_view);
  }
  return fitness_values;
}

}  // namespace ew_cmaes
