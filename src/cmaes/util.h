#pragma once
#include <blaze/math/TransposeFlag.h>
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

#include "concepts.h"
#include "random.h"
#include "types.h"

namespace ew_cmaes {

constexpr auto to_underlying(auto input) -> std::underlying_type_t<decltype(input)> {
  return static_cast<std::underlying_type_t<decltype(input)>>(input);
}

template <typename OutputType> constexpr auto as(auto input) -> OutputType { return static_cast<OutputType>(input); }

namespace math {
inline auto diag(auto size, Numeric auto scalar) {
  return blaze::generate(size, size, [scalar](auto row, auto col) { return row == col ? scalar : 0; });
}

inline auto diag(const ranges::range auto& vector) {
  return blaze::generate(vector.size(), vector.size(), [&vector](auto row, auto col) { return row == col ? vector[row] : 0; });
}

inline auto diag(const ranges::range auto& vector, auto&& op) {
  return blaze::generate(vector.size(), vector.size(), [&](auto row, auto col) { return row == col ? op(vector[row]) : 0; });
}

inline auto random(auto rows, auto cols) {
  return blaze::generate(rows, cols, [](auto r, auto c) { return random::rnorm(); });
}

}  // namespace math

auto selection_sort(auto&& fitness_vals, const auto selected_num) {
  auto zipped = fitness_vals | ranges::views::enumerate | ranges::to_vector;
  ranges::sort(zipped, [](const auto& r1, const auto& r2) { return r1.second < r2.second; });
  auto indices = zipped | ranges::views::take(selected_num) | ranges::views::keys | ranges::to<std::vector<int>>;
  auto fitness = zipped | ranges::views::take(selected_num) | ranges::views::values | ranges::to<std::vector<int>>;
  return std::make_pair(indices, fitness);
}

auto evaluate(const auto& population, auto&& fitness_fn) {
  blaze::DynamicVector<double, blaze::rowVector> fitness_values(population.columns());
  for (auto column : ranges::views::iota(0) | ranges::views::take(population.columns())) {
    const auto& col_view = blaze::column(population, column);
    fitness_values[column] = fitness_fn(col_view);
  }
  return fitness_values;
}

}  // namespace ew_cmaes
