#pragma once

#include <concepts>

namespace ew_cmaes {
template <typename T>
concept Numeric = std::integral<T> || std::floating_point<T>;
}
