#pragma once

#include <functional>

namespace ew_cmaes {
using fitness_function = std::function<double(double*, int)>;

using vec = std::vector<std::vector<double>>;
using matrix = std::vector<std::vector<double>>;

}
