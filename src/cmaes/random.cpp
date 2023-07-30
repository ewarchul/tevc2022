#include "random.h"

#include <random>

namespace {
std::random_device rd{};
std::mt19937 gen{rd()};
}  // namespace

namespace ew_cmaes::random {
auto rnorm() -> double {
  auto iso_normal_dist = std::normal_distribution<>{0, 1};
  return iso_normal_dist(gen);
}
}  // namespace ew_cmaes::random
