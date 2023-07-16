#pragma once

#include <functional>
#include <blaze/Blaze.h>
#include <blaze/math/dense/DynamicVector.h>
#include <blaze/math/dense/StaticMatrix.h>
#include <blaze/math/dense/StaticVector.h>

namespace ew_cmaes {
using fitness_function_t = std::function<double(double*, int)>;
using number_t = double;
}
