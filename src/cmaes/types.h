#pragma once

#include <functional>
#include <range/v3/range/concepts.hpp>
#include <span>
#include <blaze/Blaze.h>
#include <blaze/math/TransposeFlag.h>
#include <blaze/math/dense/DynamicMatrix.h>
#include <blaze/math/dense/DynamicVector.h>
#include <blaze/math/dense/StaticMatrix.h>
#include <blaze/math/dense/StaticVector.h>
#include <blaze/math/expressions/Vector.h>

namespace ew_cmaes {
using number_t = double;

using fitness_function_t = std::function<double(blaze::DynamicVector<double>)>;

}
