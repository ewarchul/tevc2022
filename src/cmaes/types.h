#pragma once

#include <blaze/Blaze.h>
#include <blaze/math/TransposeFlag.h>
#include <blaze/math/dense/DynamicMatrix.h>
#include <blaze/math/dense/DynamicVector.h>
#include <blaze/math/dense/StaticMatrix.h>
#include <blaze/math/dense/StaticVector.h>
#include <blaze/math/expressions/Vector.h>

#include <functional>
#include <range/v3/range/concepts.hpp>
#include <span>

namespace ew_cmaes {
using number_t = double;

using u8_t = uint8_t;
using u16_t = uint16_t;
using u32_t = uint32_t;
using u64_t = uint64_t;

using i8_t = int8_t;
using i16_t = int16_t;
using i32_t = int32_t;
using i64_t = int64_t;

using dvec_t = blaze::DynamicVector<number_t>;
using dmat_t = blaze::DynamicMatrix<number_t>;

template <std::size_t Length> using svec_t = blaze::StaticVector<number_t, Length>;

template <std::size_t Row, std::size_t Col> using smat_t = blaze::StaticMatrix<number_t, Row, Col>;

using fitness_function_t = std::function<number_t(dvec_t)>;

}  // namespace ew_cmaes
