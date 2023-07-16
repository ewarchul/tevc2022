#pragma once

#include <blaze/math/dense/DynamicMatrix.h>
#include <blaze/math/dense/StaticMatrix.h>
#include <blaze/math/dense/StaticVector.h>
#include <blaze/math/expressions/DMatGenExpr.h>

#include "types.h"
#include "util.h"
namespace ew_cmaes {

template <int Dimension> struct solutions {
  solutions() : B(math::diag(dim_, dim_)), D(math::diag(dim_, dim_)) {
/*    auto BD_t = blaze::transpose(B * D);*/
    /*C = (B * D) * BD_t;*/
  }

  unsigned dim_{Dimension};

  unsigned iter_{0};
  unsigned fevals_{0};
  double sigma_{-1};
  bool hsig_{true};

  blaze::StaticVector<number_t, Dimension> mean_;
  blaze::DynamicMatrix<number_t> population_;
  blaze::StaticVector<number_t, Dimension> best_so_far_;

  blaze::StaticMatrix<number_t, Dimension, Dimension> B;
  blaze::StaticMatrix<number_t, Dimension, Dimension> D;
  blaze::StaticMatrix<number_t, Dimension, Dimension> C;
};

}  // namespace ew_cmaes
