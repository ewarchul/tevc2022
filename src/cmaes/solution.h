#pragma once

#include <blaze/math/adaptors/symmetricmatrix/BaseTemplate.h>
#include <blaze/math/dense/DynamicMatrix.h>
#include <blaze/math/dense/DynamicVector.h>
#include <blaze/math/dense/StaticMatrix.h>
#include <blaze/math/dense/StaticVector.h>
#include <blaze/math/expressions/DMatGenExpr.h>

#include "types.h"
#include "util.h"
namespace ew_cmaes {

template <int Dimension> struct solutions {
  solutions() : B(math::diag(dim_, 1)), D(math::diag(dim_, 1)) {
    BD = B * D;
    C = BD * blaze::trans(B*D);
  }

  unsigned dim_{Dimension};

  unsigned iter_{};
  unsigned fevals_{};
  double sigma_{1};
  bool hsig_{true};

  blaze::StaticVector<number_t, Dimension> mean_{};
  blaze::DynamicMatrix<number_t> population_{};
  blaze::StaticVector<number_t, Dimension> best_so_far_{};

  blaze::StaticVector<number_t, Dimension> pcov_{};
  blaze::StaticVector<number_t, Dimension> psigma_{};

  blaze::StaticMatrix<number_t, Dimension, Dimension> B{};
  blaze::StaticMatrix<number_t, Dimension, Dimension> D{};
  blaze::StaticMatrix<number_t, Dimension, Dimension> BD{};
  blaze::SymmetricMatrix<blaze::StaticMatrix<number_t, Dimension, Dimension>> C{};
};

}  // namespace ew_cmaes
