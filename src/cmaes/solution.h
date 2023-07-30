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
  solutions() : eigen_vecs(math::diag(dim_, 1))  {
    smat_t<Dimension, Dimension> x = math::diag(dim_, 1);
    BD = eigen_vecs * x;
    cov_mat = BD * blaze::trans(eigen_vecs * x);
  }

  unsigned dim_{Dimension};

  unsigned iter_{};
  unsigned fevals_{};
  double sigma_{1};
  bool hsig_{true};

  svec_t<Dimension> mean_{};
  dmat_t population_{};
  svec_t<Dimension> best_so_far_{};

  svec_t<Dimension> pcov_{};
  svec_t<Dimension> psigma_{};

  smat_t<Dimension, Dimension> eigen_vecs{};
  smat_t<Dimension, Dimension> BD{};
  blaze::SymmetricMatrix<smat_t<Dimension, Dimension>> cov_mat{};
};

}  // namespace ew_cmaes
