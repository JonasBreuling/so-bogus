/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_EIGEN_MATRIX_TRAITS_HPP
#define BOGUS_EIGEN_MATRIX_TRAITS_HPP

#include <Eigen/Core>

#include "EigenLinearSolvers.hpp"
#include "EigenSparseLinearSolvers.hpp"

#include "../Block/ScalarBindings.hpp"

namespace bogus {

template <typename _MatrixType>
struct MatrixTraits {
  using MatrixType = _MatrixType;
  using Scalar = typename MatrixType::Scalar;

  using LUType = LU<Eigen::MatrixBase<MatrixType> >;
  using LDLTType = LDLT<Eigen::MatrixBase<MatrixType> >;

  static const MatrixType& asConstMatrix(const MatrixType& src) { return src; }
};

template <typename _Scalar, int _Options, typename _Index>
struct MatrixTraits<Eigen::SparseMatrix<_Scalar, _Options, _Index> > {
  using Scalar = _Scalar;
  using MatrixType = Eigen::SparseMatrix<Scalar, _Options, _Index>;

  using LUType = LU<Eigen::SparseMatrixBase<Eigen::SparseMatrix<Scalar, _Options, _Index> > >;
  using LDLTType = LDLT<Eigen::SparseMatrixBase<Eigen::SparseMatrix<Scalar, _Options, _Index> > >;

  static const MatrixType& asConstMatrix(const MatrixType& src) { return src; }
};

#define BOGUS_PROCESS_SCALAR(Scalar)                                                \
  template <>                                                                       \
  struct MatrixTraits<Scalar> : public MatrixTraits<Eigen::Matrix<Scalar, 1, 1> > { \
    static Eigen::Matrix<Scalar, 1, 1> asConstMatrix(const Scalar src) {            \
      Eigen::Matrix<Scalar, 1, 1> mat;                                              \
      mat(0, 0) = src;                                                              \
      return mat;                                                                   \
    }                                                                               \
  };
BOGUS_BLOCK_SCALAR_TYPES
#undef BOGUS_PROCESS_SCALAR

}  // namespace bogus

#endif
