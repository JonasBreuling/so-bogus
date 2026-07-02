/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_EIGEN_PROBLEM_TRAITS_HPP
#define BOGUS_EIGEN_PROBLEM_TRAITS_HPP

#include "EigenMatrixTraits.hpp"

namespace bogus {

template <typename Scalar_>
struct ProblemTraits {
  using Scalar = Scalar_;

  using DynVector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using DynMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  template <typename OtherMatrix>
  struct MutableClone {
    using Type = Eigen::Matrix<Scalar, OtherMatrix::RowsAtCompileTime, OtherMatrix::ColsAtCompileTime>;
  };
};

template <DenseIndexType Dimension, typename Scalar_>
struct LocalProblemTraits : public ProblemTraits<Scalar_> {
  enum { dimension = Dimension };

  using Scalar = Scalar_;
  using Vector = Eigen::Matrix<Scalar, Dimension, 1>;
  using Array = Eigen::Array<Scalar, Dimension, 1>;
  using Matrix = Eigen::Matrix<Scalar, Dimension, Dimension>;

  using LUType = typename MatrixTraits<Matrix>::LUType;
  using LDLTType = typename MatrixTraits<Matrix>::LDLTType;

  using TgMatrix = Eigen::Matrix<Scalar, Dimension - 1, Dimension - 1>;

  static Scalar np(const Vector& v) { return v[0]; }
  static Scalar& np(Vector& v) { return v[0]; }

  static typename Vector::template ConstFixedSegmentReturnType<Dimension - 1>::Type tp(const Vector& v) {
    return v.template segment<Dimension - 1>(1);
  }
  static typename Vector::template FixedSegmentReturnType<Dimension - 1>::Type tp(Vector& v) {
    return v.template segment<Dimension - 1>(1);
  }

  static typename Matrix::ColXpr nc(Matrix& m) { return m.col(0); }
  static Eigen::Block<Matrix, Dimension, Dimension - 1> tc(Matrix& m) {
    return m.template block<Dimension, Dimension - 1>(0, 1);
  }

  static Scalar& nnb(Matrix& m) { return m(0, 0); }
  static Eigen::Block<Matrix, Dimension - 1, Dimension - 1> ttb(Matrix& m) {
    return m.template block<Dimension - 1, Dimension - 1>(1, 1);
  }
  static Eigen::Block<Matrix, Dimension - 1, 1> tnb(Matrix& m) { return m.template block<Dimension - 1, 1>(1, 0); }
  static Eigen::Block<Matrix, 1, Dimension - 1> ntb(Matrix& m) { return m.template block<1, Dimension - 1>(0, 1); }
};

template <typename Scalar>
struct LocalProblemTraits<Eigen::Dynamic, Scalar>
    : public ProblemTraits<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > {
  using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using Base = ProblemTraits<Matrix>;

  enum { dimension = internal::DYNAMIC };

  using Vector = typename Base::DynVector;
  using TgMatrix = typename Base::DynMatrix;

  static Scalar np(const Vector& v) { return v[0]; }
  static Scalar& np(Vector& v) { return v[0]; }

  static typename Vector::ConstSegmentReturnType tp(const Vector& v) { return v.segment(1, v.rows() - 1); }
  static typename Vector::SegmentReturnType tp(Vector& v) { return v.segment(1, v.rows() - 1); }

  static typename Matrix::ColXpr nc(Matrix& m) { return m.col(0); }
  static Eigen::Block<Matrix> tc(Matrix& m) { return m.block(0, 1, m.rows(), m.cols() - 1); }

  static Scalar& nnb(Matrix& m) { return m(0, 0); }
  static Eigen::Block<Matrix> ttb(Matrix& m) { return m.block(1, 1, m.rows() - 1, m.cols() - 1); }
  static Eigen::Block<Matrix> tnb(Matrix& m) { return m.block(1, 0, m.rows() - 1, 1); }
  static Eigen::Block<Matrix> ntb(Matrix& m) { return m.block(0, 1, 1, m.cols() - 1); }
};

}  // namespace bogus

#endif
