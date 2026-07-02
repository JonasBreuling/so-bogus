/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_EIGEN_LINEAR_SOLVERS_HPP
#define BOGUS_EIGEN_LINEAR_SOLVERS_HPP

#include "../Utils/LinearSolverBase.hpp"

#include <Eigen/LU>
#include <Eigen/Cholesky>

namespace bogus {

template <typename Derived>
struct LinearSolverTraits<LU<Eigen::MatrixBase<Derived> > > {
  using MatrixType = typename Derived::PlainObject;
  using FactType = Eigen::FullPivLU<MatrixType>;

  template <typename RhsT>
  struct Result {
    using Type = const Eigen::Solve<FactType, RhsT>;
  };
  template <typename RhsT>
  struct Result<Eigen::MatrixBase<RhsT> > {
    using Type = typename Result<RhsT>::Type;
  };
};

template <typename Derived>
struct LU<Eigen::MatrixBase<Derived> > : public LinearSolverBase<LU<Eigen::MatrixBase<Derived> > > {
  using MatrixType = Eigen::MatrixBase<Derived>;
  using Traits = LinearSolverTraits<LU<MatrixType> >;

  LU() {}
  template <typename OtherDerived>
  explicit LU(const Eigen::MatrixBase<OtherDerived>& mat) : m_fact(mat) {}

  template <typename OtherDerived>
  LU<Eigen::MatrixBase<Derived> >& compute(const Eigen::MatrixBase<OtherDerived>& mat) {
    m_fact.compute(mat);
    return *this;
  }

  template <typename RhsT, typename ResT>
  void solve(const Eigen::MatrixBase<RhsT>& rhs, ResT& res) const {
    res = m_fact.solve(rhs);
  }

  template <typename RhsT>
  typename Traits::template Result<Eigen::MatrixBase<RhsT> >::Type solve(const Eigen::MatrixBase<RhsT>& rhs) const {
    return m_fact.solve(rhs);
  }

 private:
  typename Traits::FactType m_fact;
};

template <typename Scalar, int Rows, int Cols = Rows, int Options = 0>
struct DenseLU : public LU<Eigen::MatrixBase<Eigen::Matrix<Scalar, Rows, Cols, Options> > > {
  DenseLU() {}
  template <typename OtherDerived>
  explicit DenseLU(const Eigen::MatrixBase<OtherDerived>& mat)
      : LU<Eigen::MatrixBase<Eigen::Matrix<Scalar, Rows, Cols, Options> > >(mat) {}
};

template <typename Derived>
struct LinearSolverTraits<LDLT<Eigen::MatrixBase<Derived> > > {
  using MatrixType = typename Derived::PlainObject;
  using FactType = Eigen::LDLT<MatrixType>;

  template <typename RhsT>
  struct Result {
    using Type = const Eigen::Solve<FactType, RhsT>;
  };
  template <typename RhsT>
  struct Result<Eigen::MatrixBase<RhsT> > {
    using Type = typename Result<RhsT>::Type;
  };
};

template <typename Derived>
struct LDLT<Eigen::MatrixBase<Derived> > : public LinearSolverBase<LDLT<Eigen::MatrixBase<Derived> > > {
  using MatrixType = Eigen::MatrixBase<Derived>;
  using Traits = LinearSolverTraits<LDLT<MatrixType> >;

  LDLT() {}
  template <typename OtherDerived>
  explicit LDLT(const Eigen::MatrixBase<OtherDerived>& mat) : m_fact(mat) {}

  template <typename OtherDerived>
  LDLT<Eigen::MatrixBase<Derived> >& compute(const Eigen::MatrixBase<OtherDerived>& mat) {
    m_fact.compute(mat);
    return *this;
  }

  template <typename RhsT, typename ResT>
  void solve(const Eigen::MatrixBase<RhsT>& rhs, ResT& res) const {
    res = m_fact.solve(rhs);
  }

  template <typename RhsT>
  typename Traits::template Result<Eigen::MatrixBase<RhsT> >::Type solve(const Eigen::MatrixBase<RhsT>& rhs) const {
    return m_fact.solve(rhs);
  }

 private:
  typename Traits::FactType m_fact;
};

template <typename Scalar, int Rows, int Options = 0>
struct DenseLDLT : public LDLT<Eigen::MatrixBase<Eigen::Matrix<Scalar, Rows, Rows, Options> > > {
  DenseLDLT() {}
  template <typename OtherDerived>
  explicit DenseLDLT(const Eigen::MatrixBase<OtherDerived>& mat)
      : LDLT<Eigen::MatrixBase<Eigen::Matrix<Scalar, Rows, Rows, Options> > >(mat) {}
};

template <typename Derived, typename RhsT>
typename LinearSolverTraits<Derived>::template Result<Eigen::MatrixBase<RhsT> >::Type operator*(
    const LinearSolverBase<Derived>& solver, const Eigen::MatrixBase<RhsT>& rhs) {
  return solver.solve(rhs);
}

}  // namespace bogus

#endif
