/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_POLYNOMIAL_IMPL_HPP
#define BOGUS_POLYNOMIAL_IMPL_HPP

#ifndef BOGUS_WITHOUT_EIGEN
#include <Eigen/Eigenvalues>
#endif

#include "Polynomial.hpp"
#include "NumTraits.hpp"

namespace bogus {

namespace polynomial {

#ifndef BOGUS_WITHOUT_EIGEN

// Companion-matrix factory
//
// The companion matrix is partially constant (sub-diagonal identity block) and
// partially overwritten each call (last column ← -coefficients).  We cache the
// constant skeleton once per thread with thread_local + a lambda initialiser,
// which is both simpler and portable across all C++17 targets.

template <unsigned Dimension, typename Scalar>
struct CompanionMatrix {
  using BaseType = Eigen::Matrix<Scalar, Dimension, Dimension>;
  using ReturnType = BaseType;
  // using ReturnType = typename BaseType::MapType;
  // typedef typename BaseType::MapType ReturnType ;

  static ReturnType get() {
    // Initialised exactly once per thread; the sub-diagonal identity block
    // is constant, so we only need to set it at construction time.
    thread_local ReturnType s_matrix = []() noexcept {
      ReturnType m;
      m.template block<1, Dimension - 1>(0, 0).setZero();
      m.template block<Dimension - 1, Dimension - 1>(1, 0).setIdentity();
      return m;
    }();
    return s_matrix;
  }
};

// Real-root finder via companion-matrix eigenvalues

template <unsigned Dimension, typename Scalar>
unsigned RootsFinder<Dimension, Scalar>::getRealRoots(const Scalar *coeffs, Scalar *realRoots, RealRootsFilter filter) {
  using CM = CompanionMatrix<Dimension, Scalar>;
  typename CM::ReturnType matrix = CM::get();

  // Overwrite the last column with −coefficients.
  matrix.template block<Dimension, 1>(0, Dimension - 1) = -Eigen::Matrix<Scalar, Dimension, 1>::Map(coeffs);

  const auto &ev = Eigen::EigenSolver<typename CM::BaseType>(matrix).eigenvalues();

  unsigned count = 0;
  for (unsigned i = 0; i < Dimension; ++i) {
    if (NumTraits<Scalar>::isZero(std::imag(ev[i]))) {
      const bool discard = (filter == StrictlyPositiveRoots && std::real(ev[i]) <= 0) ||
                           (filter == StrictlyNegativeRoots && std::real(ev[i]) >= 0);
      if (!discard) realRoots[count++] = std::real(ev[i]);
    }
  }
  return count;
}

#else  // BOGUS_WITHOUT_EIGEN

template <unsigned Dimension, typename Scalar>
unsigned RootsFinder<Dimension, Scalar>::getRealRoots(const Scalar* /*coeffs*/, Scalar* /*realRoots*/,
                                                      RealRootsFilter /*filter*/) {
  assert(0 && "bogus::Polynomial::RootsFinder::getRealRoots requires Eigen");
  return 0;
}

#endif  // BOGUS_WITHOUT_EIGEN

// Degenerate (arbitrary leading coefficient) root finder

template <typename Scalar>
struct PossiblyDegenerateRootsFinder<0, Scalar> {
  static unsigned getRealRoots(const Scalar *coeffs, Scalar *realRoots, RealRootsFilter filter = AllRoots) {
    realRoots[0] = Scalar(0);
    return filter == AllRoots && NumTraits<Scalar>::isZero(coeffs[0]);
  }
};

template <unsigned Dimension, typename Scalar>
unsigned PossiblyDegenerateRootsFinder<Dimension, Scalar>::getRealRoots(Scalar *coeffs, Scalar *realRoots,
                                                                        RealRootsFilter filter) {
  if (NumTraits<Scalar>::isZero(coeffs[Dimension])) {
    return PossiblyDegenerateRootsFinder<Dimension - 1, Scalar>::getRealRoots(coeffs, realRoots, filter);
  }

  const Scalar inv = Scalar(1) / coeffs[Dimension];
  for (unsigned k = 0; k < Dimension; ++k) coeffs[k] *= inv;

  return RootsFinder<Dimension, Scalar>::getRealRoots(coeffs, realRoots, filter);
}

}  // namespace polynomial

}  // namespace bogus

#endif
