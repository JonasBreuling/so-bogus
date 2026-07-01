/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Preconditioners.hpp"
#include "../Utils/NumTraits.hpp"

#ifndef BOGUS_PRECONDITIONERS_IMPL_HPP
#define BOGUS_PRECONDITIONERS_IMPL_HPP

namespace bogus {

template <typename BlockMatrixType>
class DiagonalPreconditioner<BlockObjectBase<BlockMatrixType> > {
 public:
  using Traits = BlockMatrixTraits<BlockMatrixType>;
  using Scalar = typename Traits::Scalar;
  using BlockType = typename Traits::BlockType;
  using GlobalProblemTraits = ProblemTraits<Scalar>;
  using Vector = typename GlobalProblemTraits::DynVector;

  void setMatrix(const BlockMatrixBase<BlockMatrixType> &matrix) {
    using Index = typename Traits::Index;

    m_diagonal.resize(matrix.rows());
    for (Index i = 0, cur_row = 0; i < matrix.rowsOfBlocks(); ++i) {
      const BlockType &b = matrix.diagonal(i);
      m_diagonal.segment(cur_row, b.rows()) = b.diagonal();
      for (Index k = 0; k < (Index)b.rows(); ++k, ++cur_row) {
        if (NumTraits<Scalar>::isSquareZero(m_diagonal[cur_row])) {
          m_diagonal[cur_row] = 1;
        }
      }
    }
  }

  template <bool transpose, typename ResT, typename RhsT>
  void apply(const RhsT &rhs, ResT &res) const {
    res = rhs.cwiseQuotient(m_diagonal);
  }

 private:
  Vector m_diagonal;
};

template <typename BlockMatrixType, typename FactorizationType>
class DiagonalFactorizationPreconditioner {
 public:
  using BlockType = typename BlockMatrixTraits<BlockMatrixType>::BlockType;
  using Index = typename BlockMatrixTraits<BlockMatrixType>::Index;

  template <bool transpose, typename ResT, typename RhsT>
  void apply(const RhsT &rhs, ResT &res) const {
    static_assert(!transpose, "TRANSPOSE_MAKES_NO_SENSE_IN_THIS_CONTEXT");

    res.setZero();
    m_fact.template multiply<transpose>(rhs, res);
  }

  void setMatrix(const BlockMatrixBase<BlockMatrixType> &matrix) {
    m_fact.cloneDimensions(matrix);

    for (Index i = 0; i < matrix.rowsOfBlocks(); ++i) {
      m_fact.insertBack(i, i).compute(matrix.diagonal(i));
    }
    m_fact.finalize();
  }

 private:
  SparseBlockMatrix<FactorizationType> m_fact;
};

template <typename BlockMatrixType>
class DiagonalLUPreconditioner<BlockObjectBase<BlockMatrixType> >
    : public DiagonalFactorizationPreconditioner<
          BlockMatrixType, typename MatrixTraits<typename BlockMatrixTraits<BlockMatrixType>::BlockType>::LUType> {};

template <typename BlockMatrixType>
class DiagonalLDLTPreconditioner<BlockObjectBase<BlockMatrixType> >
    : public DiagonalFactorizationPreconditioner<
          BlockMatrixType, typename MatrixTraits<typename BlockMatrixTraits<BlockMatrixType>::BlockType>::LDLTType> {};

}  // namespace bogus

#endif
