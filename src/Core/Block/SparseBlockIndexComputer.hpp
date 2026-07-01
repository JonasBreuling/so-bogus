/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_SPARSE_BLOCK_INDEX_COMPUTER_HPP
#define BOGUS_SPARSE_BLOCK_INDEX_COMPUTER_HPP

#include "CompoundSparseBlockIndex.hpp"

namespace bogus {

// Index getter

template <typename Derived, bool Major>
struct SparseBlockIndexGetter {
  using MatrixType = SparseBlockMatrixBase<Derived>;
  using ReturnType = typename MatrixType::UncompressedIndexType;

  static ReturnType& get(MatrixType& matrix) { return matrix.m_minorIndex; }

  static const ReturnType& get(const MatrixType& matrix) { return matrix.minorIndex(); }

  static const ReturnType& getOrCompute(const MatrixType& matrix,
                                        typename MatrixType::UncompressedIndexType& tempIndex) {
    return matrix.getOrComputeMinorIndex(tempIndex);
  }
};

template <typename Derived>
struct SparseBlockIndexGetter<Derived, true> {
  using MatrixType = SparseBlockMatrixBase<Derived>;
  using ReturnType = typename MatrixType::MajorIndexType;

  static ReturnType& get(MatrixType& matrix) { return matrix.m_majorIndex; }
  static const ReturnType& get(const MatrixType& matrix) { return matrix.majorIndex(); }

  static const ReturnType& getOrCompute(const MatrixType& matrix, typename MatrixType::UncompressedIndexType&) {
    return matrix.majorIndex();
  }
};

// Index computer

template <typename MatrixType, bool ColWise, bool Transpose,
          bool Symmetric = BlockMatrixTraits<MatrixType>::is_symmetric>
struct SparseBlockIndexComputer {
  using Traits = BlockMatrixTraits<MatrixType>;
  enum { is_major = (ColWise != Transpose) == bool(Traits::is_col_major) };
  using Getter = SparseBlockIndexGetter<MatrixType, is_major>;
  using ReturnType = typename Getter::ReturnType;

  SparseBlockIndexComputer(const SparseBlockMatrixBase<MatrixType>& matrix) : m_matrix(matrix.derived()) {}

  const ReturnType& get() { return Getter::getOrCompute(m_matrix, m_aux); }

 private:
  const MatrixType& m_matrix;
  typename MatrixType::UncompressedIndexType m_aux;
};

template <typename MatrixType, bool ColWise, bool Transpose>
struct SparseBlockIndexComputer<MatrixType, ColWise, Transpose, true> {
  using Traits = BlockMatrixTraits<MatrixType>;
  enum { is_major = (ColWise == bool(Traits::is_col_major)) };
  using ReturnType = CompoundSparseBlockIndex<typename MatrixType::MajorIndexType, typename MatrixType::MinorIndexType,
                                              bool(is_major)>;

  SparseBlockIndexComputer(const SparseBlockMatrixBase<MatrixType>& matrix)
      : m_index(matrix.majorIndex(), matrix.minorIndex()) {
    assert(matrix.minorIndex().valid);
  }

  const ReturnType& get() { return m_index; }

 private:
  ReturnType m_index;
};

}  // namespace bogus

#endif
