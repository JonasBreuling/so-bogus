/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_SPARSEBLOCKMATRIX_HPP
#define BOGUS_SPARSEBLOCKMATRIX_HPP

#include "SparseBlockMatrixBase.hpp"

namespace bogus {

//! Specialization of BlockMatrixTraits for SparseBlockMatrix
template <typename BlockT, int Flags>
struct BlockMatrixTraits<SparseBlockMatrix<BlockT, Flags> >
    : public BlockMatrixTraits<BlockObjectBase<SparseBlockMatrix<BlockT, Flags> > > {
  using BaseTraits = BlockMatrixTraits<BlockObjectBase<SparseBlockMatrix<BlockT, Flags> > >;
  using Index = typename BaseTraits::Index;
  using Scalar = typename BlockTraits<BlockT>::Scalar;

  enum {
    is_symmetric = !!(Flags & flags::SYMMETRIC),
    RowsPerBlock = BlockTraits<BlockT>::RowsAtCompileTime,
    ColsPerBlock = BlockTraits<BlockT>::ColsAtCompileTime
  };

  using BlockPtr = BOGUS_DEFAULT_BLOCK_PTR_TYPE;
  using BlockType = BlockT;
  using BlockRef = BlockT&;
  using ConstBlockRef = const BlockT&;

  template <typename OtherBlockType>
  struct ResizableBlockContainer {
    using Type = typename ResizableSequenceContainer<OtherBlockType>::Type;
  };

  using BlocksArrayType = typename ResizableBlockContainer<BlockType>::Type;

  enum { is_compressed = !!(~Flags & flags::UNCOMPRESSED), is_col_major = !!(Flags & flags::COL_MAJOR), flags = Flags };

  using MajorIndexType = SparseBlockIndex<is_compressed, Index, BlockPtr>;
};

//! Sparse Block Matrix
/*!
  \tparam BlockT the type of the blocks of the matrix. Can be scalar, Eigen dense of sparse matrices,
  or basically anything provided a few functions are specialized
  \tparam Flags a combination of the values defined in \ref bogus::flags
  */
template <typename BlockT, int Flags>
class SparseBlockMatrix : public SparseBlockMatrixBase<SparseBlockMatrix<BlockT, Flags> > {
 public:
  using Base = SparseBlockMatrixBase<SparseBlockMatrix>;

  SparseBlockMatrix() : Base() {}

  template <typename Index>
  SparseBlockMatrix(Index rowsOfBlocks, Index colsOfBlocks) : Base() {
    static_assert(Base::has_fixed_size_blocks, "BLOCKS_MUST_HAVE_FIXED_DIMENSIONS");
    Base::setRows(rowsOfBlocks);
    Base::setCols(colsOfBlocks);
  }

  template <typename RhsT>
  SparseBlockMatrix(const BlockObjectBase<RhsT>& rhs) : Base() {
    Base::operator=(rhs.derived());
  }

  template <typename RhsT>
  SparseBlockMatrix& operator=(const BlockObjectBase<RhsT>& rhs) {
    return (Base::operator=(rhs.derived())).derived();
  }

  // Allocations and stuff

  using Base::m_blocks;
  using Base::m_transposeBlocks;
  using Traits = typename Base::Traits;
  using Index = typename Base::Index;
  using MajorIndexType = typename Base::MajorIndexType;
  using BlockPtr = typename Base::BlockPtr;

  template <bool EnforceThreadSafety>
  typename Base::BlockRef allocateBlock(typename Base::BlockPtr& ptr, Index, Index) {
#ifndef BOGUS_DONT_PARALLELIZE
    Lock::Guard<EnforceThreadSafety> guard(Base::m_lock);
#endif

    ptr = m_blocks.size();
    m_blocks.resize(ptr + 1);

    return m_blocks[ptr];
  }

  //! Resizes \c m_blocks
  void prealloc(std::size_t nBlocks) {
    Base::clear();
    m_blocks.resize(nBlocks);
    Base::m_minorIndex.valid = false;
  }

  //! Reserve enough memory to accomodate \p nBlocks
  void reserve(std::size_t nBlocks, std::size_t = 0) {
    m_blocks.reserve(nBlocks);
    Base::majorIndex().reserve(nBlocks);
  }

  // Inherited from Base

  template <typename BlocksType>
  void resetFor(const BlocksType& blocks) {
    Base::clear();
    reserve(blocks.size());
  }

  template <bool Transpose, typename BlocksType>
  void copyBlockShapes(const BlocksType& blocks) {
    m_blocks.resize(blocks.size());
  }

  template <bool Transpose, typename BlocksType>
  void copyTransposeBlockShapes(const BlocksType& blocks) {
    m_transposeBlocks.resize(blocks.size());
  }

  template <typename IndexType>
  void createBlockShapes(const BlockPtr nBlocks, const IndexType&, typename Traits::BlocksArrayType& blocks) {
    blocks.resize(nBlocks);
  }
};

// Specialization for block matrix of MappedSparseBlockMatrix
template <typename BlockT, int Flags>
struct BlockTraits<SparseBlockMatrix<BlockT, Flags> > {
  using BlockType = SparseBlockMatrix<BlockT, Flags>;

  using Scalar = typename BlockType::Scalar;
  using TransposeStorageType = SparseBlockMatrix<typename BlockType::TransposeBlockType, Flags ^ flags::COL_MAJOR>;

  enum {
    RowsAtCompileTime = internal::DYNAMIC,
    ColsAtCompileTime = internal::DYNAMIC,
    uses_plain_array_storage = 0,
    is_row_major = !BlockMatrixTraits<BlockType>::is_col_major,
    is_self_transpose = BlockMatrixTraits<BlockType>::is_symmetric
  };
};

}  // namespace bogus

#endif  // SPARSEBLOCKMATRIX_HH
