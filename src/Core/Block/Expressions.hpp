/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_BLOCK_EXPRESSIONS_HPP
#define BOGUS_BLOCK_EXPRESSIONS_HPP

#include "BlockObjectBase.hpp"

namespace bogus {

//! Base class for Transpose views of a BlockObjectBase
template <typename MatrixT>
struct Transpose : public BlockObjectBase<Transpose<MatrixT> > {
  using Base = BlockObjectBase<Transpose<MatrixT> >;
  using Traits = BlockMatrixTraits<Transpose<MatrixT> >;
  using PlainObjectType = typename Traits::PlainObjectType;
  using Index = typename Traits::Index;
  using Scalar = typename Base::Scalar;

  const PlainObjectType& matrix;

  Transpose(const PlainObjectType& m) : matrix(m.derived()) {}

  typename Base::ConstTransposeReturnType transpose() const { return matrix; }

  template <bool DoTranspose, typename RhsT, typename ResT>
  void multiply(const RhsT& rhs, ResT& res, Scalar alpha = 1, Scalar beta = 0) const {
    matrix.template multiply<!DoTranspose>(rhs, res, alpha, beta);
  }

  Index rows() const { return matrix.cols(); }
  Index cols() const { return matrix.rows(); }

  Index rowsOfBlocks() const { return matrix.colsOfBlocks(); }
  Index colsOfBlocks() const { return matrix.rowsOfBlocks(); }
  Index blockRows(Index row) const { return matrix.blockCols(row); }
  Index blockCols(Index col) const { return matrix.blockRows(col); }
  const Index* rowOffsets() const { return matrix.colOffsets(); }
  const Index* colOffsets() const { return matrix.rowOffsets(); }
};

template <typename MatrixT>
struct BlockMatrixTraits<Transpose<MatrixT> > {
  using OrigTraits = BlockMatrixTraits<MatrixT>;
  using Index = typename OrigTraits::Index;
  using Scalar = typename OrigTraits::Scalar;

  enum { is_transposed = 1, is_temporary = 1, is_symmetric = OrigTraits::is_symmetric };
  enum { RowsPerBlock = OrigTraits::ColsPerBlock, ColsPerBlock = OrigTraits::RowsPerBlock };

  using PlainObjectType = typename OrigTraits::PlainObjectType;
  using ConstTransposeReturnType = const PlainObjectType&;
  using TransposeObjectType = PlainObjectType;
};

template <typename ObjectT, bool IsTemporary>
struct BlockStorage {
  using ConstValue = const ObjectT&;
};
template <typename ObjectT>
struct BlockStorage<ObjectT, true> {
  using ConstValue = const ObjectT;
};

template <typename ObjectT>
struct BlockOperand {
  using ObjectType = ObjectT;
  using PlainObjectType = typename ObjectT::PlainObjectType;

  using Traits = BlockMatrixTraits<ObjectT>;
  enum { do_transpose = Traits::is_transposed };
  using Scalar = typename Traits::Scalar;

  typename BlockStorage<ObjectT, Traits::is_temporary>::ConstValue object;
  Scalar scaling;

  BlockOperand(const ObjectT& o, Scalar s = 1) : object(o), scaling(s) {}
};

template <template <typename, typename> class BlockOp, typename LhsMatrixT, typename RhsMatrixT>
struct BinaryBlockOp : public BlockObjectBase<BlockOp<LhsMatrixT, RhsMatrixT> > {
  using Base = BlockObjectBase<BlockOp<LhsMatrixT, RhsMatrixT> >;
  using PlainObjectType = typename Base::PlainObjectType;

  using Lhs = BlockOperand<LhsMatrixT>;
  using Rhs = BlockOperand<RhsMatrixT>;

  using PlainLhsMatrixType = typename Lhs::PlainObjectType;
  using PlainRhsMatrixType = typename Rhs::PlainObjectType;

  const Lhs lhs;
  const Rhs rhs;
  enum { transposeLhs = Lhs::do_transpose, transposeRhs = Rhs::do_transpose };

  BinaryBlockOp(const LhsMatrixT& l, const RhsMatrixT& r, typename Lhs::Scalar lscaling = 1,
                typename Lhs::Scalar rscaling = 1)
      : lhs(l, lscaling), rhs(r, rscaling) {}
};

template <typename LhsMatrixT, typename RhsMatrixT>
struct Product : public BinaryBlockOp<Product, LhsMatrixT, RhsMatrixT> {
  using Base = BinaryBlockOp<bogus::Product, LhsMatrixT, RhsMatrixT>;
  using Scalar = typename Base::Scalar;
  using Index = typename Base::Index;

  Product(const LhsMatrixT& l, const RhsMatrixT& r, typename Base::Lhs::Scalar lscaling = 1,
          typename Base::Lhs::Scalar rscaling = 1)
      : Base(l, r, lscaling, rscaling) {}

  typename Base::ConstTransposeReturnType transpose() const {
    return typename Base::ConstTransposeReturnType(Base::rhs.object.transpose(), Base::lhs.object.transpose(),
                                                   Base::rhs.scaling, Base::lhs.scaling);
  }

  template <bool DoTranspose, typename RhsT, typename ResT>
  void multiply(const RhsT& rhs, ResT& res, Scalar alpha = 1, Scalar beta = 0) const;

  Index rows() const { return Base::lhs.object.rows(); }
  Index cols() const { return Base::rhs.object.cols(); }

  Index rowsOfBlocks() const { return Base::lhs.object.rowsOfBlocks(); }
  Index colsOfBlocks() const { return Base::rhs.object.colsOfBlocks(); }
  Index blockRows(Index row) const { return Base::lhs.object.blockRows(row); }
  Index blockCols(Index col) const { return Base::rhs.object.blockCols(col); }
  const Index* rowOffsets() const { return Base::lhs.object.rowOffsets(); }
  const Index* colOffsets() const { return Base::rhs.object.colOffsets(); }
};

template <typename LhsMatrixT, typename RhsMatrixT>
struct BlockMatrixTraits<Product<LhsMatrixT, RhsMatrixT> > {
  using LhsTraits = BlockMatrixTraits<LhsMatrixT>;
  using RhsTraits = BlockMatrixTraits<RhsMatrixT>;

  using Index = typename LhsTraits::Index;
  using Scalar = typename LhsTraits::Scalar;

  enum { is_transposed = 0, is_temporary = 1, is_symmetric = 0 };

  using ProductType = Product<LhsMatrixT, RhsMatrixT>;

  using LhsBlockType = typename BlockMatrixTraits<typename LhsTraits::PlainObjectType>::BlockType;
  using RhsBlockType = typename BlockMatrixTraits<typename RhsTraits::PlainObjectType>::BlockType;

  using ResBlockType = typename BlockBlockProductTraits<LhsBlockType, RhsBlockType, LhsTraits::is_transposed,
                                                        RhsTraits::is_transposed>::ReturnType;

  using PlainObjectType = typename LhsTraits::PlainObjectType ::template MutableImpl<ResBlockType, false>::Type;

  enum { RowsPerBlock = LhsTraits::RowsPerBlock, ColsPerBlock = RhsTraits::ColsPerBlock };

  using ConstTransposeReturnType = Product<typename BlockOperand<RhsMatrixT>::ObjectType::TransposeObjectType,
                                           typename BlockOperand<LhsMatrixT>::ObjectType::TransposeObjectType>;
  using TransposeObjectType = ConstTransposeReturnType;
};

template <typename LhsMatrixT, typename RhsMatrixT>
struct Addition : public BinaryBlockOp<Addition, LhsMatrixT, RhsMatrixT> {
  using Base = BinaryBlockOp<bogus::Addition, LhsMatrixT, RhsMatrixT>;
  using Scalar = typename Base::Scalar;
  using Index = typename Base::Index;

  Addition(const LhsMatrixT& l, const RhsMatrixT& r, typename Base::Lhs::Scalar lscaling = 1,
           typename Base::Lhs::Scalar rscaling = 1)
      : Base(l, r, lscaling, rscaling) {}

  typename Base::ConstTransposeReturnType transpose() const {
    return typename Base::ConstTransposeReturnType(Base::lhs.object.transpose(), Base::rhs.object.transpose(),
                                                   Base::lhs.scaling, Base::rhs.scaling);
  }

  template <bool DoTranspose, typename RhsT, typename ResT>
  void multiply(const RhsT& rhs, ResT& res, Scalar alpha = 1, Scalar beta = 0) const {
    Base::lhs.object.template multiply<DoTranspose>(rhs, res, alpha * Base::lhs.scaling, beta);
    Base::rhs.object.template multiply<DoTranspose>(rhs, res, alpha * Base::rhs.scaling, 1);
  }

  Index rows() const { return Base::lhs.object.rows(); }
  Index cols() const { return Base::lhs.object.cols(); }

  Index rowsOfBlocks() const { return Base::lhs.object.rowsOfBlocks(); }
  Index colsOfBlocks() const { return Base::lhs.object.colsOfBlocks(); }
  Index blockRows(Index row) const { return Base::lhs.object.blockRows(row); }
  Index blockCols(Index col) const { return Base::lhs.object.blockCols(col); }
  const Index* rowOffsets() const { return Base::lhs.object.rowOffsets(); }
  const Index* colOffsets() const { return Base::lhs.object.colOffsets(); }
};

template <typename LhsMatrixT, typename RhsMatrixT>
struct BlockMatrixTraits<Addition<LhsMatrixT, RhsMatrixT> > {
  using OrigTraits = BlockMatrixTraits<LhsMatrixT>;
  using Index = typename OrigTraits::Index;
  using Scalar = typename OrigTraits::Scalar;

  using ResBlockType = typename BlockMatrixTraits<typename OrigTraits::PlainObjectType>::BlockType;

  using PlainObjectType = typename OrigTraits::PlainObjectType ::template MutableImpl<ResBlockType, false>::Type;

  enum {
    is_transposed = 0,
    is_temporary = 1,
    is_symmetric = (BlockMatrixTraits<LhsMatrixT>::is_symmetric && BlockMatrixTraits<RhsMatrixT>::is_symmetric)
  };
  enum { RowsPerBlock = OrigTraits::RowsPerBlock, ColsPerBlock = OrigTraits::ColsPerBlock };

  using ConstTransposeReturnType = Addition<typename BlockOperand<LhsMatrixT>::ObjectType::TransposeObjectType,
                                            typename BlockOperand<RhsMatrixT>::ObjectType::TransposeObjectType>;
  using TransposeObjectType = ConstTransposeReturnType;
};

template <typename MatrixT>
struct Scaling : public BlockObjectBase<Scaling<MatrixT> > {
  using Operand = BlockOperand<MatrixT>;
  using PlainOperandMatrixType = typename Operand::PlainObjectType;
  Operand operand;

  enum { transposeOperand = Operand::do_transpose };

  using Base = BlockObjectBase<Scaling>;

  using Scalar = typename Base::Scalar;
  using Index = typename Base::Index;

  Scaling(const MatrixT& object, const typename MatrixT::Scalar scaling) : operand(object, scaling) {}

  typename Base::ConstTransposeReturnType transpose() const {
    return typename Base::ConstTransposeReturnType(operand.object.transpose(), operand.scaling);
  }

  template <bool DoTranspose, typename RhsT, typename ResT>
  void multiply(const RhsT& rhs, ResT& res, Scalar alpha = 1, Scalar beta = 0) const {
    operand.object.template multiply<DoTranspose>(rhs, res, alpha * operand.scaling, beta);
  }

  typename Base::Index rows() const { return operand.object.rows(); }
  typename Base::Index cols() const { return operand.object.cols(); }

  Index rowsOfBlocks() const { return operand.object.rowsOfBlocks(); }
  Index colsOfBlocks() const { return operand.object.colsOfBlocks(); }
  Index blockRows(Index row) const { return operand.object.blockRows(row); }
  Index blockCols(Index col) const { return operand.object.blockCols(col); }
  const Index* rowOffsets() const { return operand.object.rowOffsets(); }
  const Index* colOffsets() const { return operand.object.colOffsets(); }
};

template <typename MatrixT>
struct BlockMatrixTraits<Scaling<MatrixT> > {
  using OrigTraits = BlockMatrixTraits<MatrixT>;
  using Index = typename OrigTraits::Index;
  using Scalar = typename OrigTraits::Scalar;

  enum { is_symmetric = OrigTraits::is_symmetric, is_transposed = 0, is_temporary = 1 };
  enum { RowsPerBlock = OrigTraits::RowsPerBlock, ColsPerBlock = OrigTraits::ColsPerBlock };

  using PlainObjectType = typename OrigTraits::PlainObjectType;

  using ConstTransposeReturnType = Scaling<typename BlockOperand<MatrixT>::ObjectType::TransposeObjectType>;
  using TransposeObjectType = ConstTransposeReturnType;
};

template <typename ObjectT>
struct BlockOperand<Scaling<ObjectT> > : public BlockOperand<ObjectT> {
  using Base = BlockOperand<ObjectT>;

  BlockOperand(const Scaling<ObjectT>& o, typename Base::Scalar s = 1)
      : Base(o.operand.object, s * o.operand.scaling) {}
  BlockOperand(const typename Base::ObjectType& o, typename Base::Scalar s = 1) : Base(o, s) {}
};

}  // namespace bogus

#endif  // EXPRESSIONS_HPP
