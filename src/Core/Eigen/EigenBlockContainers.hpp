/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef BOGUS_EIGEN_BLOCK_CONTAINERS_HPP
#define BOGUS_EIGEN_BLOCK_CONTAINERS_HPP

#include <Eigen/Core>
#include "../Utils/CppTools.hpp"

namespace bogus {

// We do not have to use the specialized allocator if the size is dynamic or not a multiple of 16 bytes
template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct ResizableSequenceContainer<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>> {
  using BlockType = Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>;

  static constexpr bool UseAlignedAllocator = !(Rows == Eigen::Dynamic || Cols == Eigen::Dynamic ||
                                                ((static_cast<std::size_t>(Rows * Cols * sizeof(Scalar)) & 0xf) != 0));

  using Type = std::conditional_t<UseAlignedAllocator, std::vector<BlockType, Eigen::aligned_allocator<BlockType>>,
                                  std::vector<BlockType>>;
};

}  // namespace bogus

#endif
