/*
 * This file is part of So-bogus, a C++ sparse block matrix library and
 * Second Order Cone solver.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * So-bogus is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * So-bogus is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with So-bogus.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef BOGUS_CUDA_LOADERS_H
#define BOGUS_CUDA_LOADERS_H

#include "mat.h"

#include <Core/Block.impl.hpp>

#define COMBINE
  
template < typename Scalar, typename Index, unsigned D,
  typename IndexT, typename BlocksT >
void load(
  const IndexT& sourceIndex, const BlocksT& sourceBlocks,
  bogus::CudaMatrix< Scalar, Index, D > &dest ) 
{
  Eigen::Matrix< Scalar, D, Eigen::Dynamic > data ;

  const std::size_t nBlocks = sourceBlocks.size() ;
  const Scalar *dataPtr  = NULL ;

#ifdef REORDER

    data.resize( D, nBlocks * D ) ;

#pragma omp parallel for
    for( Index i = 0 ; i < sourceIndex.outerSize() ; ++i )
    {
      const Index start = sourceIndex.outerIndexPtr()[ i ] ;
      const Index size  = sourceIndex.size( i ) ;

      for ( Index j = 0 ; j < size  ; ++j )
      {
        for( unsigned k = 0 ; k < D ; ++k )
        {
          data.col( D*start + j + k*size ) = sourceBlocks[ start+j ].col( k ) ;
        }
      }
    }

    dataPtr = data.data() ;
#else
    dataPtr = sourceBlocks[0].data() ;
#endif

  std::cout << " Sending " << nBlocks*D*D*sizeof( Scalar ) / (1<<20) << " MB to GPU memory..." << std::endl ;
  dest.load( nBlocks, sourceIndex.outerSize(), sourceIndex.innerSize(),
    dataPtr, sourceIndex.outerIndexPtr(), sourceIndex.innerIndexPtr() ) ;
}
  
template < typename Scalar, typename Index, unsigned D >
void load(
  const typename bogus::DualFrictionProblem<D>::WType& source,
  bogus::CudaMatrix< Scalar, Index, D > &dest ) 
{
  load( source.majorIndex(), source.blocks(), dest ) ;
}

template < typename Scalar, typename Index, unsigned D >
void load_transposed(
  const typename bogus::DualFrictionProblem<D>::WType& source,
  bogus::CudaMatrix< Scalar, Index, D > &dest ) 
{
  load( source.transposeIndex(), source.transposeBlocks(), dest ) ;
}


template < typename Scalar, typename Index, unsigned D >
void load_combined(
  const typename bogus::DualFrictionProblem<D>::WType& source,
  std::vector< Index >& rowOffsets, 
  std::vector< Index >& columns, 
  bogus::CudaMatrix< Scalar, Index, D > &dest ) 
{

  const std::size_t nBlocks = source.blocks().size() ;
  const std::size_t ntBlocks = source.transposeBlocks().size() ;
  const std::size_t nnz = nBlocks + ntBlocks ;

  Eigen::Matrix< Scalar, Eigen::Dynamic, D, Eigen::RowMajor > data ( D*nnz, D ) ;

  rowOffsets.resize( source.majorIndex().outerSize() + 1 ) ;
  columns   .resize( nnz );

#pragma omp parallel for
    for( Index i = 0 ; i < source.majorIndex().outerSize() ; ++i )
    {
      const Index start  = source.majorIndex().outerIndexPtr()[ i ] ;
      const Index tstart = source.transposeIndex().outerIndexPtr()[ i ] ;
      
      const Index size  = source.majorIndex().size( i ) ;
      const Index tsize = source.transposeIndex().size( i ) ;
 
      rowOffsets[i] = start + tstart ; ;

#ifdef REORDER
        for( unsigned k = 0 ; k < D ; ++k )
        {
          for ( Index j = 0 ; j < size  ; ++j )
          {
            data.row( D*rowOffsets[i] + j + k*( size+tsize ) ) = source.block( start+j ).col( k ) ;
          }
          for ( Index j = 0 ; j < tsize  ; ++j )
          {
            data.row( D*rowOffsets[i] + size + j + k*( size+tsize ) ) = source.transposeBlocks()[ tstart+j ].col( k ) ;
          }
      }
#else
      data.block( rowOffsets[i]*D, 0, D*size, D ) = 
        Eigen::Matrix< Scalar, Eigen::Dynamic, D, Eigen::RowMajor >::Map(
          source.block(start).data(), D*size, D ) ;

      data.block( (rowOffsets[i]+size )*D, 0, D*tsize, D ) = 
        Eigen::Matrix< Scalar, Eigen::Dynamic, D, Eigen::RowMajor >::Map(
          source.transposeBlocks()[tstart].data(), D*tsize, D ) ;
#endif
      std::copy( source.majorIndex().innerIndexPtr() + start, 
        source.majorIndex().innerIndexPtr() + start + size, 
        columns.begin() + rowOffsets[i] ) ;

      std::copy( source.transposeIndex().innerIndexPtr() + tstart, 
        source.transposeIndex().innerIndexPtr() + tstart + tsize, 
        columns.begin() + rowOffsets[i] + size ) ;
      
    }

  rowOffsets.back() = nnz ;

  std::cout << " Sending " << nnz*D*D*sizeof( Scalar ) / (1<<20) << " MB to GPU memory..." << std::endl ;
  dest.load( nnz, source.rowsOfBlocks(), source.colsOfBlocks(),
    data.data(), &rowOffsets[0], &columns[0] ) ;
}

#endif

