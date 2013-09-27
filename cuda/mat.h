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

#ifndef BOGUS_CUDA_MAT_H
#define BOGUS_CUDA_MAT_H

#include <cuda_runtime.h>
#include <vector>

//#define USE_CUSPARSE
//#define REORDER

namespace bogus {

struct CusparseData ;

struct CudaContext
{  
  CudaContext() ;
  ~CudaContext() ;

  CusparseData *cusparse ;
  std::vector< cudaStream_t > streams  ;

  private:
  CudaContext( const CudaContext& ) ;
  CudaContext& operator=( const CudaContext& ) ;
} ;

struct CudaStream
{
  CudaStream() ;
  ~CudaStream() ;
} ;

template < typename Scalar, typename Index, unsigned Dimension >
struct CudaVector
{
  CudaContext &ctx ;

  Scalar *d_data ;
  Index   rows ;

  CudaVector( CudaContext &context, Index rowsOfBlocks = 0 ) 
    : ctx( context )
  { resize( rowsOfBlocks ) ; }
  ~CudaVector() ;
  
  void resize( Index rowsOfBlocks ) ;
  
  bool valid() const { return 0 != rows ; } 

  void set( const Scalar *rhs, unsigned strteamId = 0 ) ;
  void set( Index row, const Scalar *rhs, unsigned strteamId = 0 ) ;

  void reset( unsigned strteamId = 0 ) ;
  void resetSeg( Index row, unsigned strteamId = 0 ) ;

  void get( Scalar *res, unsigned strteamId = 0 ) const ;
  void get( Index row, Scalar *res, unsigned strteamId = 0 ) const ;

private:
  CudaVector( const CudaVector& ) ;
  CudaVector &operator=( const CudaVector& ) ;

} ;

template < typename Scalar, typename Index, unsigned Dimension >
struct CudaMatrix
{
  CudaContext &ctx ;

  std::size_t nnz ;
  Index rows ;
  Index cols ;
  Index totBlocks ;
  const Index* rowOffsets;

  Scalar *d_data ;
  Index  *d_columns ;
  Index  *d_rowOffsets ;
  Index  *d_idvRowOffsets ;

  typedef CudaVector< Scalar, Index, Dimension > Vector ;

  CudaMatrix( CudaContext &context ) 
  : ctx( context ), rowOffsets( 0 ) {} 
  ~CudaMatrix() ;

  void load( 
      std::size_t numberOfNonZeros, 
      Index rowsOfBlocks, Index colsOfBlocks,
      const Scalar* dataPtr,           
      const Index*  outerIndexPtr,    
      const Index*  innerIndexPtr    
  ) ;

  bool valid() const { return 0 != rowOffsets ; } 

  void mult( const Vector& rhs, Vector& res, unsigned strteamId = 0 ) const ;
  void rowMult( Index row, const Vector& rhs, Vector& res, unsigned strteamId = 0  ) const ;

private:
  CudaMatrix( const CudaMatrix& ) ;
  CudaMatrix &operator=( const CudaMatrix& ) ;

} ;

} 

#endif

