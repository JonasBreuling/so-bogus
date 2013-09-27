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


#include "mat.h"
#include "helpers.h"

#include <iostream>
#include <assert.h>

#include <cusparse_v2.h>

#define TPB 64
#define NSTREAMS 16

//#define DATATIME


template< unsigned D, typename Scalar, typename Index>
__global__ void rowMul( 
    unsigned md, const Scalar* A, const Index* columns,
    const Scalar *x, Scalar *y )
{
 __shared__ double sh_tmp[TPB*D] ;

 int i = blockIdx.x * blockDim.x + threadIdx.x ;
 Scalar* const tmp = &sh_tmp[ threadIdx.x ] ;
 
 *tmp = 0 ;
 if( i >= md ) return ;

#ifdef REORDER
 x += columns[i/D]*D ;
 A += i ;
 
#pragma unroll
  for( int j = 0 ; j < D ; ++j )
  {
     *tmp += A[ j*md ] * x[j] ;
  }
#else
 x += columns[i/D]*D ;
 A += D*i ;

#pragma unroll
  for( int j = 0 ; j < D ; ++j )
  {
     *tmp += A[ j ] * x[j] ;
  }
#endif

  //printf(" %d %f \n ", i, *tmp ) ;

 __syncthreads() ;
 for( int offset = TPB/2 ; D*offset > threadIdx.x ; offset >>= 1 )
 {
     *tmp += tmp[ D*offset ] ;
     __syncthreads() ;
 }
  //printf(" %d %f \n ", i, *tmp ) ;


 if( threadIdx.x < D ) {
    atomicAdd( y + threadIdx.x, *tmp ) ;
 }
}

template< unsigned D, typename Scalar, typename Index>
__global__ void mul( 
    const int totBlocks,
    const Scalar* A, const Index *offsets, 
    const Index* columns, 
    const Index* blocksData,
    const Scalar *x, Scalar *y )
{
 __shared__ double sh_tmp[TPB*D] ;

 const int b = blockIdx.y * gridDim.x + blockIdx.x ;
 if( b >= totBlocks ) return ;

 const int r = blocksData[ 2*b ];
 const int i = (b - blocksData[ 2*b + 1 ]) * blockDim.x + threadIdx.x ;
;

 //if( threadIdx.x == 0 ) printf( "%d %d %d \n", b, i, r ) ;
// if( r >= n ) return ;
 //int i = x - blockOffsets[r] ;

 Scalar* const tmp = &sh_tmp[ threadIdx.x ] ;
 *tmp = 0 ;

 const int md = D*( offsets[r+1] - offsets[r] ) ;
 if( i >= md ) 
   return ;

 A += offsets[r]*D*D ;
 columns += offsets[r] ;

#ifdef REORDER
 x += columns[i/D]*D ;
 A += i ;
 
#pragma unroll
  for( int j = 0 ; j < D ; ++j )
  {
     *tmp += A[ j*md ] * x[j] ;
  }
#else
 x += columns[i/D]*D ;
 A += D*i ;

#pragma unroll
  for( int j = 0 ; j < D ; ++j )
  {
     *tmp += A[ j ] * x[j] ;
  }
#endif

  //printf(" %d %f \n ", i, *tmp ) ;

 __syncthreads() ;
#pragma unroll
 for( int offset = TPB/2 ; D*offset > threadIdx.x ; offset >>= 1 )
 {
     *tmp += tmp[ D*offset ] ;
     __syncthreads() ;
 }
  //printf(" %d %f \n ", i, *tmp ) ;


 if( threadIdx.x < D ) {
    atomicAdd( y + D*r + threadIdx.x, *tmp ) ;
 }
}


template< unsigned D, typename Scalar, typename Index>
__global__ void rowMulD( 
    unsigned m, const Scalar* A, const Index* columns,
    const Scalar *x, Scalar *y )
{
 __shared__ double sh_tmp[TPB*D] ;

 int i = blockIdx.x * blockDim.x + threadIdx.x ;

 Scalar* const tmp = sh_tmp + D*threadIdx.x ;

#pragma unroll
  for( int k = 0 ; k < D ; ++k )
  {
    tmp[k] = 0 ;
  }
 
 if( i >= m ) return ;

 x += columns[i]*D ;
 A += D*D*i ;

#pragma unroll
  for( int k = 0 ; k < D ; ++k )
  {
#pragma unroll
    for( int j = 0 ; j < D ; ++j )
    {
      tmp[k] += A[ j ] * x[j] ;
    }
    A += D ;
  }
  //printf(" %d %f \n ", i, *tmp ) ;

 __syncthreads() ;
 for( int offset = TPB/2 ; offset > threadIdx.x ; offset >>= 1 )
 {

#pragma unroll
   for( int k = 0 ; k < D ; ++k )
   {
     tmp[k] += tmp[ D*offset + k ] ;
   }
     __syncthreads() ;
 }
  //printf(" %d %f \n ", i, *tmp ) ;

 if( threadIdx.x == 0 ) {
#pragma unroll
   for( int k = 0 ; k < D ; ++k )
   {
    atomicAdd( y + k, sh_tmp[k] ) ;
   }
 }
}

namespace bogus
{

struct CusparseData
{
  cusparseHandle_t handle ;
  cusparseMatDescr_t descr ;

  CusparseData()
  {
    cusparseCreate(&handle) ;
    cusparseCreateMatDescr(&descr) ;
  }
  
  ~CusparseData()
  {
    cusparseDestroy(handle) ;
    cusparseDestroyMatDescr(descr) ;
  }
} ;

// CONTEXT

CudaContext::CudaContext()
 : cusparse( NULL ) 
{
  cudaSetDevice(0);
  CudaCheckError() ;

  streams.resize( NSTREAMS ) ;
  for( unsigned i = 0 ; i < streams.size() ; ++i )
  {
    cudaStreamCreate( &streams[i] ) ;
  }
  
#ifdef USE_CUSPARSE
  cusparse = new CusparseData() ;
#endif
}

CudaContext::~CudaContext()
{
  for( unsigned i = 0 ; i < streams.size() ; ++i )
  {
    cudaStreamDestroy( streams[i] ) ;
  }

  delete cusparse ;
  cudaDeviceReset();
}

// VECTOR

template < typename Scalar, typename Index, unsigned D >
CudaVector< Scalar, Index, D>::
~CudaVector()
{
  if( valid() )
  {
    cudaFree( d_data ) ;
  }
}

template < typename Scalar, typename Index, unsigned D >
void CudaVector< Scalar, Index, D>::
resize( Index rowsOfBlocks )
{
  if( valid() )
  {
    cudaFree( d_data ) ;
  }
  rows = rowsOfBlocks ;
  if( valid() )
  {
    CudaSafeCall( cudaMalloc( &d_data,   rows * D * sizeof(Scalar)) );
  }
}

template < typename Scalar, typename Index, unsigned D >
void CudaVector< Scalar, Index, D>::
set( const Scalar * data, unsigned sId )
{
 cudaMemcpyAsync( d_data, data, rows * D * sizeof(Scalar), cudaMemcpyHostToDevice, ctx.streams[sId] ) ;
}
template < typename Scalar, typename Index, unsigned D >
void CudaVector< Scalar, Index, D>::
set( Index row, const Scalar * data, unsigned sId )
{
 cudaMemcpyAsync( d_data + row*D, data, D * sizeof(Scalar), cudaMemcpyHostToDevice, ctx.streams[sId] ) ;
}

template < typename Scalar, typename Index, unsigned D >
void CudaVector< Scalar, Index, D>::
reset( unsigned sId) 
{
 cudaMemsetAsync( d_data, 0, rows * D * sizeof(Scalar), ctx.streams[sId] ) ;
}

template < typename Scalar, typename Index, unsigned D >
void CudaVector< Scalar, Index, D>::
resetSeg( Index row, unsigned sId ) 
{
 cudaMemsetAsync( d_data + row*D, 0, D * sizeof(Scalar), ctx.streams[sId] ) ;
}

template < typename Scalar, typename Index, unsigned D >
void CudaVector< Scalar, Index, D>::
get( Scalar * data, unsigned sId ) const
{
 cudaMemcpyAsync( data, d_data, rows * D * sizeof(Scalar), cudaMemcpyDeviceToHost, ctx.streams[sId] ) ;
}

template < typename Scalar, typename Index, unsigned D >
void CudaVector< Scalar, Index, D>::
get( Index row, Scalar * data, unsigned sId ) const
{
 cudaMemcpyAsync( data, d_data + row*D, D * sizeof(Scalar), cudaMemcpyDeviceToHost, ctx.streams[sId] ) ;
}

//MATRIX

template < typename Scalar, typename Index, unsigned D >
void CudaMatrix< Scalar, Index, D>::load( 
    std::size_t numberOfNonZeros, 
    Index rowsOfBlocks, Index colsOfBlocks,
    const Scalar* dataPtr,           
    const Index*  outerIndexPtr,    
    const Index*  innerIndexPtr    
    )
{

  assert( !valid() ) ;
 
  nnz = numberOfNonZeros ;
  rows = rowsOfBlocks ;
  cols = colsOfBlocks ;
  rowOffsets = outerIndexPtr ;
  
  assert( valid() ) ;
 CudaSafeCall( cudaMalloc( &d_data, D*D*numberOfNonZeros * sizeof(Scalar)) );

 CudaSafeCall( cudaMalloc( &d_columns, numberOfNonZeros * sizeof(Index)) );
 
 CudaSafeCall( cudaMemcpy( d_data, dataPtr, D*D*numberOfNonZeros * sizeof(Scalar), cudaMemcpyHostToDevice ) );
 CudaSafeCall( cudaMemcpy( d_columns, innerIndexPtr, numberOfNonZeros * sizeof(Index), cudaMemcpyHostToDevice  ) );

 CudaSafeCall( cudaMalloc( &d_rowOffsets, (rows+1) * sizeof(Index)) );
 CudaSafeCall( cudaMemcpy( d_rowOffsets, rowOffsets, (rows+1) * sizeof(Index), cudaMemcpyHostToDevice ) );

 std::vector< Index > blocksData ;

 totBlocks = 0 ;
 for( int i = 0 ; i < rows ; ++i )
 {
   const int nBlocks = 
     ( rowOffsets[i+1] - rowOffsets[i] ) / TPB + 1 ;

   for( int b = 0 ; b < nBlocks ; ++b )
   {
     blocksData.push_back( i ) ;
     blocksData.push_back( totBlocks ) ;
   }
   totBlocks += nBlocks ;
 }

 CudaSafeCall( cudaMalloc( &d_idvRowOffsets, nnz * sizeof(Index)) );
 CudaSafeCall( cudaMemcpy( d_idvRowOffsets, blocksData.data(), 2 * totBlocks * sizeof(Index), cudaMemcpyHostToDevice ) );

} 

template < typename Scalar, typename Index, unsigned D >
CudaMatrix< Scalar, Index, D>::~CudaMatrix( )
{
  if( valid() )
  {
    cudaFree( d_data ) ;
    cudaFree( d_columns ) ;
    if( d_rowOffsets) 
      cudaFree( d_rowOffsets ) ;
    if( d_idvRowOffsets ) 
      cudaFree( d_idvRowOffsets ) ;

    rowOffsets = 0 ;
  }
} 

  
template < typename Scalar, typename Index, unsigned D >
void CudaMatrix< Scalar, Index, D>::
rowMult( Index row, const Vector& rhs, Vector& res, unsigned sId  ) const
{
  const int m = rowOffsets[row+1] - rowOffsets[row] ;

  assert( rows == res.rows ) ;
  assert( cols == rhs.rows ) ;

#ifdef DATATIME
  const int nThreads = TPB ;
  const int nBlocks  = m/nThreads + 1 ;
  
  rowMulD< D >
  <<< nBlocks, nThreads, 0, ctx.streams[sId] >>>( 
    m, d_data + rowOffsets[row]*D*D, 
    d_columns + rowOffsets[row],
    rhs.d_data, res.d_data + D*row );
#else 

  const int nThreads = TPB*D ;
  const int nBlocks  = (D*m)/nThreads + 1 ;
  rowMul< D >
  <<< nBlocks, nThreads, 0, ctx.streams[sId] >>>( 
    m*D, d_data + rowOffsets[row]*D*D, 
    d_columns + rowOffsets[row],
    rhs.d_data, res.d_data + D*row );
#endif
  
}

template < typename Scalar, typename Index, unsigned D >
void CudaMatrix< Scalar, Index, D>::
mult( const Vector& rhs, Vector& res, unsigned sId ) const
{
  if( ! valid() ) return ;

#ifdef USE_CUSPARSE

  double alpha = 1, beta = 1 ;

  cusparseDbsrmv( ctx.cusparse->handle,
    CUSPARSE_DIRECTION_ROW, 
    CUSPARSE_OPERATION_NON_TRANSPOSE,
    rows, cols, nnz, 
    &alpha, ctx.cusparse->descr, 
    d_data, d_rowOffsets, d_columns, D, 
    rhs.d_data, &beta, res.d_data);

#else
  dim3 dimGrid ( 1<<10, totBlocks/(1<<10) + 1, 1 ) ;
  dim3 dimBlock ( TPB*D, 1, 1 ) ;

  mul< D ><<< dimGrid, dimBlock, 0, ctx.streams[sId] >>>( 
    totBlocks,
    d_data, 
    d_rowOffsets,
    d_columns, d_idvRowOffsets,
    rhs.d_data, res.d_data );

#endif
}


template struct CudaMatrix< double, int, 2u > ;
template struct CudaMatrix< double, int, 3u > ;
template struct CudaVector< double, int, 2u > ;
template struct CudaVector< double, int, 3u > ;

}

