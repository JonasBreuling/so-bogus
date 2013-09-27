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


#ifndef BOGUS_BLOCK_CUGAUSS_SEIDEL_HPP
#define BOGUS_BLOCK_CUGAUSS_SEIDEL_HPP

#include "mat.h"

#include "Core/BlockSolvers/GaussSeidelBase.hpp"

#include <vector>

namespace bogus
{

//! Cuda-based projected Gauss-Seidel iterative solver.
template < typename BlockMatrixType >
class CuGaussSeidel : public GaussSeidelBase< CuGaussSeidel, BlockMatrixType >
{
public:
	typedef GaussSeidelBase< bogus::CuGaussSeidel, BlockMatrixType > Base ;

	typedef typename Base::GlobalProblemTraits GlobalProblemTraits ;
	typedef typename GlobalProblemTraits::Scalar Scalar ;

	//! Default constructor -- you will have to call setMatrix() before using the solve() function
	CuGaussSeidel( ) : Base(),
          m_mat( m_ctx ), m_rhs( m_ctx ), m_res( m_ctx )
        { }

	void setMatrix( const BlockMatrixBase< BlockMatrixType > & matrix ) ;

	template < typename NSLaw, typename RhsT, typename ResT >
	Scalar solve( const NSLaw &law, const RhsT &b, ResT &x, bool tryZeroAsWell = true ) const ;

protected:

	typedef typename Base::Index Index ;

	using Base::m_matrix ;
	using Base::m_maxIters ;
	using Base::m_tol ;
	using Base::m_scaling ;
	using Base::m_maxThreads ;
	using Base::m_evalEvery ;
	using Base::m_skipTol ;
	using Base::m_skipIters ;
	using Base::m_localMatrices ;
	using Base::m_regularization ;

        CudaContext m_ctx ;
        CudaMatrix< Scalar, typename BlockMatrixType::Index, BlockMatrixType::RowsPerBlock >  m_mat ;
        mutable CudaVector< Scalar, typename BlockMatrixType::Index, BlockMatrixType::RowsPerBlock >  m_rhs ;
        mutable CudaVector< Scalar, typename BlockMatrixType::Index, BlockMatrixType::RowsPerBlock >  m_res ;

        // Stroage space
        std::vector< int > m_ris ;
        std::vector< int > m_cis ;
} ;



} //namespace bogus




#endif
