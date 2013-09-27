
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

#ifndef BOGUS_BLOCK_CUGAUSS_SEIDEL_IMPL_HPP
#define BOGUS_BLOCK_CUGAUSS_SEIDEL_IMPL_HPP

#include "loaders.hpp"

#include "Core/BlockSolvers/GaussSeidelBase.impl.hpp"

#include <omp.h>

namespace bogus {


template < typename BlockMatrixType >
void CuGaussSeidel< BlockMatrixType >::setMatrix( const BlockMatrixBase< BlockMatrixType > & M )
{
	m_matrix = &M ;

        m_rhs.resize( M.rowsOfBlocks() ) ;
        m_res.resize( M.rowsOfBlocks() ) ;
        m_res.reset() ;
        load_combined( M, m_ris, m_cis, m_mat ) ;

	Base::updateLocalMatrices() ;
}

template < typename BlockMatrixType >
template < typename NSLaw, typename RhsT, typename ResT >
typename CuGaussSeidel< BlockMatrixType >::Scalar CuGaussSeidel< BlockMatrixType >::solve( const NSLaw &law,
  const RhsT &b, ResT &x, bool ) const
{
	assert( m_matrix ) ;
	typedef typename NSLaw::Traits LocalProblemTraits ;
	const Segmenter< NSLaw::dimension, const RhsT, typename BlockMatrixType::Index >
			bSegmenter( b, m_matrix->rowOffsets() ) ;
	Segmenter< NSLaw::dimension, ResT, typename BlockMatrixType::Index >
			xSegmenter( x, m_matrix->rowOffsets() ) ;

        m_rhs.set( x.data() ) ;
	typename GlobalProblemTraits::DynVector y ( x.rows() ), x_best( x ) ;
        
        m_mat.mult( m_rhs, m_res ) ;
        m_res.get( y.data() ) ;
        m_res.reset(15) ;
        y += b ;

	Scalar err_best = Base::eval( law, y, x      ) ;

	this->m_callback.trigger( 0, err_best ) ;

	const Index n = m_matrix->rowsOfBlocks() ;
	std::vector< unsigned char > skip( n, 0 ) ;

	unsigned GSIter ;
	for( GSIter = 1 ; GSIter <= m_maxIters ; ++GSIter )
	{
          typename LocalProblemTraits::Vector lb, lr, lx, ldx ;

          #pragma omp parallel for private( lb, lr, lx, ldx )
          for( std::ptrdiff_t i = 0 ; i < n ; ++ i )
          {

            if( skip[i] ) {
              --skip[i] ;
              continue ;
            }

            const unsigned sid = omp_get_thread_num() ;
            
            m_mat.rowMult( i, m_rhs, m_res, sid ) ;

            lx = xSegmenter[ i ] ;
            lb = bSegmenter[ i ] - m_localMatrices[i] * lx ;

            ldx = -lx ;
            
            m_res.get( i, lr.data(), sid ) ;
            m_res.resetSeg(i, sid) ;

            const bool ok = law.solveLocal( i, m_localMatrices[i], lb+lr, lx, m_scaling[ i ] ) ;
            ldx += lx ;

            if( !ok ) { ldx *= .5 ; }
            xSegmenter[ i ] += ldx ;
            m_rhs.set( i, xSegmenter[ i ].data(), sid+8 ) ;

            if( m_scaling[ i ] * m_scaling[ i ] * ldx.squaredNorm() < m_skipTol ||
                m_scaling[ i ] * m_scaling[ i ] *  lx.squaredNorm() < m_tol )
            {
              skip[i] = m_skipIters ;
            }
          }

          if( 0 == ( GSIter % m_evalEvery ) )
          {
            m_mat.mult( m_rhs, m_res ) ;
            m_res.get( y.data() ) ;
            m_res.reset(15) ;
            y += b ;

            const double err = Base::eval( law, y, x ) ;

            this->m_callback.trigger( GSIter, err ) ;

            if( err < m_tol )
            {
              err_best = err ;
              break ;
            }

            if( err < err_best )
            {
              x_best = x ;
              err_best = err ;
            }
          }

        }

	if( GSIter > m_maxIters ) x = x_best ;

	return err_best ;

}


} //namespace bogus


#endif
