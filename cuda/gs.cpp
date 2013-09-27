
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

#include <Interfaces/MecheInterface.hpp>
#include <Interfaces/FrictionProblem.hpp>
#include <Interfaces/FrictionProblem.impl.hpp>

#include <Core/Utils/Timer.hpp>

#include "loaders.hpp"

#include "gs.hpp"
#include "gs.impl.hpp"

void ackCurrentResidual( unsigned GSIter, double err )
{
	std::cout << "Finished iteration " << GSIter
				 << " with residual " << err
				 << std::endl ;
}

int main( int argc, const char* argv[] )
{
  const unsigned D = 3 ;

  if( argc < 2 )
  {
    std::cerr << " Please provide a problem data file " << std::endl ;
    std::cerr << " Syntax: " << argv[0] << " dataFile "
      << " [ maxThreads ] [ tol ] [ maxIters ] [ staticPb ] [ regul ] [ useInfNorm ] "
      << std::endl ;
    return 1 ;
  }
  
  //Eigen::initParallel()
  Eigen::MatrixXf A = Eigen::MatrixXf::Zero(1,1); A = A*A;

  bogus::MecheFrictionProblem mfp ;

  double * r = NULL ;
  if(  ! mfp.fromFile( argv[1], r ) ) 
    return 1 ;


  const int maxThreads    = argc > 2 ? std::atoi( argv[2] ) : 0 ;
  const double tol        = argc > 3 ? std::strtod( argv[3], NULL ) : 0 ;
  const int maxIters      = argc > 4 ? std::atoi( argv[4] ) : 0 ;
  const int staticPb      = argc > 5 ? std::atoi( argv[5] ) : 0 ;
  const double regul      = argc > 6 ? std::strtod( argv[6], NULL ) : 0 ;
  const int useInfNorm    = argc > 7 ? std::atoi( argv[7] ) : 0 ;

  std::cout << "Computing dual.." << std::endl ;
  mfp.computeDual( regul ) ;
  
  mfp.dual().W.cacheTranspose() ;
  
  const int n = mfp.dual().W.rowsOfBlocks() ;
  Eigen::VectorXd rhs( D*n ), res( D*n ) ;
  rhs.setOnes() ;

/* Matrix-vector prodcut tests
 *
 
  bogus::Timer timer ;
  Eigen::VectorXd cpuRes = mfp.dual().W*rhs ;
  std::cout << timer.elapsed() << " s" << std::endl ;
  
  bogus::CudaContext ctx ;
  typedef bogus::CudaMatrix< double, int, D > CuMat ;
  typedef typename CuMat::Vector CuVec ;

  CuMat mat( ctx ), tmat( ctx ) ;

#ifdef COMBINE
  std::cout << "Recombining matrix.." << std::endl ;
  std::vector<int> ri,ci;
  load_combined( mfp.dual().W, ri, ci, mat ) ;
#else
  std::cout << "Reordering matrix.." << std::endl ;
  load( mfp.dual().W, mat ) ;
  std::cout << "Making transpose.." << std::endl ;
  load_transposed( mfp.dual().W, tmat ) ;
#endif
  
  std::cout << "Creating rhs and res vecs" << std::endl ;

  CuVec cuRhs( ctx, n ), cuRes( ctx, n ) ;

  cuRhs.set( rhs.data() ) ;
  cuRes.reset() ;

  std::cout << "Testing mv kernel.." << std::endl ;
  timer.reset() ;
  mat.mult( cuRhs, cuRes ) ;
  tmat.mult( cuRhs, cuRes ) ;
  std::cout << timer.elapsed() << " s" << std::endl ;
  
  cuRes.get( res.data() ) ;
  std::cout << timer.elapsed() << " s" << std::endl ;
  
  std::cout << "Done" << std::endl ;
  
  
  std::cout << (cpuRes - res).squaredNorm() << std::endl ;
  
  timer.reset() ;
#pragma omp parallel for
  for( unsigned i = 0 ; i < n ; ++i )
  {
    Eigen::Vector3d res3 ;
     mfp.dual().W.splitRowMultiply( i, rhs, res3 ) ;
  }
  std::cout << timer.elapsed() << " s" << std::endl ;

 */

  bogus::CuGaussSeidel< typename bogus::DualFrictionProblem<D>::WType > cugs ;

  if( tol != 0. ) cugs.setTol( tol );
  if( maxIters != 0 ) cugs.setMaxIters( maxIters );

  cugs.setMaxThreads( maxThreads ) ;
  cugs.setAutoRegularization( regul ) ;
  cugs.useInfinityNorm( useInfNorm ) ;
  cugs.setSkipTol( std::sqrt( cugs.tol() ) ) ;


  cugs.callback().connect( &ackCurrentResidual ) ;
  cugs.setMatrix( mfp.dual().W ) ;
  bogus::friction_problem::solve( mfp.dual(), cugs, rhs.data(), staticPb ) ;

  return 0 ;

}

