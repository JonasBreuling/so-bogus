// #include "../src/Core/Block.impl.hpp"
#include "Core/Block.impl.hpp"
#include "Core/Utils/Signal.hpp"

#include "Core/Block.impl.hpp"
#include "Core/BlockSolvers/GaussSeidel.impl.hpp"
#include "Core/BlockSolvers/LCPLaw.impl.hpp"
#include "Core/BlockSolvers/PyramidLaw.impl.hpp"

#include "Extra/SecondOrder.impl.hpp"

#include <Eigen/LU>
#include <Eigen/Cholesky>

#include <string>
#include <iostream>

class ResidualInfo {

public:
	explicit ResidualInfo( bool verbose = false )
		: m_verbose( verbose )
	{ }

	void setVerbose( bool verbose ) {
		m_verbose = verbose ;
	}

	void ack( unsigned iter, double err ) {
		if( m_verbose ) {
			std::cout << m_meth << ": \t" << iter << "\t ==> " << err << std::endl ;
		}
	}

	void setMethodName( const std::string& meth ) {
		m_meth = meth ;
	}

	void bindTo( bogus::Signal<unsigned, double> &signal ) {
		signal.connect( *this, &ResidualInfo::ack );
	}

private:

	bool m_verbose ;
	std::string m_meth ;

};

class SmallFrictionPb {
public:

	SmallFrictionPb() {

		const unsigned dofs[2] = { 4, 2 } ;

		MassMat.setRows( 2, dofs ) ;
		MassMat.setCols( 2, dofs ) ;

		MassMat.insertBackAndResize( 0, 0 ) << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 ;
		MassMat.insertBackAndResize( 1 ,1 ) << 2, 0, 0, 2 ;

		MassMat.finalize() ;

		InvMassMat.cloneStructure( MassMat ) ;
		for( unsigned i = 0 ; i < MassMat.nBlocks() ; ++i )
		{
			InvMassMat.block( i ) = MassMat.block( i ).inverse() ;
		}


		H.setCols( 2, dofs ) ;
		H.setRows( 2 ) ;
		H.insertBackAndResize( 0, 0 ) << 3, 5, 7, -1, 2, 4, 5, 4, 1, 6, 9, 1 ;
		H.insertBackAndResize( 0, 1 ) << -3, -6, -6, -5, -5, -8 ;
		H.insertBackAndResize( 1, 0 ) << 3, 7, 1, -1, 3, 6, 3, 2, 4, 4, 7, -1 ;

		GradBlockT H2B( 3, dofs[0] ) ;
		H2B << 6, 7, 10, 9, 5, 8, 11, 8, 4, 9, 12, 7 ;
		H.block( 2 ) -= H2B ;

		Eigen::Matrix3d E ;
		E << 0, 0, 1, 1, 0, 0, 0, 1, 0 ;

		H.block( 2 ) = E.transpose() * H.block( 2 ) ;
		H.finalize() ;

		//std::cout << MassMat << std::endl ;
		//	std::cout << InvMassMat << std::endl ;
		//	std::cout << H << std::endl ;

		//	std::cout << W << std::endl ;
		//	W.cacheTranspose();

		f.resize( 6 ) ;
		f << 1, 2, 3, 4, 5, 6 ;
		w.resize( 6 ) ;
		w << 2, 1, 2, 1, 3, 3 ;

		mu.resize(2) ;
		mu <<  0.5, 0.7 ;

		sol.resize( 6 ) ;
		sol << 0.0152695, 0.0073010, 0.0022325, 0.0, 0.0, 0.0 ;
	}

	typedef bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::UNCOMPRESSED > MType ;
	typedef Eigen::Matrix< double, 3, Eigen::Dynamic > GradBlockT ;
	typedef bogus::SparseBlockMatrix< GradBlockT > HType ;

	MType MassMat ;
	MType InvMassMat ;

	HType H ;
	Eigen::VectorXd f ;
	Eigen::VectorXd w ;

	Eigen::VectorXd mu ;

	Eigen::VectorXd sol ;

	// virtual void TearDown() {}

};

// int main(int argc, char** argv) {
int main(int, char**) {

    SmallFrictionPb prob;

    {
        ResidualInfo ri ;

        Eigen::VectorXd b = prob.w - prob.H * ( prob.InvMassMat * prob.f );


        typedef bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::SYMMETRIC > WType;
        WType W = prob.H * prob.InvMassMat * prob.H.transpose();


        Eigen::VectorXd x( W.rows() ) ;
        double res = -1 ;

        bogus::GaussSeidel< WType > gs( W );
        ri.bindTo( gs.callback() );

        x.setOnes() ;
        ri.setMethodName( "GS_Hyb" );
        res = gs.solve( bogus::SOCLaw< 3u, double, true, bogus::local_soc_solver::Hybrid >( 2, prob.mu.data() ), b, x ) ;

        std::cout << "res: " << res << std::endl;
        std::cout << "x: " << x.transpose() << std::endl;
        std::cout << "sol: " << prob.sol.transpose() << std::endl;
        // ASSERT_LT( res, 1.e-8 ) ;
        // ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;

        // x.setOnes() ;
        // ri.setMethodName( "GS_PureN" );
        // res = gs.solve( bogus::SOCLaw< 3u, double, true, bogus::local_soc_solver::PureNewton >( 2, mu.data() ), b, x ) ;
        // ASSERT_LT( res, 1.e-8 ) ;
        // ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;

        // x.setOnes() ;
        // ri.setMethodName( "GS_PureE" );
        // res = gs.solve( bogus::SOCLaw< 3u, double, true, bogus::local_soc_solver::PureEnumerative >( 2, mu.data() ), b, x ) ;
        // ASSERT_LT( res, 1.e-8 ) ;
        // ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;

        // x.setOnes() ;
        // ri.setMethodName( "GS_RevHyb" );
        // res = gs.solve( bogus::SOCLaw< 3u, double, true, bogus::local_soc_solver::RevHybrid >( 2, mu.data() ), b, x ) ;
        // ASSERT_LT( res, 1.e-8 ) ;
        // ASSERT_TRUE( sol.isApprox( x, 1.e-4 ) ) ;

        // x.setOnes() ;
        // ri.setMethodName( "GS_Def" );
        // gs.setTol( 1.e-8 );
        // res = gs.solve( bogus::SOC3D( 2, mu.data() ), b, x ) ;
        // ASSERT_LT( res, 1.e-8 ) ;
    }
    return 0;
}