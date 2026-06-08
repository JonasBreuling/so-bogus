/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_NoneLAW_HPP
#define BOGUS_NoneLAW_HPP

#include <cmath>

namespace bogus
{

//! Dummy local solver that can be used within GaussSeidel and ProjectedGradient solvers
/*!
  For testing purposes only.
	\tparam Scalar the scalar type
	\tparam Dimension the dimension of the blocks of the global matrix
  */
template < typename Scalar >
class NoneLaw
{
public:
	enum{ dimension = 1 } ;

	typedef LocalProblemTraits< dimension, Scalar > Traits ;

	//! Constructor
	NoneLaw( ) {}

	//! Projects x on \f$ R \f$. This is a none operation.
	void projectOnConstraint( const unsigned problemIndex, typename Traits::Vector &x ) const {};

	//! \return \f$ 0 \f$
	Scalar eval( const unsigned problemIndex,
	             const typename Traits::Vector &x,
	             const typename Traits::Vector &y ) const
	{
		return (Scalar) 0;
	}

	//! Solves the local problem (which is always solved).
	bool solveLocal(
	        const unsigned problemIndex,
	        const typename Traits::Matrix &A,
	        const typename Traits::Vector &b,
	        typename Traits::Vector &x,
	        const Scalar scaling
	        ) const
	{
		return true;
	}

	//! This NSLaw is always associated, so dualityCOV is null.
	template< typename Segment >
	void dualityCOV( const unsigned , const Segment& ,
	                 typename Traits::Vector& s ) const
	{ }
	// { s->setZero() ; }
} ;


}

#endif
