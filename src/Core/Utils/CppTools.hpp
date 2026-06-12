/*
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * Copyright 2013 Gilles Daviet <gdaviet@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOGUS_CPP_TOOLS_HPP
#define BOGUS_CPP_TOOLS_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include <type_traits>

namespace bogus
{

#ifndef BOGUS_HAS_CPP11
#define BOGUS_HAS_CPP11 (__cplusplus >= 201103L)
#endif

template<typename T, typename = void>
struct HasReturnType : std::false_type {};

template<typename T>
struct HasReturnType<T, std::void_t<typename T::ReturnType>>
    : std::true_type {};

template<typename T, typename = void>
struct HasConstTransposeReturnType : std::false_type {};

template<typename T>
struct HasConstTransposeReturnType<
    T,
    std::void_t<typename T::ConstTransposeReturnType>>
    : std::true_type {};

// TODO: This is never used!
template<typename T, typename = void>
struct HasBase : std::false_type {};

template<typename T>
struct HasBase<T, std::void_t<typename T::Base>>
    : std::true_type {};

// Static assertions

template < bool Assertion >
struct StaticAssert
{
	enum {
		BLOCKS_MUST_BE_SQUARE_OR_HAVE_DYNAMIC_DIMENSIONS,
		BLOCKS_MUST_HAVE_FIXED_DIMENSIONS,
		MATRICES_ORDERING_IS_INCONSISTENT,
		TRANSPOSE_MAKES_NO_SENSE_IN_THIS_CONTEXT,
		TRANSPOSE_IS_NOT_DEFINED_FOR_THIS_BLOCK_TYPE,
		OPERANDS_HAVE_INCONSISTENT_FLAGS,
		UNORDERED_INSERTION_WITH_COMPRESSED_INDEX,
		NOT_IMPLEMENTED
	} ;
} ;

template < >
struct StaticAssert< false >
{
} ;

#define BOGUS_STATIC_ASSERT( test, message ) (void) ::bogus::StaticAssert< test >::message

//! Const mapped array, used for Mapped Block Matrices
template< typename Element >
class ConstMappedArray
{
public:

	typedef ConstMappedArray< Element > Type ;
	enum { is_mutable = 0 } ;

	ConstMappedArray ( )
	    : m_data( 0 ), m_size( 0 )
	{}

	ConstMappedArray ( const Element* data, std::size_t size )
	    : m_data( data ), m_size( size )
	{}

	void setData( const Element* data, std::size_t size )
	{
		m_data = data ;
		m_size = size ;
	}

	const Element* data() const { return m_data ; }

	inline std::size_t size() const { return m_size ; }
	inline bool empty() const { return 0 == m_size ; }

	const Element* begin() const { return data() ; }
	const Element* end() const { return data() + size() ; }

	const Element& operator[]( std::size_t idx ) const
	{ return m_data[idx] ; }

	inline void resize( std::size_t s)
	{ assert( !m_data || m_size == s ) ; (void) s ; }
	inline void reserve( std::size_t )
	{ }
	inline void clear( )
	{ }
	inline void assign( std::size_t s, const Element&)
	{ assert( !m_data || m_size == s ) ; (void) s ; }

private:
	const Element* m_data ;
	std::size_t    m_size ;
} ;

template< typename Scalar>
struct ConstantArray {
	const Scalar s ;
	explicit ConstantArray( Scalar s_=1 ) : s(s_) {}

	inline bool  has_element( int ) const { return true ; }
	inline Scalar    element( int ) const { return s ; }
	inline Scalar operator[]( int ) const { return s ; }
};
template< typename Scalar>
inline ConstantArray<Scalar> make_constant_array( Scalar s) { return ConstantArray<Scalar>(s) ; }

} //namespace bogus


#endif



