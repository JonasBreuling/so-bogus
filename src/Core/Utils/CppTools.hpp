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

namespace bogus {

// Type-trait helpers

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

// Static assertions
#define BOGUS_STATIC_ASSERT( test, message ) \
    static_assert( (test), #message )

//! Const mapped array (used for Mapped Block Matrices)
template< typename Element >
class ConstMappedArray
{
public:

	using Type = ConstMappedArray< Element >;
	static constexpr bool is_mutable = false;

	ConstMappedArray() : m_data( nullptr ), m_size( 0 ) {}

	ConstMappedArray ( const Element* data, std::size_t size )
	    : m_data( data ), m_size( size ) {}

	void setData( const Element* data, std::size_t size )
	{
		m_data = data ;
		m_size = size ;
	}

	const Element*     data()  const noexcept { return m_data ; }
	inline std::size_t size()  const noexcept { return m_size ; }
	inline bool        empty() const noexcept { return m_size == 0 ; }

	const Element* begin() const noexcept { return m_data ; }
	const Element* end()   const noexcept { return m_data + m_size ; }

	const Element& operator[]( std::size_t idx ) const { return m_data[idx]; }

	inline void resize( std::size_t s) { assert( !m_data || m_size == s ) ; (void)s ; }
	inline void reserve( std::size_t ) {}
	inline void clear() {}
	inline void assign( std::size_t s, const Element&) { assert( !m_data || m_size == s ) ; (void)s ; }

private:
	const Element* m_data;
	std::size_t    m_size;
} ;

//! Constant-array helper (avoids heap allocation for scalar scaling)
template< typename Scalar>
struct ConstantArray
{
	const Scalar s ;
	explicit constexpr ConstantArray( Scalar s_ = 1 ) : s(s_) {}

	inline constexpr bool  has_element( int ) const noexcept { return true ; }
	inline constexpr Scalar    element( int ) const noexcept { return s ; }
	inline constexpr Scalar operator[]( int ) const noexcept { return s ; }
};
template< typename Scalar>
inline ConstantArray<Scalar> make_constant_array( Scalar s) { return ConstantArray<Scalar>(s) ; }

} //namespace bogus


#endif



