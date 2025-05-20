#pragma once

#include "common_macros.h"
#include <vector>
#include <cmath>

namespace usdot {

/**
 * Tridiagonal symmetric matrix
 * 
*/
template<class T>
class TridiagonalSymmetricMatrix {
public:
    using    PI                        = std::size_t;
    using    TV                        = std::vector<T>;
        
    /**/     TridiagonalSymmetricMatrix( PI size, T default_bands_value );
    /**/     TridiagonalSymmetricMatrix( PI size );
    /**/     TridiagonalSymmetricMatrix();

    void     clear_values              ();
    void     resize                    ( PI size );
    
    T        secure_get                ( PI r, PI c ) const;
    const T& operator()                ( PI r, PI c ) const;
    T&       operator()                ( PI r, PI c );
    
    T_U void display                   ( U &ds ) const;
    PI       size                      () const;
          
    int      inplace_ldlt_decomposition();
    TV       solve_using_ldlt          ( const TV &x ) const;

    TV       values;                   ///< M[0,0], M[0,1], M[1,1], M[1,2], ...
};

#define DTP template<class T>
#define UTP TridiagonalSymmetricMatrix<T>

DTP UTP::TridiagonalSymmetricMatrix( PI size, T default_value ) : values( size ? size * 2 - 1 : 0, default_value ) {
}

DTP UTP::TridiagonalSymmetricMatrix( PI size ) : values( size ? size * 2 - 1 : 0 ) {
}

DTP UTP::TridiagonalSymmetricMatrix() {
}

DTP void UTP::resize( PI size ) {
    values.resize( size ? size * 2 - 1 : 0 );
}

DTP T UTP::secure_get( PI r, PI c ) const {
    using namespace std; 
    if ( r <= c + 1 && c <= r + 1 )
        return operator()( r, c );
    return 0;
}

DTP const T& UTP::operator()( PI r, PI c ) const {
    return values[ r + c ];
}

DTP T &UTP::operator()( PI r, PI c ) {
    return values[ r + c ]; 
}

DTP T_U void UTP::display( U &ds ) const {
    ds.start_array();
    for( PI r = 0; r < size(); ++r ) {
        ds.start_array();
        for( PI c = 0; c < size(); ++c )
            ds << secure_get( r, c );
        ds.end_array();
    }
    ds.end_array();
}

DTP typename UTP::PI UTP::size() const {
    return ( values.size() + 1 ) / 2;
}

DTP void UTP::clear_values() {
    std::fill( values.begin(), values.end(), 0 );
}

// DTP UTP::TV UTP::SymmetricBandMatrix<T>::mul( const Vec<T> &x ) const {
//     if ( nb_rows() == 0 )
//         return {};

//     if ( nb_rows() == 1 )
//         return { data[ 0 ] * x[ 0 ] };

//     Vec<T> res( FromSize(), nb_rows() );
//     const PI l = nb_rows() - 1;

//     res[ 0 ] = data[ 0 ] * x[ 0 ] + data[ 1 ] * x[ 1 ];
//     for( PI r = 1; r < l; ++r )
//         res[ r ] = data[ 2 * r - 1 ] * x[ r - 1 ] + data[ 2 * r + 0 ] * x[ r + 0 ] + data[ 2 * r + 1 ] * x[ r + 1 ];
//     res[ l ] = data[ 2 * l - 1 ] * x[ l - 1 ] + data[ 2 * l + 0 ] * x[ l + 0 ];

//     return res;
// }

template<class T>
int TridiagonalSymmetricMatrix<T>::inplace_ldlt_decomposition() {
    for( PI i = 1; i < values.size(); i += 2 ) {
        const T p = values[ i - 1 ];
        if ( p == 0 )
            return 1;
        const T l = values[ i - 0 ] / p;
        
        values[ i + 0 ] = l;
        values[ i + 1 ] -= l * l * p;
    }
    if ( values.back() == 0 )
        return 1;
    return 0;
}

template<class T>
inline typename UTP::TV UTP::solve_using_ldlt( const TV &x ) const {
    using std::pow;

    const PI n = size();
    if ( n == 0 )
        return {};

    // solve with L
    TV res( n );
    res[ 0 ] = x[ 0 ];
    for( PI i = 1; i < n; ++i )
        res[ i ] = x[ i ] - res[ i - 1 ] * values[ 2 * i - 1 ];

    // solve with D and Lt
    res[ n - 1 ] /= values[ 2 * ( n - 1 ) ];
    for( PI i = n - 1; i--; )
        res[ i ] = res[ i ] / values[ 2 * i ] - res[ i + 1 ] * values[ 2 * i + 1 ];

    return res;
}

#undef DTP
#undef UTP

} // namespace usdot
