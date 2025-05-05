#pragma once

#include "common_macros.h"
#include <vector>

namespace usdot {

/**
 * Tridiagonal symmetric matrix + 1 line (that can be used for instance for a lagrange multiplier)
 * 
*/
template<class T>
class TridiagonalPlus1LineSymmetricMatrix {
public:
    using    TV                                 = std::vector<T>;
    using    PI                                 = std::size_t;
        
    /**/     TridiagonalPlus1LineSymmetricMatrix( PI tridiag_size );
    /**/     TridiagonalPlus1LineSymmetricMatrix();

    void     clear_tridiag_values               ();
               
    const T& tridiag_value                      ( PI r, PI c ) const;
    T&       tridiag_value                      ( PI r, PI c );
    const T& line_value                         ( PI c ) const; ///< 
    T&       line_value                         ( PI c ); ///< 

    PI       full_matrix_size                   () const { return tridiag_size() + 1; }
    PI       tridiag_size                       () const;
                
    T        secure_get                         ( PI r, PI c ) const;
    T_U void display                            ( U &ds ) const;
                   
    void     inplace_ldlt_decomposition         ();
    TV       solve_using_ldlt                   ( const TV &tridiag_rhs, T line_rhs ) const;
         
    TV       tridiag_values;                    ///< M[0,0], M[0,1], M[1,1], M[1,2], ...
    TV       line_values;                       ///< 
};

#define DTP template<class T>
#define UTP TridiagonalPlus1LineSymmetricMatrix<T>

DTP UTP::TridiagonalPlus1LineSymmetricMatrix( PI tridiag_size ) : tridiag_values( tridiag_size ? tridiag_size * 2 - 1 : 0 ), line_values( tridiag_size + 1 ) {
}

DTP UTP::TridiagonalPlus1LineSymmetricMatrix() {
}

DTP T UTP::secure_get( PI r, PI c ) const {
    if ( r > tridiag_size() || c > tridiag_size() )
        return 0;

    if ( r < tridiag_size() && c < tridiag_size() ) {
        if ( r <= c + 1 && c <= r + 1 )
            return tridiag_value( r, c );
        return 0;
    }

    return line_values[ std::min( r, c ) ];
}

DTP const T& UTP::tridiag_value( PI r, PI c ) const {
    return tridiag_values[ r + c ];
}

DTP T &UTP::tridiag_value( PI r, PI c ) {
    return tridiag_values[ r + c ]; 
}

DTP const T& UTP::line_value( PI c ) const {
    return line_values[ c ];
}

DTP T &UTP::line_value( PI c ) {
    return line_values[ c ]; 
}

DTP T_U void UTP::display( U &ds ) const {
    ds.start_array();
    for( PI r = 0; r < full_matrix_size(); ++r ) {
        ds.start_array();
        for( PI c = 0; c < full_matrix_size(); ++c )
            ds << secure_get( r, c );
        ds.end_array();
    }
    ds.end_array();
}

DTP typename UTP::PI UTP::tridiag_size() const {
    return ( tridiag_values.size() + 1 ) / 2;
}

DTP void UTP::clear_tridiag_values() {
    std::fill( tridiag_values.begin(), tridiag_values.end(), 0 );
}

template<class T>
void TridiagonalPlus1LineSymmetricMatrix<T>::inplace_ldlt_decomposition() {
    const PI s = tridiag_size();
    if ( s == 0 )
        return;

    for( PI i = 1; i < tridiag_values.size(); i += 2 ) {
        const T p = tridiag_values[ i - 1 ];
        const T l = tridiag_values[ i - 0 ] / p;

        tridiag_values[ i + 0 ] = l;
        tridiag_values[ i + 1 ] -= l * l * p;
    }

    T d = line_values[ s ];
    for( PI c = 0; c < s; ++c ) {
        T p = line_values[ c ];
        if ( c )
            p -= tridiag_values[ 2 * c - 1 ] * tridiag_values[ 2 * c - 2 ] * line_values[ c - 1 ];
        p /= tridiag_values[ 2 * c ];

        d -= p * p * tridiag_values[ 2 * c ];
        line_values[ c ] = p;
    }
    line_values[ s ] = d;
}

template<class T>
inline typename UTP::TV UTP::solve_using_ldlt( const TV &x, T lag_rhs ) const {
    using std::pow;

    const PI n = tridiag_size();
    if ( n == 0 )
        return {};

    // solve with L
    TV res( n );
    res[ 0 ] = x[ 0 ];
    for( PI i = 1; i < n; ++i )
        res[ i ] = x[ i ] - res[ i - 1 ] * tridiag_values[ 2 * i - 1 ];
    for( PI i = 0; i < n; ++i )
        lag_rhs -= res[ i ] * line_values[ i ];

    // solve with D and Lt
    T l = lag_rhs / line_values[ n ];
    for( PI i = n; i--; )
        res[ i ] = res[ i ] / tridiag_values[ 2 * i ] - res[ i + 1 ] * tridiag_values[ 2 * i + 1 ] - line_values[ i ] * l;

    return res;
}

#undef DTP
#undef UTP

} // namespace usdot
