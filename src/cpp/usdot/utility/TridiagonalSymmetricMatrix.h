#pragma once

#include "common_macros.h"
#include <vector>

namespace usdot {

/**
 * Tridiagonal symmetric matrix + optional room for lines (like lagrange multipliers)
 * 
*/
template<class T>
class TridiagonalSymmetricMatrix {
public:
    using    PI                        = std::size_t;
    using    TV                        = std::vector<T>;
        
    /**/     TridiagonalSymmetricMatrix( PI size, PI nb_lags, T default_bands_value, T default_lags_value );
    /**/     TridiagonalSymmetricMatrix( PI size, PI nb_lags = 0 );
    /**/     TridiagonalSymmetricMatrix();

    void     clear_tridiag_values      ();
    void     set_lag_values            ( T value );
      
    const T& tridiag_value             ( PI r, PI c ) const;
    T&       tridiag_value             ( PI r, PI c );
    const T& lag_value                 ( PI l, PI u ) const; ///< l = num multiplier, u = num unknown
    T&       lag_value                 ( PI l, PI u ); ///< l = num multiplier, u = num unknown
       
    PI       full_matrix_size          () const { return tridiag_size() + nb_lags(); }
    PI       tridiag_size              () const;
    PI       nb_lags                   () const;
       
    T        secure_get                ( PI r, PI c ) const;
    T_U void display                   ( U &ds ) const;
          
    void     inplace_ldlt_decomposition();
    TV       solve_using_ldlt          ( const TV &x, const TV &lag_rhs ) const;

    TV       tridiag_values;           ///< M[0,0], M[0,1], M[1,1], M[1,2], ...
    TV       lag_values;
};

#define DTP template<class T>
#define UTP TridiagonalSymmetricMatrix<T>

DTP UTP::TridiagonalSymmetricMatrix( PI tridiag_size, PI nb_lags, T default_bands_value, T default_lags_value ) : tridiag_values( tridiag_size ? tridiag_size * 2 - 1 : 0, default_bands_value ), lag_values( ( tridiag_size + nb_lags ) * nb_lags, default_lags_value ) {
    if ( default_lags_value )
        for( PI r = 0; r < nb_lags; ++r )
            for( PI c = tridiag_size; c < tridiag_size + nb_lags; ++c )
                lag_value( r, c ) = 0;
}

DTP UTP::TridiagonalSymmetricMatrix( PI tridiag_size, PI nb_lags ) : tridiag_values( tridiag_size ? tridiag_size * 2 - 1 : 0 ), lag_values( tridiag_size * nb_lags ) {
    for( PI r = 0; r < nb_lags; ++r )
        for( PI c = tridiag_size; c < tridiag_size + nb_lags; ++c )
            lag_value( r, c ) = 0;
}

DTP UTP::TridiagonalSymmetricMatrix() {
}

DTP T UTP::secure_get( PI r, PI c ) const {
    using namespace std; 
    if ( r < tridiag_size() && c < tridiag_size() ) {
        if ( r <= c + 1 && c <= r + 1 )
            return tridiag_value( r, c );
        return 0;
    }
    if ( c < full_matrix_size() && r - tridiag_size() < nb_lags() )
        return lag_values[ ( r - tridiag_size() ) * full_matrix_size() + c ];
    if ( r < full_matrix_size() && c - tridiag_size() < nb_lags() )
        return lag_values[ ( c - tridiag_size() ) * full_matrix_size() + r ];
    return 0;
}

DTP const T& UTP::tridiag_value( PI r, PI c ) const {
    return tridiag_values[ r + c ];
}

DTP T &UTP::tridiag_value( PI r, PI c ) {
    return tridiag_values[ r + c ]; 
}

DTP const T& UTP::lag_value( PI l, PI c ) const {
    return lag_values[ l * full_matrix_size() + c ];
}

DTP T &UTP::lag_value( PI l, PI c ) {
    return lag_values[ l * full_matrix_size() + c ]; 
}

DTP T_U void UTP::display( U &ds ) const {
    ds.start_array();
    for( PI r = 0; r < tridiag_size() + nb_lags(); ++r ) {
        ds.start_array();
        for( PI c = 0; c < tridiag_size() + nb_lags(); ++c )
            ds << secure_get( r, c );
        ds.end_array();
    }
    ds.end_array();
}

DTP UTP::PI UTP::tridiag_size() const {
    return ( tridiag_values.size() + 1 ) / 2;
}

DTP UTP::PI UTP::nb_lags() const {
    return lag_values.size() / tridiag_size();
}

DTP void UTP::clear_tridiag_values() {
    std::fill( tridiag_values.begin(), tridiag_values.end(), 0 );
}

DTP void UTP::set_lag_values( T value ) {
    std::fill( lag_values.begin(), lag_values.end(), value );
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
void TridiagonalSymmetricMatrix<T>::inplace_ldlt_decomposition() {
    const PI t = tridiag_size();
    const PI l = nb_lags();
    if ( t == 0 )
        return;

    for( PI i = 1; i < tridiag_values.size(); i += 2 ) {
        const T p = tridiag_values[ i - 1 ];
        const T l = tridiag_values[ i - 0 ] / p;

        tridiag_values[ i + 0 ] = l;
        tridiag_values[ i + 1 ] -= l * l * p;
    }

    for( PI r = 0; r < l; ++r ) {
        T s = 0;
        for( PI c = 0; c < t; ++c ) {
            T p = lag_value( r, c );
            if ( c )
                p -= tridiag_values[ 2 * c - 1 ] * tridiag_values[ 2 * c - 2 ] * lag_value( r, c - 1 );
            p /= tridiag_values[ 2 * c ];

            s -= p * p * tridiag_values[ 2 * c ];
            lag_value( r, c ) = p;
        }
        lag_value( r, t ) = s;
    }
}

template<class T>
inline UTP::TV UTP::solve_using_ldlt( const TV &x, const TV &lag_rhs ) const {
    using std::pow;

    const PI n = tridiag_size();
    if ( n == 0 )
        return {};

    // solve with L
    TV res( n );
    res[ 0 ] = x[ 0 ];
    for( PI i = 1; i < n; ++i )
        res[ i ] = x[ i ] - res[ i - 1 ] * tridiag_values[ 2 * i - 1 ];

    // solve with D and Lt
    res[ n - 1 ] /= tridiag_values[ 2 * ( n - 1 ) ];
    for( PI i = n - 1; i--; ) {
        // ASSERT( chol.data[ 2 * i ] );
        res[ i ] = res[ i ] / tridiag_values[ 2 * i ] - res[ i + 1 ] * tridiag_values[ 2 * i + 1 ];
    }

    return res;
}

#undef DTP
#undef UTP

} // namespace usdot
