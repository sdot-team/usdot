#pragma once

#include <tl/support/containers/Vec.h>
#include <tl/support/Displayer.h>
#include <tl/support/ASSERT.h>
#include <cmath>

/**
 * Symetric diagonal + 1 band
 * 
*/
template<class T>
class TridiagonalMatrix {
public:
    using    TV                 = TL_NAMESPACE::Vec<T>;
        
    /**/     TridiagonalMatrix  ( PI size, T default_value ) : data( FromSizeAndItemValue(), 3 * size, default_value ) {  }
    /**/     TridiagonalMatrix  ( PI size ) : data( FromSize(), 3 * size ) {  }
    /**/     TridiagonalMatrix  () {}

    T        secure_get         ( PI r, PI c ) const { using namespace std; return r < nb_rows() && r < nb_cols() && r <= c + 1 && c <= r + 1 ? operator()( r, c ) : T( 0 ); }
    const T& operator()         ( PI r, PI c ) const { return data[ 2 * r + c ]; }
    T&       operator()         ( PI r, PI c ) { return data[ 2 * r + c ]; }

    TV       diagonal           () const { TV res{ FromSize(), nb_rows() }; for( PI i = 0; i < nb_rows(); ++i ) res[ i ] = operator()( i, i ); return res; }           
    void     display            ( Displayer &ds ) const { ds.start_array(); for( PI r = 0; r < nb_rows(); ++r ) { ds.start_array(); for( PI c = 0; c < nb_cols(); ++c ) ds << secure_get( r, c ); ds.end_array(); } ds.end_array(); }
    PI       nb_rows            () const { return data.size() / 3; }
    PI       nb_cols            () const { return nb_rows(); }
   
    void     fill_with          ( T value ) { data.fill_with( value ); }
    // auto  cholesky           () const -> TridiagonalMatrix;
    // TV    solve              ( const TV &x ) const;
    // TV    mul                ( const TV &x ) const;

    TV       data;
};

// template<class T>
// inline Vec<T> TridiagonalMatrix<T>::mul( const Vec<T> &x ) const {
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

// template<class T>
// TridiagonalMatrix<T> TridiagonalMatrix<T>::cholesky() const {
//     const PI n = nb_rows();
//     if ( n == 0 )
//         return {};

//     TridiagonalMatrix<T> chol( FromSize(), n );
//     chol.data[ 0 ] = data[ 0 ];
//     for( PI i = 1; i < data.size(); i += 2 ) {
//         const T p = chol.data[ i - 1 ];
//         const T l = data[ i - 0 ] / p;

//         chol.data[ i + 0 ] = l;
//         chol.data[ i + 1 ] = data[ i + 1 ] - pow( l, 2 ) * p;
//     }

//     return chol;
// }

// template<class T>
// inline TridiagonalMatrix<T>::TV TridiagonalMatrix<T>::solve( const TV &x ) const {
//     using std::pow;

//     const PI n = nb_rows();
//     if ( n == 0 )
//         return {};

//     // decomposition
//     TridiagonalMatrix chol = cholesky();

//     // solve with L
//     TV res( FromSize(), n );
//     res[ 0 ] = x[ 0 ];
//     for( PI i = 1; i < n; ++i )
//         res[ i ] = x[ i ] - res[ i - 1 ] * chol.data[ 2 * i - 1 ];

//     // solve with D and Lt
//     res[ n - 1 ] /= chol.data[ 2 * ( n - 1 ) ];
//     for( PI i = n - 1; i--; ) {
//         // ASSERT( chol.data[ 2 * i ] );
//         res[ i ] = res[ i ] / chol.data[ 2 * i ] - res[ i + 1 ] * chol.data[ 2 * i + 1 ];
//     }

//     return res;
// }
