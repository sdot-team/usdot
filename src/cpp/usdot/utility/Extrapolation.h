#pragma once
 
#include <tl/support/containers/Vec.h>
#include <eigen3/Eigen/Dense>

namespace usdot {

/// assuming ws = value at 0, 1, ..., ws.size() - 1
template<class TF>
Vec<TF> extrapolation( Vec<Vec<TF>> &ws, TF x, int n = -1 ) {
    if ( n < 0 )
        n = ws.size();

    using TM = Eigen::Matrix<TF,Eigen::Dynamic,Eigen::Dynamic>;
    using TV = Eigen::Matrix<TF,Eigen::Dynamic,1>;
    using std::pow;

    TM M( n, n );
    TV V( n );
    for( int r = 0; r < n; ++r )
        for( int c = 0; c < n; ++c )
            M.coeffRef( r, c ) = pow( TF( r ), c );

    Vec<TF> res( FromReservationSize(), ws.size() );
    Eigen::FullPivLU<TM> lu( M );
    for( PI i = 0; i < ws[ 0 ].size(); ++i ) {
        for( int c = 0; c < n; ++c )
            V[ c ] = ws[ c ][ i ];
        auto R = lu.solve( V );

        TF v = 0;
        for( int c = 0; c < n; ++c )
            v += pow( x, c ) * R[ c ];
        res << v;

    }

    return res;
}

/// assuming ws = value at 0, 1, ..., ws.size() - 1
template<class TX,class TF>
TF extrapolation( const TX &ws, TF x, int n = -1 ) {
    if ( n < 0 )
        n = ws.size();

    using TM = Eigen::Matrix<TF,Eigen::Dynamic,Eigen::Dynamic>;
    using TV = Eigen::Matrix<TF,Eigen::Dynamic,1>;
    using std::pow;

    TM M( n, n );
    TV V( n );
    for( int r = 0; r < n; ++r )
        for( int c = 0; c < n; ++c )
            M.coeffRef( r, c ) = pow( TF( r ), c );

    Vec<TF> res( FromReservationSize(), ws.size() );
    Eigen::FullPivLU<TM> lu( M );

    for( int c = 0; c < n; ++c )
        V[ c ] = ws[ c ];
    auto R = lu.solve( V );

    TF v = 0;
    for( int c = 0; c < n; ++c )
        v += pow( x, c ) * R[ c ];

    return v;
}

} // namespace usdot
