#pragma once

#include <tl/support/operators/sgn.h>
#include <tl/support/ASSERT.h>
#include <tl/support/TODO.h>
#include <stdexcept>
#include <cmath>

template<class TF>
TF dichotomy( auto &&func, TF tol, TF beg_x, TF end_x, TF beg_y, TF end_y ) {
    ASSERT( sgn( beg_y ) != sgn( end_y ) );

    //
    for( int num_iter = 0; ; ++num_iter ) {
        const TF mid_x = ( beg_x * end_y - end_x * beg_y ) / ( end_y - beg_y );
        if ( num_iter == 50000 )
            throw std::runtime_error( "dichotomy pb" );

        const TF mid_y = func( mid_x );
        if ( abs( mid_y ) < tol || beg_x == end_x )
            return mid_x;

        if ( ( mid_y < 0 ) ^ ( beg_y > end_y ) ) {
            beg_x = mid_x;
            beg_y = mid_y;
        } else {
            end_x = mid_x;
            end_y = mid_y;
        }
    }

    return {};
}

template<class TF>
TF dichotomy( auto &&func, TF tol, TF beg_x, TF end_x ) {
    //
    TF beg_y = func( beg_x );
    if ( abs( beg_y ) < tol )
        return beg_x;
    
    TF end_y = func( end_x );
    if ( abs( end_y ) < tol )
        return end_x;

    return dichotomy( func, tol, beg_x, end_x, beg_y, end_y );
}

template<class TF>
TF dichotomy_growing_from_zero( auto &&func, TF tol, TF mid ) {
    // beg, end
    TF beg_x = mid, beg_y = func( mid );
    TF end_x = mid, end_y = beg_y;
    if ( end_y < 0 ) {
        do {
            end_y = func( end_x *= 2 );
        } while ( end_y < 0 );
    } else {
        do {
            beg_y = func( beg_x /= 2 );
        } while ( beg_y > 0 );
    }

    return dichotomy( func, tol, beg_x, end_x, beg_y, end_y );
}
