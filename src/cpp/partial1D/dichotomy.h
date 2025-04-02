#pragma once

#include <tl/support/ASSERT.h>
#include <tl/support/TODO.h>

template<class TF>
TF dichotomy( auto &&func, TF tol, TF beg_x, TF end_x ) {
    //
    TF beg_y = func( beg_x );
    if ( abs( beg_y ) < tol )
        return beg_x;

    TF end_y = func( end_x );
    if ( abs( end_y ) < tol )
        return end_x;

    ASSERT( ( beg_y > 0 ) != ( end_y > 0 ) );

    //
    for( int num_iter = 0; ; ++num_iter ) {
        const TF mid_x = ( beg_x * end_y - end_x * beg_y ) / ( end_y - beg_y );
        if ( num_iter == 50 )
            return mid_x;

        const TF mid_y = func( mid_x );
        if ( abs( mid_y ) < tol || beg_x == end_x )
            return mid_x;

        if ( mid_y < 0 ) {
            beg_x = mid_x;
            beg_y = mid_y;
        } else {
            end_x = mid_x;
            end_y = mid_y;
        }
    }
}

template<class TF>
TF dichotomy( auto &&func, TF tol, TF mid ) {
    // // beg, end
    // TF beg_x = prev_x, beg_y = mid_y;
    // TF end_x = prev_x, end_y = mid_y;
    // if ( mid_y < 0 ) {
    //     do {
    //         end_x *= 2;
    //         end_y = integral( center - end_x, center + end_x ) - target_mass;
    //     } while ( end_y < 0 );
    // } else {
    //     do {
    //         beg_x /= 2;
    //         beg_y = integral( center - beg_x, center + beg_x ) - target_mass;
    //     } while ( beg_y > 0 );
    // }
    TODO;
}