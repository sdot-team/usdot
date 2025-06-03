#pragma once

#include <stdexcept>
#include <iostream>
#include <cmath>

namespace usdot {

template<class Func,class TF>
TF dichotomy( const Func &func, TF tol, TF beg_x, TF end_x, TF beg_y, TF end_y, int max_iter = 50000 ) {
    using namespace std;

    if ( beg_y * end_y > 0 )
        return abs( beg_y ) < abs( end_y ) ? beg_x : end_x;

    //
    for( int num_iter = 0; ; ++num_iter ) {
        // const TF mid_x = ( beg_x + end_x ) / 2; // ( beg_x * end_y - end_x * beg_y ) / ( end_y - beg_y );
        const TF f = max( 0.2, min( 0.8, abs( beg_y ) / ( abs( end_y ) + abs( beg_y ) ) ) );
        const TF mid_x = ( 1 - f ) * beg_x + f * end_x;
        const TF mid_y = func( mid_x );
        if ( abs( mid_y ) <= tol || beg_x == end_x )
            return mid_x;

        if ( num_iter == max_iter ) {
            std::cout << " beg_x:" << beg_x << " end_x:" << end_x << " beg_y:" << beg_y << " mid_y:" << mid_y << " end_y:" << end_y << std::endl;
            throw runtime_error( "dichotomy pb" );
        }


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

template<class Func,class TF>
TF dichotomy( const Func &func, TF tol, TF beg_x, TF end_x ) {
    return dichotomy( func, tol, beg_x, end_x, func( beg_x ), func( end_x ) );
}

template<class Func,class TF>
TF dichotomy_unbounded( TF beg_x, TF err_tol, const Func &func ) {
    TF beg_y = func( beg_x );
    if ( abs( beg_y ) <= err_tol )
        return beg_x;

    // bounds
    TF end_x;
    TF end_y;
    for( TF a = 1; ; a *= 2 ) {
        end_x = beg_x + a;
        end_y = func( end_x );
        if ( beg_y * end_y <= 0 )
            break;

        end_x = beg_x - a;
        end_y = func( end_x );
        if ( beg_y * end_y <= 0 )
            break;
    }
        
    return dichotomy( func, err_tol, beg_x, end_x, beg_y, end_y );
}

template<class Func,class TF>
TF dichotomy_growing_from_zero( const Func &func, TF tol, TF mid ) {
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

} // namespace usdot
