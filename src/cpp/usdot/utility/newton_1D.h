#pragma once

#include "linspace.h"
#include "glot.h"
#include <stdexcept>
#include <cmath>

namespace usdot {

template<class TF,class Func>
TF newton_1D( TF old_x, TF min_x, TF max_x, TF error_target, const Func &func, int max_iter = 50000 ) {
    using namespace std;

    // check that we can find the 0
    auto vb = func( min_x );
    auto ve = func( max_x );
    if ( vb.first * ve.first > 0 ) {
        glot( linspace( min_x, max_x, 1000 ), [&]( TF x ) { return func( x ).first; } );
        assert( 0 );
        // auto vm = func( ( min_x + max_x ) / 2 );
        // if ( abs( vm.first ) < abs( vb.first ) && abs( vm.first ) < abs( ve.first ) ) {
        // }

        // return abs( vb.first ) < abs( ve.first ) ? min_x : max_x;
    }

    // first call
    old_x = max( min_x, min( max_x, old_x ) );
    auto vd = func( old_x );

    // early return ?
    TF old_error = vd.first;
    if ( abs( old_error ) < error_target )
        return old_x;

    //
    for( int num_iter = 0; ; ++num_iter ) {
        if ( num_iter >= max_iter )
            throw runtime_error( "newton 1D pb (max iter reached)" );
        
        const TF inc = old_error / vd.second;
        for( TF a = 1; ; ++num_iter, a /= 2 ) {
            if ( num_iter >= max_iter || a < 1e-6 ) {
                glot( linspace( min_x, max_x, 10000 ), [&]( TF x ) { return func(x).first; } );
                throw runtime_error( "newton 1D pb (bad direction)" );
            }

            TF new_x = max( min_x, min( max_x, old_x - a * inc ) );
            if ( new_x == old_x )
                return new_x;
            
            vd = func( new_x );

            if ( abs( vd.first ) < error_target )
                return new_x;
            
            if ( abs( vd.first ) < abs( old_error ) ) {
                old_error = vd.first;
                old_x = new_x;
                break;
            }
        }
    }
}

} // namespace usdot
