#pragma once

// #include "linspace.h"
// #include "linspace.h"
#include "dichotomy.h"
#include <stdexcept>
#include <cmath>

namespace usdot {

template<class TF,class Func>
TF newton_1D_unbounded( TF start_x, TF error_target, const Func &func, int max_iter = 20 ) {
    using namespace std;

    // first call
    TF old_x = start_x;
    auto vd = func( old_x );

    // already in the error target ?
    TF old_error = vd.first;
    if ( abs( old_error ) < error_target )
        return old_x;

    // else, make iterations
    for( int num_iter = 0; ; ++num_iter ) {
        const TF inc = old_error / vd.second;
        for( TF a = 1; ; ++num_iter, a /= 2 ) {
            if ( num_iter >= max_iter || a < 1e-9 || vd.second == 0 )
                return dichotomy_unbounded( start_x, error_target, [&]( TF x ) { return func( x ).first; } );

            TF new_x = old_x - a * inc;
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
