#pragma once

#include <stdexcept>
#include <cmath>

namespace usdot {

template<class TF,class Func>
TF newton_1D( TF start_x, TF min_x, TF max_x, TF error_target, const Func &func, int max_iter = 50000 ) {
    using namespace std;

    // first call
    start_x = max( min_x, min( max_x, start_x ) );
    auto vd = func( start_x );

    // early return ?
    TF old_error = vd.first;
    if ( abs( old_error ) < error_target )
        return start_x;

    //
    for( int num_iter = 0; ; ++num_iter ) {
        if ( num_iter >= max_iter )
            throw runtime_error( "newton pb" );
        
        const TF inc = old_error / vd.second;
        for( TF a = 1; ; ++num_iter, a /= 2 ) {
            if ( num_iter >= max_iter || a < 1e-12 )
                throw runtime_error( "newton pb" );

            TF new_x = max( min_x, min( max_x, start_x - a * inc ) );
            if ( new_x == start_x )
                return new_x;
            
            vd = func( new_x );

            if ( abs( vd.first ) < error_target )
                return new_x;
            
            if ( abs( vd.first ) < abs( old_error ) ) {
                old_error = vd.first;
                start_x = new_x;
                break;
            }
        }
    }
}

} // namespace usdot
