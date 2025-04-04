#pragma once

#include "WeightInitializer.h"
#include "partial1D/dichotomy.h"
#include <limits>

namespace usdot {
    
#define DTP template<class TF,class Density>
#define UTP WeightInitializer<TF,Density>

DTP UTP::WeightInitializer( Sys &sys ) : sys( sys ), last_ag( nullptr ) {
}

DTP void UTP::run() {
    using namespace std;

    for( TI n = 0; n < sys.nb_diracs(); ++n ) {
        const TF coord = sys.sorted_dirac_positions[ n ];
        const TF mass = sys.sorted_dirac_masses[ n ];
        const TF r = optimal_radius( coord, mass );
        append_cell( n, coord - r, coord + r );
    }

    // // get weights
    // for( Agglomeration *item = last; item; item = item->prev ) {
    //     if ( item->len_n() == 0 )
    //         continue;

    //     // data from the last cell
    //     const PI n = item->end_n - 1;
    //     TF c = seed_coords[ n ];

    //     TF end_y = item->end_y(), beg_y = end_y - mass_ratios[ n ];
    //     TF end_x = get_x( end_y ), beg_x = get_x( beg_y );

    //     // compute the weight from the required radius
    //     TF w = 0;
    //     seed_weights[ item->end_n - 1 ] = w;

    //     //  
    //     TF w_to_add = max( pow( end_x - c, 2 ), pow( c - beg_x, 2 ) );
    //     for( PI n = item->end_n; --n > item->beg_n; ) {
    //         // update bounds
    //         end_y = std::exchange( beg_y, max( 0, beg_y - mass_ratios[ n ] ) );
    //         end_x = std::exchange( beg_x, get_x( beg_y ) );

    //         // update w
    //         TF d = std::exchange( c, seed_coords[ n - 1 ] );
    //         w -= ( d + c - 2 * end_x ) * ( d - c );
            
    //         const TF min_w = max( pow( c - beg_x, 2 ), pow( end_x - c, 2 ) );
    //         if ( min_w > w )
    //             w_to_add = max( w_to_add, min_w - w );

    //         // store the result
    //         seed_weights[ n - 1 ] = w;
    //     }

    //     if ( w_to_add )
    //         for( PI n = item->beg_n; n < item->end_n; ++n )
    //             seed_weights[ n ] += w_to_add;
    // }
}

DTP TF UTP::optimal_radius( TF pos, TF mass ) const {
    auto err = [&]( TF r ) {
        return sys.density->integral( pos - r, pos + r ) - mass;
    };
    return dichotomy_growing_from_zero( err, mass * 1e-2, TF( 1 ) );
}

DTP void UTP::append_cell( PI n, TF beg_x, TF end_x ) {
    // non touching cell ?
    if ( last_ag == nullptr || last_ag->end_x < beg_x ) {
        Ag *item = pool.create<Ag>();
        item->beg_x = beg_x;
        item->end_x = end_x;
        item->beg_n = n + 0;
        item->end_n = n + 1;

        item->prev = last_ag;
        last_ag = item;
        
        return;
    }

    // else
    const TF mass_to_distribute = sys.density->integral( beg_x, last_ag->end_x );
    if ( mass_to_distribute ) {
        const TF cn = TF( 1 ) / ( last_ag->len_n() + 1 );
        const TF cp = 1 - cn;
        auto err = [&]( const TF r ) {
            return sys.density->integral( last_ag->beg_x - r * cp, last_ag->beg_x ) + sys.density->integral( end_x, end_x + r * cn ) - mass_to_distribute;
        };
        const TF new_r = dichotomy_growing_from_zero( err, 1e-2 * sys.sorted_dirac_masses[ n ], last_ag->end_x - beg_x );
        last_ag->end_x = end_x + new_r * cn;
        last_ag->beg_x -= new_r * cp;
    }

    //  append to the previous block
    last_ag->end_n = n + 1;

    // while we have an Agglomerate to merge
    while ( true ) {
        Ag *prev = last_ag->prev;
        if ( prev == nullptr || prev->end_x < last_ag->beg_x )
            break;

        // append the last set
        const TF mass_to_distribute = sys.density->integral( last_ag->beg_x, prev->end_x );
        if ( mass_to_distribute ) {
            const TF cn = prev->len_n() / TF( last_ag->len_n() + prev->len_n() );
            const TF cp = 1 - cn;
            auto err = [&]( const TF r ) {
                return sys.density->integral( prev->beg_x - r * cp, prev->beg_x ) + sys.density->integral( last_ag->end_x, last_ag->end_x + r * cn ) - mass_to_distribute;
            };
            const TF new_r = dichotomy_growing_from_zero( err, 1e-2 * sys.sorted_dirac_masses[ n ], last_ag->end_x - beg_x );
            last_ag->end_x = end_x + new_r * cn;
            last_ag->beg_x -= new_r * cp;
        }
    
        //  append to the previous block
        prev->end_n = last_ag->end_n;
    }
}

#undef DTP
#undef UTP

} // namespace usdot
