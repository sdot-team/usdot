#pragma once

#include "partial1D/dichotomy.h"
#include "WeightInitializer.h"

namespace usdot {
    
#define DTP template<class TF,class Density>
#define UTP WeightInitializer<TF,Density>

DTP UTP::WeightInitializer( Sys &sys ) : sys( sys ), last_ag( nullptr ) {
    sys.density->get_inv_cdf( inv_cdf_values, mul_coeff, 500 );    
    dirac_masses = mul_coeff * sys.sorted_dirac_masses;
}

DTP TF UTP::inv_cdf( TF u ) const {
    const TI m = inv_cdf_values.size() - 1;
    u = max( 0, min( m, u ) );

    PI n = min( m - 1, PI( u ) );
    TF f = u - n;
    
    return inv_cdf_values[ n + 0 ] * ( 1 - f ) + inv_cdf_values[ n + 1 ] * f;
}

DTP TF UTP::cdf( TF x ) const {
    return mul_coeff * sys.density->primitive( x );
}

DTP void UTP::run() {
    using namespace std;

    // push cells without taking care of the surrounding and start the first phase of agglomeration
    Ag *last_ag_to_optimize = nullptr;
    last_ag = nullptr;
    for( TI n = 0; n < sys.nb_diracs(); ++n ) {
        const TF x = sys.sorted_dirac_positions[ n ];
        const TF m = dirac_masses[ n ];
        const TF c = cdf( x );

        // find the cell position
        auto err = [&]( const TF b ) {
            return inv_cdf( b ) + inv_cdf( b + m ) - 2 * x;
        };        
        const TF b = dichotomy( err, 1e-4 * m, c - m, c + m );

        // create or place in an agglomerate
        if ( last_ag == nullptr || last_ag->end_u <= b ) { // non touching cell ?
            Ag *item = pool.create<Ag>();
            item->beg_n = n + 0;
            item->end_n = n + 1;
            item->beg_u = b;
            item->end_u = b + m;

            item->prev = last_ag;
            last_ag = item;
        } else {
            if ( last_ag_to_optimize != last_ag ) {
                last_ag->prev_opt = last_ag_to_optimize;
                last_ag_to_optimize = last_ag;
                last_ag->test_width = 0;
            }
            last_ag->test_width += last_ag->end_u - b;
            last_ag->end_n = n + 1;
        }
    }

    // optimize and merge the touching agglomerate
    while ( last_ag_to_optimize ) {
        // change the positions
        for( Ag *item = last_ag_to_optimize; item; item = item->prev_opt )
            optimize_agglomerate( item );

        // merge the touching agglomerates
        last_ag_to_optimize = nullptr;
        for( Ag *item = last_ag; item; item = item->prev ) {
            while ( Ag *prev = item->prev ) {
                const TF delta = prev->end_u - item->beg_u;
                if ( delta <= 0 )
                    break;

                if ( last_ag_to_optimize != item ) {
                    item->prev_opt = last_ag_to_optimize;
                    last_ag_to_optimize = item;
                    item->test_width = 0;
                }
                item->test_width += delta;
                item->beg_n = prev->beg_n;
                item->beg_u = prev->beg_u;
                item->prev = prev->prev;
            }
        }
    }

    // get the weights
    for( Ag *item = last_ag; item; item = item->prev ) {
        TF d0 = sys.sorted_dirac_positions[ item->beg_n ];
        TF u0 = item->beg_u;
        TF x0 = inv_cdf( u0 );
        TF w0 = pow( d0 - x0, 2 );

        if ( item->beg_n + 1 == item->end_n ) { // cell with only 1 dirac
            sys.sorted_dirac_weights[ item->beg_n ] = w0;
        } else { // several cells in item
            for( PI n = item->beg_n; ; ++n ) {
                sys.sorted_dirac_weights[ n ] = w0;
                if ( n + 1 == item->end_n )
                    break;
                                
                const TF d1 = sys.sorted_dirac_positions[ n + 1 ];
                const TF u1 = u0 + dirac_masses[ n ];
                const TF x1 = inv_cdf( u1 );
                const TF w1 = w0 + ( d1 - d0 ) * ( d0 + d1 - 2 * x1 );

                d0 = d1;
                u0 = u1;
                x0 = x1;
                w0 = w1;
            }
        }
    }
}

DTP void UTP::optimize_agglomerate( Ag *item ) {
    using namespace std;

    TF m = 0;
    for( PI n = item->beg_n; n < item->end_n; ++n )
        m += dirac_masses[ n ];

    auto err = [&]( TF b ) {
        TF res = 0;
        TF xb = inv_cdf( b );
        for( PI n = item->beg_n; n < item->end_n; ++n ) {
            b += dirac_masses[ n ];

            const TF xd = sys.sorted_dirac_positions[ n ];
            const TF xe = inv_cdf( b );

            res += pow( xe - xd, 2 ) - pow( xb - xd, 2 );

            xb = xe;
        }

        return res;
    };

    const TF c = 1e-6 * pow( sys.sorted_dirac_masses[ item->beg_n ], 2 );
    const TF b = dichotomy( err, c, item->beg_u - 2 * item->test_width, item->beg_u );
    item->end_u = b + m;
    item->beg_u = b;
}

#undef DTP
#undef UTP

} // namespace usdot
