#pragma once

#include "utility/newton_1D.h"
#include "WeightInitializer.h"
#include <limits>

namespace usdot {
    
#define DTP template<class TF,class Density>
#define UTP WeightInitializer<TF,Density>

DTP UTP::WeightInitializer( Sys &sys ) : last_ag( nullptr ), sys( sys ) {
    coeff_ext_density = 1e-2;
}

DTP void UTP::run() {
    make_isolated_aggregates();
    for( nb_iterations = 0; last_ag_to_optimize; ++nb_iterations ) {
        optimize_the_new_aggregates();
        merge_touching_aggregates();
    }

    for( Ag *item = last_ag; item; item = item->prev )
        P( item->beg_u, item->end_u, item->len_n() );

    set_the_weights();
}

DTP void UTP::make_isolated_aggregates() {
    using namespace std;

    // push cells without taking care of the surrounding and start the first phase of agglomeration
    const TF u_tol = 100 * numeric_limits<TF>::epsilon();
    const TF x_tol = sys.x_tol();
    last_ag_to_optimize = nullptr;
    last_ag = nullptr;
    for( PI n = 0; n < sys.nb_sorted_diracs(); ++n ) {
        const TF x = sys.sorted_dirac_positions[ n ];
        const TF m = sys.sorted_dirac_masses[ n ];

        // find the cell position
        const TF b = min( TF( 1 - m ), max( TF( 0 ), newton_1D_unbounded<TF>( sys.density->cdf( x ) - m / 2, x_tol, [&]( TF b ) -> std::pair<TF,TF> {
            const TF err = sys.density->inv_cdf( b, coeff_ext_density ) + sys.density->inv_cdf( b + m, coeff_ext_density ) - 2 * x;
            const TF der = sys.density->der_inv_cdf( b, coeff_ext_density ) + sys.density->der_inv_cdf( b + m, coeff_ext_density );
            return { err, der };
        } ) ) );
        
        // create or place in an agglomerate
        const TF e = b + m;
        if ( last_ag == nullptr || last_ag->end_u + u_tol <= b ) { // non touching cell ?
            Ag *item = pool.create<Ag>();
            item->total_mass = m;
            item->beg_n = n + 0;
            item->end_n = n + 1;
            item->beg_u = b;
            item->end_u = e;

            item->prev = last_ag;
            last_ag = item;
        } else {
            if ( last_ag_to_optimize != last_ag ) {
                last_ag->prev_opt = last_ag_to_optimize;
                last_ag_to_optimize = last_ag;
            }
            if ( last_ag->end_u < e )
                last_ag->end_u = e;
            last_ag->total_mass += m;
            last_ag->end_n = n + 1;
        }
    }
}

DTP void UTP::merge_touching_aggregates() {
    using namespace std;
   
    const TF eps_u = 100 * numeric_limits<TF>::epsilon();
    last_ag_to_optimize = nullptr;
    for( Ag *item = last_ag; item; item = item->prev ) {
        while ( Ag *prev = item->prev ) {
            const TF delta = prev->end_u - item->beg_u;
            if ( delta + eps_u <= 0 )
                break;

            if ( last_ag_to_optimize != item ) {
                item->prev_opt = last_ag_to_optimize;
                last_ag_to_optimize = item;
            }
            item->total_mass += prev->total_mass;
            item->beg_n = prev->beg_n;
            item->beg_u = prev->beg_u;
            item->prev = prev->prev;
        }
    }
}

DTP void UTP::set_the_weights() {
    using namespace std;

    const TF eps_u = 100 * numeric_limits<TF>::epsilon();
    for( Ag *item = last_ag; item; item = item->prev ) {
        TF d0 = sys.sorted_dirac_positions[ item->beg_n ];
        TF u0 = item->beg_u;
        TF x0 = sys.density->inv_cdf( u0 );
        TF w0 = pow( d0 - x0, 2 );

        if ( item->beg_n + 1 == item->end_n ) { // cell with only 1 dirac
            if ( item->beg_u <= eps_u ) {
                if ( item->end_u >= 1 - eps_u ) {
                    const TF mi = min( sys.density->min_x(), sys.sorted_dirac_positions.front() );
                    const TF ma = max( sys.density->max_x(), sys.sorted_dirac_positions.back() );
                    w0 = pow( ma - mi, 2 ); // a large enough w0
                } else {
                    const PI en = item->end_n - 1;
                    const TF d1 = sys.sorted_dirac_positions[ en ];
                    const TF u1 = item->end_u;
                    const TF x1 = sys.density->inv_cdf( u1 );
                    w0 = pow( d1 - x1, 2 );
                }
            }
            sys.sorted_dirac_weights[ item->beg_n ] = w0;
            // P( w0, sys.density->integral( d0 - sqrt( w0 ), d0 + sqrt( w0 ) ) );
        } else { // several cells in item
            if ( item->beg_u <= eps_u ) {
                if ( item->end_u >= 1 - eps_u ) {
                    const TF mi = min( sys.density->min_x(), sys.sorted_dirac_positions.front() );
                    const TF ma = max( sys.density->max_x(), sys.sorted_dirac_positions.back() );
                    w0 = pow( ma - mi, 2 ); // a large enough w0
                    for( PI n = item->beg_n; ; ++n ) {
                        sys.sorted_dirac_weights[ n ] = w0;
                        if ( n + 1 == item->end_n )
                            break;
                                        
                        const TF d1 = sys.sorted_dirac_positions[ n + 1 ];
                        const TF u1 = u0 + sys.sorted_dirac_masses[ n ];
                        const TF x1 = sys.density->inv_cdf( u1 );
                        const TF w1 = w0 + ( d1 - d0 ) * ( d0 + d1 - 2 * x1 );
    
                        d0 = d1;
                        u0 = u1;
                        x0 = x1;
                        w0 = w1;
                    }
                } else {
                    PI n = item->end_n - 1;

                    TF d1 = sys.sorted_dirac_positions[ n ];
                    TF u1 = item->end_u;
                    TF x1 = sys.density->inv_cdf( u1 );
                    TF w1 = pow( d1 - x1, 2 );

                    for( ; ; --n ) {
                        sys.sorted_dirac_weights[ n ] = w1;
                        if ( n == item->beg_n )
                            break;
                                        
                        const TF d0 = sys.sorted_dirac_positions[ n - 1 ];
                        const TF u0 = u1 - sys.sorted_dirac_masses[ n ];
                        const TF x0 = sys.density->inv_cdf( u0 );
                        w0 = w1 - ( d1 - d0 ) * ( d0 + d1 - 2 * x0 );
    
                        d1 = d0;
                        u1 = u0;
                        x1 = x0;
                        w1 = w0;
                    }
                }
                // P( w0, sys.density->integral( d0 - sqrt( w0 ), x0 ) );
            } else {
                for( PI n = item->beg_n; ; ++n ) {
                    sys.sorted_dirac_weights[ n ] = w0;
                    if ( n + 1 == item->end_n )
                        break;
                                    
                    const TF d1 = sys.sorted_dirac_positions[ n + 1 ];
                    const TF u1 = u0 + sys.sorted_dirac_masses[ n ];
                    const TF x1 = sys.density->inv_cdf( u1 );
                    const TF w1 = w0 + ( d1 - d0 ) * ( d0 + d1 - 2 * x1 );

                    d0 = d1;
                    u0 = u1;
                    x0 = x1;
                    w0 = w1;
                }
                // P( w0, sys.density->integral( x0, d0 + sqrt( w0 ) ) );
            }
        }
    }
}

DTP void UTP::optimize_the_new_aggregates() {
    using namespace std;
    // if ( int( item->beg_u ) == 2208 ) {
    //     // glot( linspace( min_x, max_x, 1000 ), [&]( TF x ) { return func(x).first; } );
    //     sys.plot();
    //     assert( 0 );
    // }

    TF x_tol = sys.x_tol();
    for( Ag *item = last_ag_to_optimize; item; item = item->prev_opt ) {
        // new position
        const TF b = min( TF( 1 - item->total_mass ), max( TF( 0 ), newton_1D_unbounded<TF>( item->beg_u, x_tol, [&]( TF b ) {
            TF d0 = sys.density->der_inv_cdf( b, coeff_ext_density );
            TF x0 = sys.density->inv_cdf( b, coeff_ext_density );
            TF err = 0, der = 0;
            for( PI n = item->beg_n; n < item->end_n; ++n ) {
                b += sys.sorted_dirac_masses[ n ];

                const TF d1 = sys.density->der_inv_cdf( b, coeff_ext_density );
                const TF x1 = sys.density->inv_cdf( b, coeff_ext_density );
                const TF xd = sys.sorted_dirac_positions[ n ];

                err += pow( x1 - xd, 2 ) - pow( x0 - xd, 2 ); ///< TODO: x error (not x^2)
                der += 2 * ( d1 * ( x1 - xd ) - d0 * ( x0 - xd ) );

                d0 = d1;
                x0 = x1;
            }

            return std::pair<TF,TF>{ err, der };
        } ) ) );
        
        item->end_u = b + item->total_mass;
        item->beg_u = b;
    }
}

#undef DTP
#undef UTP

} // namespace usdot
