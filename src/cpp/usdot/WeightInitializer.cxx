#pragma once

#include "utility/newton_1D.h"
#include "WeightInitializer.h"
#include <algorithm>
#include <stdexcept>
#include <limits>

namespace usdot {
    
#define DTP template<class TF,class Density>
#define UTP WeightInitializer<TF,Density>

DTP UTP::WeightInitializer( Sys &sys ) : last_ag( nullptr ), sys( sys ) {
    sys.density->get_inv_cdf( inv_cdf_values, mul_coeff, 100 );
    
    dirac_masses.resize( sys.nb_sorted_diracs() );
    for( PI n = 0; n < sys.nb_sorted_diracs(); ++n )
        dirac_masses[ n ] = mul_coeff * sys.sorted_dirac_masses[ n ];
}

DTP void UTP::run() {
    make_isolated_aggregates();

    // P( inv_cdf_values.size() - 1 );
    // for( Ag *item = last_ag; item; item = item->prev ) {
    //     P( item->beg_u, item->end_u );
    // }
    
    for( nb_iterations = 0; last_ag_to_optimize; ++nb_iterations ) {
        for( Ag *item = last_ag_to_optimize; item; item = item->prev_opt )
            optimize_aggregate( item, sys.x_tol() / 10 );
        merge_touching_aggregates();
    }

    set_the_weights();
}

DTP TF UTP::der_inv_cdf( TF u ) const {
    using namespace std;
    
    PI n = min( inv_cdf_values.size() - 2, PI( u ) );
    return inv_cdf_values[ n + 1 ] - inv_cdf_values[ n + 0 ];
}

DTP TF UTP::inv_cdf( TF u ) const {
    using namespace std;

    const PI m = inv_cdf_values.size() - 1;
    u = max( TF( 0 ), min( TF( m ), u ) );

    PI n = min( m - 1, PI( u ) );
    TF f = u - n;
    
    return inv_cdf_values[ n + 0 ] * ( 1 - f ) + inv_cdf_values[ n + 1 ] * f;
}

DTP TF UTP::cdf( TF x ) const {
    return mul_coeff * sys.density->primitive( x );
}

DTP void UTP::make_isolated_aggregates() {
    using namespace std;

    // push cells without taking care of the surrounding and start the first phase of agglomeration
    const TF eps_u = max_u() * 100 * numeric_limits<TF>::epsilon();
    last_ag_to_optimize = nullptr;
    const TF x_tol = sys.x_tol() / 10;
    last_ag = nullptr;
    for( PI n = 0; n < sys.nb_sorted_diracs(); ++n ) {
        const TF x = sys.sorted_dirac_positions[ n ];
        const TF m = dirac_masses[ n ];
        const TF c = cdf( x );

        // find the cell position
        const TF b = newton_1D<TF>( c, 0, inv_cdf_values.size() - 1 - m, x_tol, [&]( TF b ) {
            const TF e = inv_cdf( b ) + inv_cdf( b + m ) - 2 * x;
            const TF d = der_inv_cdf( b ) + der_inv_cdf( b + m );
            return std::pair<TF,TF>{ e, d };
        } );

        // auto err = [&]( const TF b ) {
        //     return inv_cdf( b ) + inv_cdf( b + m ) - 2 * x;
        // };        
        // const TF v = min( c + 2 * m, inv_cdf_values.size() - 1 ) - m;
        // const TF u = max( c - m, 0 );
        // const TF b = dichotomy( err, 1e-6 * m, u, v );

        // create or place in an agglomerate
        if ( last_ag == nullptr || last_ag->end_u + eps_u <= b ) { // non touching cell ?
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
            if ( last_ag->end_u < b + m )
                last_ag->end_u = b + m;
            last_ag->end_n = n + 1;
        }
    }
}

DTP void UTP::merge_touching_aggregates() {
    using namespace std;
   
    const TF eps_u = max_u() * 100 * numeric_limits<TF>::epsilon();
    last_ag_to_optimize = nullptr;
    for( Ag *item = last_ag; item; item = item->prev ) {
        while ( Ag *prev = item->prev ) {
            const TF delta = prev->end_u - item->beg_u;
            if ( delta + eps_u <= 0 )
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

DTP void UTP::set_the_weights() {
    using namespace std;

    const TF eps_u = max_u() * 100 * numeric_limits<TF>::epsilon();
    for( Ag *item = last_ag; item; item = item->prev ) {
        TF d0 = sys.sorted_dirac_positions[ item->beg_n ];
        TF u0 = item->beg_u;
        TF x0 = inv_cdf( u0 );
        TF w0 = pow( d0 - x0, 2 );

        if ( item->beg_n + 1 == item->end_n ) { // cell with only 1 dirac
            if ( item->beg_u <= eps_u ) {
                if ( item->end_u >= max_u() - eps_u ) {
                    w0 = pow( sys.density->width(), 2 ); // a large enough w0
                } else {
                    const PI en = item->end_n - 1;
                    const TF d1 = sys.sorted_dirac_positions[ en ];
                    const TF u1 = item->end_u;
                    const TF x1 = inv_cdf( u1 );
                    w0 = pow( d1 - x1, 2 );
                }
            }
            sys.sorted_dirac_weights[ item->beg_n ] = w0;
            // P( w0, sys.density->integral( d0 - sqrt( w0 ), d0 + sqrt( w0 ) ) );
        } else { // several cells in item
            if ( item->beg_u <= eps_u ) {
                if ( item->end_u >= max_u() - eps_u ) {
                    w0 = pow( sys.density->width(), 2 ); // a large enough w0
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
                } else {
                    PI n = item->end_n - 1;

                    TF d1 = sys.sorted_dirac_positions[ n ];
                    TF u1 = item->end_u;
                    TF x1 = inv_cdf( u1 );
                    TF w1 = pow( d1 - x1, 2 );

                    for( ; ; --n ) {
                        sys.sorted_dirac_weights[ n ] = w1;
                        if ( n == item->beg_n )
                            break;
                                        
                        const TF d0 = sys.sorted_dirac_positions[ n - 1 ];
                        const TF u0 = u1 - dirac_masses[ n ];
                        const TF x0 = inv_cdf( u0 );
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
                    const TF u1 = u0 + dirac_masses[ n ];
                    const TF x1 = inv_cdf( u1 );
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

DTP void UTP::optimize_aggregate( Ag *item, TF x_tol ) {
    using namespace std;

    // mass of the aggregate
    TF m = 0;
    for( PI n = item->beg_n; n < item->end_n; ++n )
        m += dirac_masses[ n ];

    // Newton
    const TF b = newton_1D<TF>( item->beg_u, 0, inv_cdf_values.size() - 1 - m, x_tol, [&]( TF b ) {
        TF e = 0, d = 0;
        TF x0 = inv_cdf( b );
        TF d0 = der_inv_cdf( b );
        for( PI n = item->beg_n; n < item->end_n; ++n ) {
            b += dirac_masses[ n ];

            const TF x1 = inv_cdf( b );
            const TF d1 = der_inv_cdf( b );
            const TF xd = sys.sorted_dirac_positions[ n ];

            e += pow( x1 - xd, 2 ) - pow( x0 - xd, 2 );
            d += 2 * ( d1 * ( x1 - xd ) - d0 * ( x0 - xd ) );

            d0 = d1;
            x0 = x1;
        }

        // const TF e = inv_cdf( b ) + inv_cdf( b + m ) - 2 * x;
        // const TF d = der_inv_cdf( b ) + der_inv_cdf( b + m );
        return std::pair<TF,TF>{ e, d };
    } );
      
    item->end_u = b + m;
    item->beg_u = b;

    // for( Ag *item = last_ag; item; item = item->prev ) {
    //     if ( item->len_n() == 1 )
    //         continue;

    //     const TF di = sys.sorted_dirac_positions[ item->beg_n ];
    //     TF u0 = item->beg_u;
    //     TF u1 = u0 + dirac_masses[ item->beg_n ];
    //     TF x0 = inv_cdf( u0 );
    //     TF x1 = inv_cdf( u1 );
    //     TF f0 = pow( x1 - di, 2 ) - pow( x0 - di, 2 );
    //     for( PI n = item->beg_n + 1; n < item->end_n; ++n ) {
    //         u0 = u1;
    //         x0 = x1;
            
    //         u1 = u0 + dirac_masses[ n ];
    //         x1 = inv_cdf( u1 );

    //         const TF di = sys.sorted_dirac_positions[ n ];
    //         const TF f1 = pow( x1 - di, 2 ) - pow( x0 - di, 2 );

    //         P( f1 - f0 );
    //     }
    // }


    // TF m = 0;
    // for( PI n = item->beg_n; n < item->end_n; ++n )
    //     m += dirac_masses[ n ];

    // auto err = [&]( TF b ) {
    //     TF res = 0;
    //     TF xb = inv_cdf( b );
    //     for( PI n = item->beg_n; n < item->end_n; ++n ) {
    //         b += dirac_masses[ n ];

    //         const TF xd = sys.sorted_dirac_positions[ n ];
    //         const TF xe = inv_cdf( b );

    //         res += pow( xe - xd, 2 ) - pow( xb - xd, 2 );

    //         xb = xe;
    //     }

    //     return res;
    // };

    // const TF c = 1e-6 * pow( sys.sorted_dirac_masses[ item->beg_n ], 2 );
    // const TF v = min( item->beg_u + m, inv_cdf_values.size() - 1 ) - m;
    // const TF u = max( item->beg_u - 2 * item->test_width, 0 );
    // glot( Vec<TF>::linspace( u, v, 1000 ), err );
    // TODO;

    // const TF b = dichotomy( err, c, u, v );
    // item->end_u = b + m;
    // item->beg_u = b;
}

#undef DTP
#undef UTP

} // namespace usdot
