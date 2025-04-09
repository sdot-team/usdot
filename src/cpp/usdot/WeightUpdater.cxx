#pragma once

#include "utility/dichotomy.h"
#include "WeightUpdater.h"
#include <algorithm>
#include <stdexcept>
#include <limits>

namespace usdot {

#define DTP template<class TF,class Density>
#define UTP WeightUpdater<TF,Density>

DTP UTP::WeightUpdater( Sys &sys ) : sys( sys ) {
}

DTP void UTP::get_system( VC &connected_cells, VP &opt_weights, TF &max_a, TF &cell_error, int &has_bad_cell, const VF &sorted_dirac_weights ) const {
    using namespace std;

    //
    opt_weights.resize( sys.nb_sorted_diracs() );
    connected_cells.clear();
    has_bad_cell = 0;
    cell_error = 0;
    max_a = 1;
    
    TF i0 = numeric_limits<TF>::lowest(); // prev cut position
    TF d0 = sys.sorted_dirac_positions[ 0 ]; // prev dirac position
    TF w0 = sys.sorted_dirac_weights[ 0 ]; // prev weight
    TF dp = 0; // prev prev dirac position
    TF wp = 0; // prev prev weight

    for( PI n = 0; n < sys.nb_sorted_diracs(); ++n ) {
        // negative weight ?
        if ( w0 < 0 ) {
            has_bad_cell = 1;
            return;
        }

        // next intersection
        const TF d1 = n + 1 < sys.nb_sorted_diracs() ? sys.sorted_dirac_positions[ n + 1 ] : numeric_limits<TF>::max();
        const TF w1 = n + 1 < sys.nb_sorted_diracs() ? sys.sorted_dirac_weights[ n + 1 ] : 0;
        const TF i1 = ( d0 + d1 + ( w0 - w1 ) / ( d1 - d0 ) ) / 2;
        const TF r0 = sqrt( w0 );
        const TF b0 = d0 - r0;
        const TF b1 = d0 + r0;

        // negative cell ?
        if ( i1 < i0 ) {
            has_bad_cell = 2;
            return;
        }

        // void cell
        if ( i1 < b0 || i0 > b1 ) {
            has_bad_cell = 3;
            return;
        }

        if ( b0 > i0 ) { // ball cut on the left
            if ( b1 < i1 ) { // BB
                // optimal_r = a * ( target_area - integral( b0, b1 ) ) / ( value( b0 ) + value( b1 ) ) + rd
                const TF sva = sys.density->value( b0 ) + sys.density->value( b1 );
                const TF mde = sys.sorted_dirac_masses[ n ] / sys.density->width();
                const TF itg = sys.density->integral( b0, b1 );
                const TF err = sys.sorted_dirac_masses[ n ] - itg;
                connected_cells.push_back( { n, n } );
                cell_error += pow( err, 2 );

                opt_weights[ n ] = { 
                    .c0 = r0,
                    .c1 = 0,
                    .c2 = 0,
                    .ca = err / max( mde, sva )
                };
            } else { // BI
                // opt_weights[ n ] is equal to r^2
                connected_cells.push_back( { n, n } );
                opt_weights[ n ] = { 
                    .c0 = 0,
                    .c1 = 0,
                    .c2 = 1,
                    .ca = 0
                };

                const TF err = sys.density->integral( b0, i1 ) - sys.sorted_dirac_masses[ n ];
                cell_error += pow( err, 2 );

                // mass of current cell will precribe opt_weights[ n + 1 ]
                //   ( 1 - a ) * integral( i0, i1 ) + a * target_area = integral( i0, i1 ) + v0 * ( r - r0 ) + v1 * 0.5 * ( r^2 - w0 - ow + w1 ) / ( d1 - d0 )
                //   ( 1 - a ) * integral( i0, i1 ) + a * target_area = integral( i0, i1 ) + v0 * ( r - r0 ) + ( r^2 - w0 - ow + w1 ) / C
                //   a * ( target_area - integral( i0, i1 ) ) - v0 * ( r - r0 ) = ( r^2 - w0 - ow + w1 ) / C
                //   ( a * ( target_area - integral( i0, i1 ) ) - v0 * ( r - r0 ) ) * C = r^2 - w0 - ow + w1
                //   ow = r^2 + w1 - w0 + ( a * ( integral( i0, i1 ) - target_area ) + v0 * ( r - r0 ) ) * C
                const TF dv0 = sys.density->value( b0 );
                const TF dv1 = sys.density->value( i1 );
                // if ( dv0 == 0 || dv1 == 0 ) {
                //     has_bad_cell = 5;
                //     return;
                // }

                const TF mde = sys.sorted_dirac_masses[ n ] / sys.density->width();
                const TF v0 = max( mde, dv0 );
                const TF v1 = max( mde, dv1 );

                const TF D = 2 * ( d1 - d0 ) / v1;
                opt_weights[ n + 1 ] = {
                    .c0 = w1 - w0 - v0 * r0 * D,
                    .c1 = v0 * D,
                    .c2 = 1,
                    .ca = err * D
                };
            }
        } else { // interface on the left
            if ( b1 < i1 ) { // IB
                // const TF v1 = sys.density->value( b1 );
                // if ( v1 == 0 ) {
                //     has_bad_cell = 6;
                //     return;
                // }

                const TF err = sys.density->integral( i0, b1 ) - sys.sorted_dirac_masses[ n ];
                cell_error += pow( err, 2 );

                connected_cells.back()[ 1 ] = n;
            } else { // II
                const TF err = sys.density->integral( i0, i1 ) - sys.sorted_dirac_masses[ n ];
                cell_error += pow( err, 2 );

                // ( 1 - a ) * integral( i0, i1 ) + a * target_area = integral( i0, i1 ) + ( ow0 - w0 - owp + wp ) * C
                //                                                                       + ( ow0 - w0 - ow1 + w1 ) * D
                // a * ( target_area - integral( i0, i1 ) ) = ( ow0 - w0 - owp + wp ) * C
                //                                          + ( ow0 - w0 - ow1 + w1 ) * D
                // ow1 * D = a * ( integral( i0, i1 ) - target_area ) + ( ow0 - w0 - owp + wp ) ) * C
                //                                                    + ( ow0 - w0 + w1 ) * D
                const TF dv0 = sys.density->value( i0 );
                const TF dv1 = sys.density->value( i1 );
                // if ( dv1 == 0 ) {
                //     has_bad_cell = 7;
                //     return;
                // }

                const TF mde = sys.sorted_dirac_masses[ n ] / sys.density->width();
                const TF v0 = max( mde, dv0 );
                const TF v1 = max( mde, dv1 );

                const Poly &owp = opt_weights[ n - 1 ];
                const Poly &ow0 = opt_weights[ n - 0 ];
                const TF C = 0.5 * v0 / ( d0 - dp );
                const TF D = 0.5 * v1 / ( d1 - d0 );
                const TF E = C / D;
                opt_weights[ n + 1 ] = {
                    .c0 = ow0.c0 + E * ( ow0.c0 - owp.c0 - w0 + wp ) - w0 + w1,
                    .c1 = ow0.c1 + E * ( ow0.c1 - owp.c1 ),
                    .c2 = ow0.c2 + E * ( ow0.c2 - owp.c2 ),
                    .ca = ow0.ca + E * ( ow0.ca - owp.ca ) + err / D
                };
            }
        }

        //
        i0 = i1;
        dp = d0;
        d0 = d1;
        wp = w0;
        w0 = w1;
    }
}

DTP typename UTP::VF UTP::ce_errors( const VP &polys, PI nb, PI ne, TF a, TF r ) const {
    VF res;
    res << bi_error( polys, nb, a, r );
    for( PI n = nb + 1; n < ne; ++n )
        res << ii_error( polys, n, a, r );    
    return res;
}

DTP TF UTP::ii_error( const VP &polys, PI n, TF a, TF r ) const {
    using namespace std;

    const TF dp = sys.sorted_dirac_positions[ n - 1 ];
    const TF d0 = sys.sorted_dirac_positions[ n + 0 ];
    const TF d1 = sys.sorted_dirac_positions[ n + 1 ];

    const TF old_wp = sys.sorted_dirac_weights[ n - 1 ];
    const TF old_w0 = sys.sorted_dirac_weights[ n + 0 ];
    const TF old_w1 = sys.sorted_dirac_weights[ n + 1 ];
    const TF old_i0 = ( dp + d0 + ( old_wp - old_w0 ) / ( d0 - dp ) ) / 2;
    const TF old_i1 = ( d0 + d1 + ( old_w0 - old_w1 ) / ( d1 - d0 ) ) / 2;

    const TF t_mass = ( 1 - a ) * sys.density->integral( old_i0, old_i1 ) + a * sys.sorted_dirac_masses[ n ];

    const TF new_wp = polys[ n - 1 ].c0 + polys[ n - 1 ].c1 * r + polys[ n - 1 ].c2 * pow( r, 2 ) + polys[ n - 1 ].ca * a;
    const TF new_w0 = polys[ n + 0 ].c0 + polys[ n + 0 ].c1 * r + polys[ n + 0 ].c2 * pow( r, 2 ) + polys[ n + 0 ].ca * a;
    const TF new_w1 = polys[ n + 1 ].c0 + polys[ n + 1 ].c1 * r + polys[ n + 1 ].c2 * pow( r, 2 ) + polys[ n + 1 ].ca * a;
    const TF new_i0 = ( dp + d0 + ( new_wp - new_w0 ) / ( d0 - dp ) ) / 2;
    const TF new_i1 = ( d0 + d1 + ( new_w0 - new_w1 ) / ( d1 - d0 ) ) / 2;

    const TF n_mass = sys.density->integral( old_i0, old_i1 ) + sys.density->value( old_i0 ) * ( old_i0 - new_i0 ) + sys.density->value( old_i1 ) * ( new_i1 - old_i1 );

    return t_mass - n_mass;
}

DTP TF UTP::bi_error( const VP &polys, PI n, TF a, TF r ) const {
    using namespace std;

    const TF d0 = sys.sorted_dirac_positions[ n + 0 ];
    const TF d1 = sys.sorted_dirac_positions[ n + 1 ];

    const TF old_w0 = sys.sorted_dirac_weights[ n + 0 ];
    const TF old_w1 = sys.sorted_dirac_weights[ n + 1 ];
    const TF old_r0 = sqrt( old_w0 );
    const TF old_b0 = d0 - old_r0;
    const TF old_i1 = ( d0 + d1 + ( old_w0 - old_w1 ) / ( d1 - d0 ) ) / 2;

    const TF t_mass = ( 1 - a ) * sys.density->integral( old_b0, old_i1 ) + a * sys.sorted_dirac_masses[ n ];

    const TF new_w0 = polys[ n + 0 ].c0 + polys[ n + 0 ].c1 * r + polys[ n + 0 ].c2 * pow( r, 2 ) + polys[ n + 0 ].ca * a;
    const TF new_w1 = polys[ n + 1 ].c0 + polys[ n + 1 ].c1 * r + polys[ n + 1 ].c2 * pow( r, 2 ) + polys[ n + 1 ].ca * a;
    const TF new_r0 = sqrt( new_w0 );
    const TF new_b0 = d0 - new_r0;
    const TF new_i1 = ( d0 + d1 + ( new_w0 - new_w1 ) / ( d1 - d0 ) ) / 2;

    const TF n_mass = sys.density->integral( old_b0, old_i1 ) + sys.density->value( old_b0 ) * ( old_b0 - new_b0 ) + sys.density->value( old_i1 ) * ( new_i1 - old_i1 );

    return t_mass - n_mass;
}

DTP TF UTP::best_r_for_ib( const VP &polys, PI n, TF a ) const {
    using namespace std;

    const TF dp = sys.sorted_dirac_positions[ n - 1 ];
    const TF d0 = sys.sorted_dirac_positions[ n + 0 ];

    const TF old_wp = sys.sorted_dirac_weights[ n - 1 ];
    const TF old_w0 = sys.sorted_dirac_weights[ n + 0 ];
    const TF old_i0 = ( dp + d0 + ( old_wp - old_w0 ) / ( d0 - dp ) ) / 2;
    const TF old_rd = sqrt( old_w0 );
    const TF old_b1 = d0 + old_rd;

    const TF old_in = sys.density->integral( old_i0, old_b1 );
    const TF t_mass = ( 1 - a ) * old_in + a * sys.sorted_dirac_masses[ n ];

    const TF v0 = sys.density->value( old_i0 );
    const TF v1 = sys.density->value( old_b1 );

    const TF pp0 = polys[ n - 1 ].c0, pp1 = polys[ n - 1 ].c1, pp2 = polys[ n - 1 ].c2, ppa = polys[ n - 1 ].ca;
    const TF p00 = polys[ n + 0 ].c0, p01 = polys[ n + 0 ].c1, p02 = polys[ n + 0 ].c2, p0a = polys[ n + 0 ].ca;

    const TF ter = sys.target_max_mass_error * sys.sorted_dirac_masses[ n ] * 1e-5;

    auto err = [&]( const TF r ) {
        const TF new_wp = pp0 + pp1 * r + pp2 * pow( r, 2 ) + ppa * a;
        const TF new_w0 = p00 + p01 * r + p02 * pow( r, 2 ) + p0a * a;
        const TF new_i0 = ( dp + d0 + ( new_wp - new_w0 ) / ( d0 - dp ) ) / 2;
        const TF new_rd = sqrt( max( TF( 0 ), new_w0 ) );
        const TF new_b1 = d0 + new_rd;

        const TF n_mass = old_in + v0 * ( old_i0 - new_i0 ) + v1 * ( new_b1 - old_b1 );
        
        return n_mass - t_mass;
    };

    if ( p02 ) {
        const TF b = pow( p01, 2 ) - 4 * p02 * ( p00 + p0a * a );
        if ( b >= 0 ) {
            const TF xa = ( - p01 - sqrt( b ) ) / ( 2 * p02 );
            const TF xb = ( - p01 + sqrt( b ) ) / ( 2 * p02 );
            const TF x0 = max( TF( 0 ), min( xa, xb ) );
            const TF x1 = max( TF( 0 ), max( xa, xb ) );
            if ( p02 < 0 ) { // solution must be in [ x0, x1 ]
                if ( x0 > 0 && abs( err( x0 ) ) <= ter )
                    return x0;
                if ( x1 > 0 && abs( err( x1 ) ) <= ter )
                    return x1;
                if ( err( max( TF( 0 ), x0 ) ) * err( x1 ) > 0 )
                    return -2;
                return dichotomy( err, ter, max( TF( 0 ), x0 ), x1 );
            }
            
            // solution may be in [ 0, x0 ]
            if ( x0 > 0 ) {
                if ( abs( err( 0 ) ) <= ter )
                    return 0;
                if ( abs( err( x0 ) ) <= ter )
                    return x0;
                if ( err( 0 ) * err( x0 ) <= 0 )
                    return dichotomy( err, ter, TF( 0 ), x0 );
            }

            // else, solution must be in [ x1, inf ]
            if ( x1 > 0 && abs( err( x1 ) ) <= ter )
                return x1;
            if ( err( x1 ) > 0 )
                return -3;
 
            TF max_r = x1 ? 2 * x1 : 1;
            while ( err( max_r ) < 0 ) {
                if ( max_r > sys.density->width() )
                    return -4;
                max_r *= 2;
            }
            // P( __LINE__, err( x1 ), err( max_r ), max_r, ter );
            return dichotomy( err, ter, x1, max_r );
        }

        //
        if ( p00 < 0 )
            return -5;

        if ( abs( err( 0 ) ) <= ter )
            return 0;
        if ( err( 0 ) > 0 )
            return -6;

        TF max_r = old_rd ? 2 * old_rd : 1;
        while ( err( max_r ) < 0 ) {
            if ( max_r > sys.density->width() )
                return -7;
            max_r *= 2;
        }
        // P( __LINE__, max_r );
        return dichotomy( err, ter, TF( 0 ), max_r );
    }

    throw runtime_error( "TODO" );
}

DTP int UTP::get_weights_for( VF &new_dirac_weights, const VC &connected_cells, const VP &polys, TF a ) {
    new_dirac_weights.resize( sys.nb_sorted_diracs() );
    for( std::array<PI,2> inds : connected_cells ) {
        if ( inds[ 0 ] == inds[ 1 ] ) { 
            new_dirac_weights[ inds[ 0 ] ] = pow( polys[ inds[ 0 ] ].c0 + a * polys[ inds[ 0 ] ].ca, 2 );
        } else {
            // solve
            // P( ce_errors( polys, inds[ 0 ], inds[ 1 ], 0.5, 0.5 ) );
            TF r = best_r_for_ib( polys, inds[ 1 ], a );
            if ( r < 0 )
                return - r;

            for( PI n = inds[ 0 ]; n <= inds[ 1 ]; ++n )
                new_dirac_weights[ n ] = polys[ n ].c0 + polys[ n ].c1 * r + polys[ n ].c2 * pow( r, 2 ) + polys[ n ].ca * a;
        }
    }

    return 0;
}

DTP void UTP::run() {
    if ( sys.nb_sorted_diracs() == 0 )
        return;

    // system content
    VC new_connected_cells;
    VC connected_cells;
    VP new_opt_weights;
    TF new_cell_error;
    int has_bad_cell;
    VP opt_weights;
    TF cell_error;
    TF max_a;

    // history
    // Vec<Vec<Vec<TF,2>>> bnds;

    // get the first system
    get_system( connected_cells, opt_weights, max_a, cell_error, has_bad_cell, sys.sorted_dirac_weights );
    if ( has_bad_cell )
        throw std::runtime_error( "bad initialization" );
    if ( sys.verbosity >= 3 && sys.stream )
        *sys.stream << "cell_error:" << cell_error << " (first iteration)" << std::endl;
    nb_iterations = 0;
    if ( cell_error < sys.target_max_mass_error / sys.nb_sorted_diracs() )
        return;

    // iterate
    VF new_sorted_dirac_weights;
    for( ; ; ++nb_iterations ) {
        if ( nb_iterations == 1000 )
            throw std::runtime_error( "too much iterations" );
        // bnds << normalized_cell_boundaries();

        for( TF a = max_a; ; a *= 0.8 ) {
            if ( a < 1e-20 ) {
                // plot_bnds_evolution( bnds );
                throw std::runtime_error( "bad direction" );
            }

            int esw = get_weights_for( new_sorted_dirac_weights, connected_cells, opt_weights, a );
            if ( esw == 0 ) {
                get_system( new_connected_cells, new_opt_weights, max_a, new_cell_error, has_bad_cell, new_sorted_dirac_weights );
                if ( has_bad_cell == 0 && new_cell_error < cell_error ) {
                    std::swap( new_sorted_dirac_weights, sys.sorted_dirac_weights );
                    std::swap( new_connected_cells, connected_cells );
                    std::swap( new_opt_weights, opt_weights );
                    std::swap( new_cell_error, cell_error );

                    if ( sys.verbosity > 0 && sys.stream )
                        *sys.stream << "cell_error:" << cell_error << " a:" << a << std::endl;

                    if ( sqrt( cell_error ) < sys.target_max_mass_error / sys.nb_sorted_diracs() )
                        return;
                    break;
                }

                // P( has_bad_cell, new_cell_error, cell_error, a );
            }
            // P( esw );
        }
    }
}

#undef DTP
#undef UTP

} // namespace usdot
