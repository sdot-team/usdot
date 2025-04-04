#pragma once

#include <tl/support/operators/argmin.h>
#include <tl/support/ASSERT.h>
#include <tl/support/P.h>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <limits>

#include "ConvGridSolver.h"
#include "dichotomy.h"

namespace usdot {

    
#define DTP template<class TF>
#define UTP ConvGridSolver<TF>

DTP UTP::ConvGridSolver( ConvGridSolverInput<TF> &&input ) {
    using namespace std;

    // check inputs
    if ( input.dirac_positions.empty() )
        return;

    if ( input.density_values.size() <= 1 )
        throw runtime_error( "one needs at least 2 density value" );

    // values that can be copied
    current_filter_value = input.starting_filter_value;
    target_filter_value = input.target_filter_value;
    global_mass_ratio = input.global_mass_ratio;
    target_mass_error = input.target_mass_error;
    throw_if_error = input.throw_if_error;
    end_x_density = input.end_x_density;
    beg_x_density = input.beg_x_density;
    multithread = input.multithread;
    verbosity = input.verbosity;

    // values that can be moved
    original_density_values = std::move( input.density_values );

    // sorted_dirac_nums
    const PI nb_diracs = input.dirac_positions.size();
    sorted_dirac_nums.resize( nb_diracs );
    iota( sorted_dirac_nums.begin(), sorted_dirac_nums.end(), PI( 0 ) );
    sort( sorted_dirac_nums.begin(), sorted_dirac_nums.end(), [&]( PI a, PI b ) {
        return input.dirac_positions[ a ] < input.dirac_positions[ b ];
    } );

    // sorted_dirac_positions, sorted_dirac_weights, sorted_dirac_masses, sum_of_dirac_masses
    const TF mul_dirac = ( original_density_values.size() - 1 ) / ( end_x_density - beg_x_density );
    //const TF weight = 1 * pow( TF( original_density_values.size() - 1 ) / nb_diracs, 2 );
    const TF sep = input.min_dirac_separation * mul_dirac;
    sorted_dirac_positions.resize( nb_diracs );
    sorted_dirac_masses.resize( nb_diracs );
    sum_of_dirac_masses = 0;
    for( PI i = 0; i < nb_diracs; ++i ) {
        sorted_dirac_positions[ i ] = ( input.dirac_positions[ sorted_dirac_nums[ i ] ] - beg_x_density ) * mul_dirac;
        if ( i && sorted_dirac_positions[ i ] < sorted_dirac_positions[ i - 1 ] + sep ) 
            sorted_dirac_positions[ i ] = sorted_dirac_positions[ i - 1 ] + sep;
    
        const TF dirac_mass = sorted_dirac_nums[ i ] < input.dirac_masses.size() ? input.dirac_masses[ sorted_dirac_nums[ i ] ] : 1;
        sorted_dirac_masses[ i ] = dirac_mass;
        sum_of_dirac_masses += dirac_mass;
    }

    // filter
    initialize_filter_value( input.starting_filter_value );

    // weights
    initialize_weights();
}

DTP void UTP::initialize_weights() {
    // const TF base_weight = pow( 0.5 * global_mass_ratio * TF( original_density_values.size() ) / nb_diracs, 2 );
    // sorted_dirac_weights.resize( nb_diracs );
    // for( PI i = 0; i < nb_diracs; ++i ) {
    //     sorted_dirac_weights[ i ] = sorted_dirac_nums[ i ] < input.dirac_weights.size() ? 
    //         input.dirac_weights[ sorted_dirac_nums[ i ] ] * pow( mul_dirac, 2 ) : 
    //         base_weight;
    // }
    TODO;
}

DTP PI UTP::nb_diracs() const {
    return sorted_dirac_positions.size();
}

DTP void UTP::initialize_filter_value( TF filter_value ) {
    current_filter_value = filter_value;

    PI mul_x = 1;
    TF cut_ratio = 1e-6, target_mass = sum_of_dirac_masses / global_mass_ratio;
    density = DensityPtr{ new GridDensity<TF>( original_density_values, filter_value, mul_x, cut_ratio, target_mass ) };
}

DTP void UTP::go_to_filter_value( TF filter_value ) {
    // Vec old_weights = sorted_dirac_weights;

    initialize_filter_value( filter_value );
}

DTP void UTP::for_each_normalized_cell( auto &&func ) const {
    using namespace std;

    TF i0 = numeric_limits<TF>::lowest();
    TF d0 = sorted_dirac_positions[ 0 ];
    TF w0 = sorted_dirac_weights[ 0 ];
    TF od = numeric_limits<TF>::max();
    for( PI n = 1; n < nb_diracs(); ++n ) {
        const TF d1 = sorted_dirac_positions[ n ];
        const TF w1 = sorted_dirac_weights[ n ];

        const TF i1 = ( d0 + d1 + ( w0 - w1 ) / ( d1 - d0 ) ) / 2;
        func( d0, w0, n - 1, od, d1 - d0, i0, i1 );

        od = d1 - d0;
        i0 = i1;
        d0 = d1;
        w0 = w1;
    }

    func( d0, w0, nb_diracs() - 1, od, numeric_limits<TF>::max(), i0, numeric_limits<TF>::max() );
}

DTP void UTP::for_each_normalized_cell_mass( auto &&func ) const {
    using namespace std;
    for_each_normalized_cell( [&]( TF dirac_position, TF dirac_weight, PI num_dirac, TF, TF, TF c0, TF c1 ) {
        const TF rd = sqrt( dirac_weight );
        const TF b0 = dirac_position - rd;
        const TF b1 = dirac_position + rd;

        if ( c0 > c1 ) {
            const TF a0 = min( b1, c0 );
            const TF a1 = max( b0, c1 );
            return func( num_dirac, density->integral( a0, a1 ), a1 > a0 );
        }

        const TF a0 = max( b0, c0 );
        const TF a1 = min( b1, c1 );
        return func( num_dirac, density->integral( a0, a1 ), a0 > a1 );
    } );
}

DTP Vec<TF> UTP::normalized_cell_masses() const {
    TV res( FromSize(), nb_diracs() );
    for_each_normalized_cell_mass( [&]( PI num_dirac, TF mass, bool ) {
        res[ num_dirac ] = mass;
    } );
    return res;
}

DTP std::vector<std::array<TF,2>> UTP::normalized_cell_boundaries() const {
    std::vector<std::array<TF,2>> res( nb_diracs() );
    for_each_normalized_cell( [&]( TF dirac_position, TF dirac_weight, PI num_dirac, TF, TF, TF c0, TF c1 ) {
        const TF rd = sqrt( dirac_weight );
        if ( c0 > c1 ) {
            const TF x0 = max( dirac_position - rd, c1 );
            const TF x1 = min( dirac_position + rd, c0 );
            res[ num_dirac ] = { x0, x1 };
        } else {
            const TF x0 = max( dirac_position - rd, c0 );
            const TF x1 = min( dirac_position + rd, c1 );
            res[ num_dirac ] = { x0, x1 };
        }        
    } );
    return res;
}

DTP UTP::TV UTP::normalized_cell_barycenters() const {
    TV res( FromSize(), nb_diracs() );
    for_each_normalized_cell( [&]( TF dirac_position, TF dirac_weight, PI num_dirac, TF, TF, TF c0, TF c1 ) {
        const TF rd = sqrt( dirac_weight );
        if ( c0 > c1 ) {
            const TF x0 = max( dirac_position - rd, c1 );
            const TF x1 = min( dirac_position + rd, c0 );
            const TF it = density->integral( x0, x1 );
            res[ num_dirac ] = it ? - density->x_integral( x0, x1 ) / it : 0;
        } else {
            const TF x0 = max( dirac_position - rd, c0 );
            const TF x1 = min( dirac_position + rd, c1 );
            const TF it = density->integral( x0, x1 );
            res[ num_dirac ] = it ? density->x_integral( x0, x1 ) / it : 0;
        }        
    } );
    return res;
}

DTP UTP::TV UTP::cell_barycenters() const {
    const TF mul = ( end_x_density - beg_x_density ) / ( original_density_values.size() - 1 );
    TV bar = normalized_cell_barycenters();

    TV res( FromSize(), nb_diracs() );
    for( PI i = 0; i < nb_diracs(); ++i )
        res[ sorted_dirac_nums[ i ] ] = beg_x_density + mul * bar[ i ];
    return res;
}

DTP TF UTP::density_value( TF pos ) const {
    const TF mul = ( original_density_values.size() - 1 ) / ( end_x_density - beg_x_density );
    return density->value( mul * ( pos - beg_x_density ) );
}

DTP UTP::TV UTP::dirac_positions() const {
    const TF mul = ( end_x_density - beg_x_density ) / ( original_density_values.size() - 1 );
    TV res( FromSize(), nb_diracs() );
    for( PI i = 0; i < nb_diracs(); ++i )
        res[ sorted_dirac_nums[ i ] ] = beg_x_density + mul * sorted_dirac_positions[ i ];
    return res;
}

DTP TF UTP::normalized_error() const {
    using namespace std;

    TF res = 0;
    for_each_normalized_cell_mass( [&]( PI i, TF v, bool bad_cell ) {
        res += pow( v - sorted_dirac_masses[ i ], 2 );
    } );

    return sqrt( res );
}

DTP void UTP::plot( Str filename ) const {
    std::ofstream fs( filename );

    TF de = ( end_x_density - beg_x_density ) / 3;
    TV xs = Vec<TF>::linspace( beg_x_density - de, end_x_density + de, 10000 );
    TV ys;
    for( TF x : xs )
        ys << density_value( x );

    // for( PI i = 0; i < original_density_values.size(); ++i ) {
    //     const TF x = beg_x_density + ( end_x_density - beg_x_density ) * i / ( original_density_values.size() - 1 );
    //     ys.push_back( original_density_values[ i ] );
    //     xs.push_back( x );
    // }
    // ys.push_back( original_density_values.back() );
    // xs.push_back( end_x_density );



    TV bx = cell_barycenters();
    TV by( FromSizeAndItemValue(), bx.size(), -0.1 );

    TV dx = dirac_positions();
    TV dy( FromSizeAndItemValue(), dx.size(), -0.2 );

    fs << "from matplotlib import pyplot\n";
    fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
    fs << "pyplot.plot( " << to_string( bx ) << ", " << to_string( by ) << ", '+' )\n";
    fs << "pyplot.plot( " << to_string( dx ) << ", " << to_string( dy ) << ", '+' )\n";
    fs << "pyplot.show()\n";
}

template<class TF>
void plot_bnds_evolution( const Vec<Vec<Vec<TF,2>>> &bnds ) {
    auto ys = Vec<TF>::linspace( 0, 1, bnds.size() );
    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    for( PI i = 0; i < bnds[ 0 ].size(); ++i ) {
        for( PI j = 0; j < bnds[ 0 ][ 0 ].size(); ++j ) {
            Vec<TF> xs( FromSizeAndFunctionOnIndex(), bnds.size(), [&]( PI n ) { return bnds[ n ][ i ][ j ]; } );
            fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
        }
    }
    fs << "pyplot.show()\n";
}

// DTP int UTP::update_weights() {
//     Vec<Vec<Vec<TF,2>>> bnds;
//     for( PI num_iter = 0; ; ++num_iter ) {
//         bnds << normalized_cell_boundaries(); 
//         if ( num_iter == 100 ) {
//             plot_bnds_evolution( bnds );
//             throw std::runtime_error( "too many iterations" );
//         }

//         // line search
//         Vec<TF> base_dirac_weights = sorted_dirac_weights;
//         TF base_error = normalized_error();
//         Vec<TF> dir = newton_dir();

//         TF new_error = 0;
//         for( TF a = 1;; a /= 2 ) {
//             if ( a == 0 ) {
//                 plot_bnds_evolution( bnds );
//                 throw std::runtime_error( "bad direction" );
//             }
//             sorted_dirac_weights = base_dirac_weights + a * dir;
//             new_error = normalized_error();
//             if ( new_error < base_error ) {
//                 P( a, new_error );
//                 break;
//             }
//         } 

//         //
//         if ( new_error < target_l2_error ) {
//             plot_bnds_evolution( bnds );
//             P( num_iter );
//             return 0;
//         }
//     }

//     return 0;
// }

DTP void UTP::get_system( Vec<Vec<PI,2>> &connected_cells, Vec<Poly> &opt_weights, TF &max_a, TF &cell_error, int &has_bad_cell, const TV &sorted_dirac_weights ) const {
    using namespace std;

    //
    opt_weights.resize( nb_diracs() );
    connected_cells.clear();
    has_bad_cell = 0;
    cell_error = 0;
    max_a = 1;
    
    TF i0 = numeric_limits<TF>::lowest(); // prev cut position
    TF d0 = sorted_dirac_positions[ 0 ]; // prev dirac position
    TF w0 = sorted_dirac_weights[ 0 ]; // prev weight
    TF dp = 0; // prev prev dirac position
    TF wp = 0; // prev prev weight

    PI bi_ind; // index of the bi cell
    for( PI n = 0; n < nb_diracs(); ++n ) {
        // negative weight ?
        if ( w0 < 0 ) {
            has_bad_cell = 1;
            return;
        }

        // next intersection
        const TF d1 = n + 1 < nb_diracs() ? sorted_dirac_positions[ n + 1 ] : numeric_limits<TF>::max();
        const TF w1 = n + 1 < nb_diracs() ? sorted_dirac_weights[ n + 1 ] : 0;
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
                const TF sva = density->value( b0 ) + density->value( b1 );
                const TF itg = density->integral( b0, b1 );
                connected_cells.push_back_br( n, n );
                // if ( sva == 0 ) {
                //     has_bad_cell = 4;
                //     return;
                // }
                const TF mde = sorted_dirac_masses[ n ] / original_density_values.size();

                const TF err = sorted_dirac_masses[ n ] - itg;
                cell_error += pow( err, 2 );

                opt_weights[ n ] = { 
                    .c0 = r0,
                    .c1 = 0,
                    .c2 = 0,
                    .ca = err / max( mde, sva )
                };
            } else { // BI
                // opt_weights[ n ] is equal to r^2
                connected_cells.push_back_br( n, n );
                opt_weights[ n ] = { 
                    .c0 = 0,
                    .c1 = 0,
                    .c2 = 1,
                    .ca = 0
                };

                const TF err = density->integral( b0, i1 ) - sorted_dirac_masses[ n ];
                cell_error += pow( err, 2 );

                // mass of current cell will precribe opt_weights[ n + 1 ]
                //   ( 1 - a ) * integral( i0, i1 ) + a * target_area = integral( i0, i1 ) + v0 * ( r - r0 ) + v1 * 0.5 * ( r^2 - w0 - ow + w1 ) / ( d1 - d0 )
                //   ( 1 - a ) * integral( i0, i1 ) + a * target_area = integral( i0, i1 ) + v0 * ( r - r0 ) + ( r^2 - w0 - ow + w1 ) / C
                //   a * ( target_area - integral( i0, i1 ) ) - v0 * ( r - r0 ) = ( r^2 - w0 - ow + w1 ) / C
                //   ( a * ( target_area - integral( i0, i1 ) ) - v0 * ( r - r0 ) ) * C = r^2 - w0 - ow + w1
                //   ow = r^2 + w1 - w0 + ( a * ( integral( i0, i1 ) - target_area ) + v0 * ( r - r0 ) ) * C
                const TF dv0 = density->value( b0 );
                const TF dv1 = density->value( i1 );
                // if ( dv0 == 0 || dv1 == 0 ) {
                //     has_bad_cell = 5;
                //     return;
                // }

                const TF mde = sorted_dirac_masses[ n ] / original_density_values.size();
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
                // const TF v1 = density->value( b1 );
                // if ( v1 == 0 ) {
                //     has_bad_cell = 6;
                //     return;
                // }

                const TF err = density->integral( i0, b1 ) - sorted_dirac_masses[ n ];
                cell_error += pow( err, 2 );

                connected_cells.back()[ 1 ] = n;
            } else { // II
                const TF err = density->integral( i0, i1 ) - sorted_dirac_masses[ n ];
                cell_error += pow( err, 2 );

                // ( 1 - a ) * integral( i0, i1 ) + a * target_area = integral( i0, i1 ) + ( ow0 - w0 - owp + wp ) * C
                //                                                                       + ( ow0 - w0 - ow1 + w1 ) * D
                // a * ( target_area - integral( i0, i1 ) ) = ( ow0 - w0 - owp + wp ) * C
                //                                          + ( ow0 - w0 - ow1 + w1 ) * D
                // ow1 * D = a * ( integral( i0, i1 ) - target_area ) + ( ow0 - w0 - owp + wp ) ) * C
                //                                                    + ( ow0 - w0 + w1 ) * D
                const TF dv0 = density->value( i0 );
                const TF dv1 = density->value( i1 );
                // if ( dv1 == 0 ) {
                //     has_bad_cell = 7;
                //     return;
                // }

                const TF mde = sorted_dirac_masses[ n ] / original_density_values.size();
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

DTP Vec<TF> UTP::ce_errors( const Vec<Poly> &polys, PI nb, PI ne, TF a, TF r ) const {
    Vec<TF> res;
    res << bi_error( polys, nb, a, r );
    for( PI n = nb + 1; n < ne; ++n )
        res << ii_error( polys, n, a, r );    
    return res;
}

DTP TF UTP::ii_error( const Vec<Poly> &polys, PI n, TF a, TF r ) const {
    using namespace std;

    const TF dp = sorted_dirac_positions[ n - 1 ];
    const TF d0 = sorted_dirac_positions[ n + 0 ];
    const TF d1 = sorted_dirac_positions[ n + 1 ];

    const TF old_wp = sorted_dirac_weights[ n - 1 ];
    const TF old_w0 = sorted_dirac_weights[ n + 0 ];
    const TF old_w1 = sorted_dirac_weights[ n + 1 ];
    const TF old_i0 = ( dp + d0 + ( old_wp - old_w0 ) / ( d0 - dp ) ) / 2;
    const TF old_i1 = ( d0 + d1 + ( old_w0 - old_w1 ) / ( d1 - d0 ) ) / 2;

    const TF t_mass = ( 1 - a ) * density->integral( old_i0, old_i1 ) + a * sorted_dirac_masses[ n ];

    const TF new_wp = polys[ n - 1 ].c0 + polys[ n - 1 ].c1 * r + polys[ n - 1 ].c2 * pow( r, 2 ) + polys[ n - 1 ].ca * a;
    const TF new_w0 = polys[ n + 0 ].c0 + polys[ n + 0 ].c1 * r + polys[ n + 0 ].c2 * pow( r, 2 ) + polys[ n + 0 ].ca * a;
    const TF new_w1 = polys[ n + 1 ].c0 + polys[ n + 1 ].c1 * r + polys[ n + 1 ].c2 * pow( r, 2 ) + polys[ n + 1 ].ca * a;
    const TF new_i0 = ( dp + d0 + ( new_wp - new_w0 ) / ( d0 - dp ) ) / 2;
    const TF new_i1 = ( d0 + d1 + ( new_w0 - new_w1 ) / ( d1 - d0 ) ) / 2;

    const TF n_mass = density->integral( old_i0, old_i1 ) + density->value( old_i0 ) * ( old_i0 - new_i0 ) + density->value( old_i1 ) * ( new_i1 - old_i1 );

    return t_mass - n_mass;
}

DTP TF UTP::bi_error( const Vec<Poly> &polys, PI n, TF a, TF r ) const {
    using namespace std;

    const TF d0 = sorted_dirac_positions[ n + 0 ];
    const TF d1 = sorted_dirac_positions[ n + 1 ];

    const TF old_w0 = sorted_dirac_weights[ n + 0 ];
    const TF old_w1 = sorted_dirac_weights[ n + 1 ];
    const TF old_r0 = sqrt( old_w0 );
    const TF old_b0 = d0 - old_r0;
    const TF old_i1 = ( d0 + d1 + ( old_w0 - old_w1 ) / ( d1 - d0 ) ) / 2;

    const TF t_mass = ( 1 - a ) * density->integral( old_b0, old_i1 ) + a * sorted_dirac_masses[ n ];

    const TF new_w0 = polys[ n + 0 ].c0 + polys[ n + 0 ].c1 * r + polys[ n + 0 ].c2 * pow( r, 2 ) + polys[ n + 0 ].ca * a;
    const TF new_w1 = polys[ n + 1 ].c0 + polys[ n + 1 ].c1 * r + polys[ n + 1 ].c2 * pow( r, 2 ) + polys[ n + 1 ].ca * a;
    const TF new_r0 = sqrt( new_w0 );
    const TF new_b0 = d0 - new_r0;
    const TF new_i1 = ( d0 + d1 + ( new_w0 - new_w1 ) / ( d1 - d0 ) ) / 2;

    const TF n_mass = density->integral( old_b0, old_i1 ) + density->value( old_b0 ) * ( old_b0 - new_b0 ) + density->value( old_i1 ) * ( new_i1 - old_i1 );

    return t_mass - n_mass;
}

DTP TF UTP::best_r_for_ib( const Vec<Poly> &polys, PI n, TF a ) const {
    using namespace std;

    const TF dp = sorted_dirac_positions[ n - 1 ];
    const TF d0 = sorted_dirac_positions[ n + 0 ];

    const TF old_wp = sorted_dirac_weights[ n - 1 ];
    const TF old_w0 = sorted_dirac_weights[ n + 0 ];
    const TF old_i0 = ( dp + d0 + ( old_wp - old_w0 ) / ( d0 - dp ) ) / 2;
    const TF old_rd = sqrt( old_w0 );
    const TF old_b1 = d0 + old_rd;

    const TF old_in = density->integral( old_i0, old_b1 );
    const TF t_mass = ( 1 - a ) * old_in + a * sorted_dirac_masses[ n ];

    const TF v0 = density->value( old_i0 );
    const TF v1 = density->value( old_b1 );

    const TF pp0 = polys[ n - 1 ].c0, pp1 = polys[ n - 1 ].c1, pp2 = polys[ n - 1 ].c2, ppa = polys[ n - 1 ].ca;
    const TF p00 = polys[ n + 0 ].c0, p01 = polys[ n + 0 ].c1, p02 = polys[ n + 0 ].c2, p0a = polys[ n + 0 ].ca;

    const TF ter = target_mass_error * sorted_dirac_masses[ n ] * 1e-5;

    auto err = [&]( const TF r ) {
        const TF new_wp = pp0 + pp1 * r + pp2 * pow( r, 2 ) + ppa * a;
        const TF new_w0 = p00 + p01 * r + p02 * pow( r, 2 ) + p0a * a;
        const TF new_i0 = ( dp + d0 + ( new_wp - new_w0 ) / ( d0 - dp ) ) / 2;
        const TF new_rd = sqrt( max( 0, new_w0 ) );
        const TF new_b1 = d0 + new_rd;

        const TF n_mass = old_in + v0 * ( old_i0 - new_i0 ) + v1 * ( new_b1 - old_b1 );
        
        return n_mass - t_mass;
    };

    if ( p02 ) {
        const TF b = pow( p01, 2 ) - 4 * p02 * ( p00 + p0a * a );
        if ( b >= 0 ) {
            const TF xa = ( - p01 - sqrt( b ) ) / ( 2 * p02 );
            const TF xb = ( - p01 + sqrt( b ) ) / ( 2 * p02 );
            const TF x0 = max( 0, min( xa, xb ) );
            const TF x1 = max( 0, max( xa, xb ) );
            if ( p02 < 0 ) { // solution must be in [ x0, x1 ]
                if ( x0 > 0 && abs( err( x0 ) ) <= ter )
                    return x0;
                if ( x1 > 0 && abs( err( x1 ) ) <= ter )
                    return x1;
                if ( sgn( err( x0 ) ) == sgn( err( x1 ) ) )
                    return -2;
                return dichotomy( err, ter, max( 0, x0 ), x1 );
            }
            
            // solution may be in [ 0, x0 ]
            if ( x0 > 0 ) {
                if ( abs( err( 0 ) ) <= ter )
                    return 0;
                if ( abs( err( x0 ) ) <= ter )
                    return x0;
                if ( sgn( err( 0 ) ) != sgn( err( x0 ) ) )
                    return dichotomy( err, ter, TF( 0 ), x0 );
            }

            // else, solution must be in [ x1, inf ]
            if ( x1 > 0 && abs( err( x1 ) ) <= ter )
                return x1;
            if ( err( x1 ) > 0 )
                return -3;
 
            TF max_r = x1 ? 2 * x1 : 1;
            while ( err( max_r ) < 0 ) {
                if ( max_r > original_density_values.size() )
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
            if ( max_r > original_density_values.size() )
                return -7;
            max_r *= 2;
        }
        // P( __LINE__, max_r );
        return dichotomy( err, ter, TF( 0 ), max_r );
    }

    TODO;
    return 0;
}

DTP int UTP::get_weights_for( Vec<TF> &new_dirac_weights, const Vec<Vec<PI,2>> &connected_cells, const Vec<Poly> &polys, TF a ) {
    new_dirac_weights.resize( nb_diracs() );
    for( Vec<PI,2> inds : connected_cells ) {
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

DTP void UTP::solve() {
    if ( nb_diracs() == 0 )
        return;

    // system content
    Vec<Vec<PI,2>> new_connected_cells;
    Vec<Vec<PI,2>> connected_cells;
    Vec<Poly> new_opt_weights;
    Vec<Poly> opt_weights;
    int has_bad_cell;
    TF new_cell_error;
    TF cell_error;
    TF max_a;

    // history
    // Vec<Vec<Vec<TF,2>>> bnds;

    // get the first system
    get_system( connected_cells, opt_weights, max_a, cell_error, has_bad_cell, sorted_dirac_weights );
    if ( has_bad_cell )
        throw std::runtime_error( "bad initialization" );
    // P( cell_error );
    if ( cell_error < target_mass_error / nb_diracs() )
        return;

    // iterate
    Vec<TF> new_sorted_dirac_weights;
    for( PI num_iter = 0; ; ++num_iter ) {
        if ( num_iter == 1000 )
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
                    std::swap( new_sorted_dirac_weights, sorted_dirac_weights );
                    std::swap( new_connected_cells, connected_cells );
                    std::swap( new_opt_weights, opt_weights );
                    std::swap( new_cell_error, cell_error );
    
                    // P( cell_error, a );
                    if ( sqrt( cell_error ) < target_mass_error / nb_diracs() ) {
                        P( num_iter );
                        return;
                    }
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
