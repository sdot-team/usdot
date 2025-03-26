#pragma once

// #include <tl/support/operators/self_max.h>
// #include <tl/support/operators/norm_2.h>
// #include <tl/support/operators/argmin.h>
// #include <tl/support/operators/max.h>
// #include <tl/support/operators/sum.h>
#include <algorithm>
#include <stdexcept>
#include <tl/support/ASSERT.h>
#include <tl/support/P.h>
#include <numeric>
#include <fstream>
#include <limits>

// #include "partial1D/Extrapolation.h"
#include "GradientSolver.h"
#include "glot.h"

namespace usdot {

    
#define DTP template<class TF>
#define UTP GradientSolver<TF>

DTP UTP::GradientSolver( GradientSolverInput<TF> &&input ) {
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
    target_l2_error = input.target_l2_error;
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
    const TF weight = 1 * pow( TF( original_density_values.size() - 1 ) / nb_diracs, 2 );
    const TF sep = input.min_dirac_separation * mul_dirac;
    sorted_dirac_positions.resize( nb_diracs );
    sorted_dirac_weights.resize( nb_diracs );
    sorted_dirac_masses.resize( nb_diracs );
    sum_of_dirac_masses = 0;
    for( PI i = 0; i < nb_diracs; ++i ) {
        sorted_dirac_positions[ i ] = ( input.dirac_positions[ sorted_dirac_nums[ i ] ] - beg_x_density ) * mul_dirac;
        if ( i && sorted_dirac_positions[ i ] < sorted_dirac_positions[ i - 1 ] + sep ) 
            sorted_dirac_positions[ i ] = sorted_dirac_positions[ i - 1 ] + sep;

        const TF dirac_weight = sorted_dirac_nums[ i ] < input.dirac_weights.size() ? input.dirac_weights[ sorted_dirac_nums[ i ] ] * pow( mul_dirac, 2 ) : weight;
        sorted_dirac_weights[ i ] = dirac_weight;
    
        const TF dirac_mass = sorted_dirac_nums[ i ] < input.dirac_masses.size() ? input.dirac_masses[ sorted_dirac_nums[ i ] ] : 1;
        sorted_dirac_masses[ i ] = dirac_mass;
        sum_of_dirac_masses += dirac_mass;
    }

    // density_gaussian_width
    initialize_filter_value( input.starting_filter_value );
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

DTP void UTP::for_each_normalized_system_item( auto &&func ) const {
    using namespace std;

    for_each_normalized_cell( [&]( TF dirac_position, TF dirac_weight, PI num_dirac, TF d0, TF d1, TF c0, TF c1 ) {
        const TF v = sorted_dirac_masses[ num_dirac ];
        if ( dirac_weight <= 0 )
            return func( num_dirac, 0, 0, v, true );

        const TF rd = sqrt( dirac_weight );
        const TF b0 = dirac_position - rd;
        const TF b1 = dirac_position + rd;

        // "wrong" bounds ordering
        if ( c0 > c1 ) {
            // void cell
            if ( c0 < b0 || c1 > b1 )
                return func( num_dirac, 0, 0, v, true );

            // ball cut on the left
            if ( b0 > c1 ) {
                // ball cut on the right
                if ( b1 < c0 ) {
                    return func( num_dirac, 
                        ( density->value( b0 ) + density->value( b1 ) ) / ( 2 * sqrt( dirac_weight ) ),
                        0, 
                        - density->integral( b0, b1 ),
                        false
                    );
                }

                // interface on the right
                return func( num_dirac, 
                    ( density->value( c0 ) / d0 - density->value( b0 ) / sqrt( dirac_weight ) ) / 2,
                    0,
                    - density->integral( b0, c0 ),
                    false
                );
            }
            
            // interface on the left
            if ( b1 < c0 ) { // ball cut on the right
                const TF cl = density->value( c1 ) / ( 2 * d1 );
                return func( num_dirac, 
                    cl - density->value( b1 ) / ( 2 * sqrt( dirac_weight ) ),
                    - cl,
                    - density->integral( c1, b1 ),
                    false
                );
            }
            
            // interface on the right
            const TF cl = density->value( c1 ) / ( 2 * d1 );
            const TF cr = density->value( c0 ) / ( 2 * d0 );
            return func( num_dirac, 
                - ( cl + cr ),
                cr,
                - density->integral( c1, c0 ),
                false
            );
        }

        // void cell
        if ( c1 < b0 || c0 > b1 )
            return func( num_dirac, 0, 0, v, true );

        if ( b0 > c0 ) { // ball cut on the left
            // ball cut on the right
            if ( b1 < c1 ) { 
                return func( num_dirac, 
                    ( density->value( b0 ) + density->value( b1 ) ) / ( 2 * sqrt( dirac_weight ) ),
                    0,
                    density->integral( b0, b1 ),
                    false
                );
            }

            // interface on the right
            const TF cr = density->value( c1 ) / ( 2 * d1 );
            return func( num_dirac, 
                density->value( b0 ) / ( 2 * sqrt( dirac_weight ) ) + cr,
                - cr,
                density->integral( b0, c1 ),
                false
            );
        } 
        
        // interface on the left
        if ( b1 < c1 ) { // ball cut on the right
            const TF cl = density->value( c0 ) / ( 2 * d0 );
            return func( num_dirac, 
                cl + density->value( b1 ) / ( 2 * sqrt( dirac_weight ) ),
                0,
                density->integral( c0, b1 ),
                false
            );
        }
        
        const TF cl = density->value( c0 ) / ( 2 * d0 );
        const TF cr = density->value( c1 ) / ( 2 * d1 );
        const TF in = density->integral( c0, c1 );
        return func( num_dirac, cl + cr, - cr, in, false );
    } );
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
    TV res( nb_diracs() );
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

    TV res( nb_diracs() );
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
    TV res( nb_diracs() );
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

    TV xs;
    TV ys;
    for( PI i = 0; i < original_density_values.size(); ++i ) {
        const TF x = beg_x_density + ( end_x_density - beg_x_density ) * i / ( original_density_values.size() - 1 );
        ys.push_back( original_density_values[ i ] );
        xs.push_back( x );
    }
    ys.push_back( original_density_values.back() );
    xs.push_back( end_x_density );

    TV bx = cell_barycenters();
    TV by( bx.size(), -0.1 );

    TV dx = dirac_positions();
    TV dy( dx.size(), -0.2 );

    fs << "from matplotlib import pyplot\n";
    fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
    fs << "pyplot.plot( " << to_string( bx ) << ", " << to_string( by ) << ", '+' )\n";
    fs << "pyplot.plot( " << to_string( dx ) << ", " << to_string( dy ) << ", '+' )\n";
    fs << "pyplot.show()\n";
}

// DTP std::pair<typename UTP::TV,TF> UTP::newton_dir() const {
//     bool has_bad_cell = false;
//     TV x( nb_diracs() );
//     TV e( nb_diracs() );
//     TF prev_dll = 0;
//     TF prev_v = 0;
//     TF error = 0;
//     for_each_normalized_system_item( [&]( PI index, TF m0, TF m1, TF v, bool bad_cell ) {
//         v -= sorted_dirac_masses[ index ];
//         has_bad_cell |= bad_cell;
//         error += pow( v, 2 );

//         const TF d = m0 - prev_dll;
//         const TF l = m1 / d;
//         e[ index ] = l;

//         prev_dll = d * l * l;

//         const TF v0 = v - prev_v;
//         x[ index ] = v0 / d;
//         prev_v = l * v0;
//     } );

//     for( PI i = nb_diracs() - 1; i--; )
//         x[ i ] -= e[ i ] * x[ i + 1 ];

//     return { x, has_bad_cell ? std::numeric_limits<TF>::max() : error };
// }

DTP int UTP::update_weights() {
    //  
    TV g0( FromSize(), nb_diracs() );
    for_each_normalized_cell_mass( [&]( PI index, TF mass, bool bad_cell ) {
        g0[ index ] = sorted_dirac_masses[ index ] - mass;
    } );

    Vec<TV> gradients{ std::move( g0 ) };
    for( PI i = 1; i < floor( log2( nb_diracs() ) ); ++i ) {
        TV ng( FromSize(), nb_diracs() );
        TF f = 1 - pow( 0.1, i );
        TF prev = 0;
        for( PI j = 0; j < nb_diracs(); ++j ) {
            const TF n = f * prev + ( 1 - f ) * gradients[ 0 ][ j ];
            ng[ j ] = n;
            prev = n;
        }
        prev = 0;
        for( PI j = nb_diracs(); j--; ) {
            const TF n = f * prev + ( 1 - f ) * ng[ j ];
            ng[ j ] = n;
            prev = n;
        }

        const TF m = max( ng );
        ng /= m;

        gradients << std::move( ng );
    }

    P( gradients.size() );

    glot_vec( Vec<TF>::linspace( 0, 1, nb_diracs() ), gradients[ 0 ], gradients[ 1 ], gradients[ 5 ], gradients[ 6 ] );
    
    return 0;
}

DTP void UTP::solve() {
    if ( nb_diracs() == 0 )
        return;

    // solve
    update_weights();

}

#undef DTP
#undef UTP


} // namespace usdot
