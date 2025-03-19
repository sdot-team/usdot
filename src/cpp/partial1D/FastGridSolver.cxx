#pragma once

#include <tl/support/operators/self_max.h>
#include <tl/support/operators/norm_2.h>
#include <tl/support/operators/max.h>
#include <tl/support/operators/sum.h>
#include <tl/support/ASSERT.h>
#include <tl/support/P.h>
#include <numbers>
#include <fstream>
#include <limits>

// #include "SymmetricBandMatrix.h"
#include "FastGridSolver.h"
#include "ThreadPool.h"

namespace usdot {

    
#define DTP template<class TF>
#define UTP FastGridSolver<TF>

DTP UTP::FastGridSolver( FastGridSolverInput<TF> &&input ) {
    using namespace std;

    // check inputs
    if ( input.dirac_positions.empty() )
        return;

    if ( input.density_values.empty() )
        throw runtime_error( "one needs at least 1 density value" );

    // values that can be copied
    target_contrast_ratio = input.target_contrast_ratio;
    global_mass_ratio = input.global_mass_ratio;
    target_l2_error = input.target_l2_error;
    throw_if_error = input.throw_if_error;
    multithread = input.multithread;
    verbosity = input.verbosity;

    // values that can be moved
    original_density_values = std::move( input.density_values );

    // density bounds
    end_x_density = input.end_x_density ? *input.end_x_density : 1;
    beg_x_density = input.beg_x_density ? *input.beg_x_density : 0;

    // sorted_dirac_nums
    const PI nb_diracs = input.dirac_positions.size();
    sorted_dirac_nums = { FromSizeAndFunctionOnIndex(), nb_diracs, []( auto i ) { return i; } };
    sort( sorted_dirac_nums.begin(), sorted_dirac_nums.end(), [&]( PI a, PI b ) {
        return input.dirac_positions[ a ] < input.dirac_positions[ b ];
    } );

    // sorted_dirac_positions, sorted_dirac_weights, sorted_dirac_masses, sum_of_dirac_masses
    const TF mul_dirac = original_density_values.size() / ( end_x_density - beg_x_density );
    const TF weight = 2 * pow( TF( original_density_values.size() ) / nb_diracs, 2 );
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
    normalized_density_gaussian_width = max( sorted_dirac_positions.back (), end_x_density ) - 
                                        min( sorted_dirac_positions.front(), beg_x_density ) ;

    // density_values
    max_of_original_density_values = max( original_density_values );
    set_density_contrast( input.starting_contrast_ratio );
}

DTP PI UTP::nb_diracs() const {
    return sorted_dirac_positions.size();
}

DTP void UTP::for_each_normalized_cell_mt( auto &&func ) const {
    using namespace std;

    PI nb_jobs = 1 * thread_pool.nb_threads();
    thread_pool.execute( nb_jobs, [&]( PI num_job, PI num_thread ) {
        const PI beg = ( num_job + 0 ) * nb_diracs() / nb_jobs;
        const PI end = ( num_job + 1 ) * nb_diracs() / nb_jobs;
        if ( beg == end )
            return;

        //
        TF i0 = numeric_limits<TF>::lowest(); // previous interface
        TF o0 = numeric_limits<TF>::max(); // previous distance
        TF d0 = sorted_dirac_positions[ beg ];
        TF w0 = sorted_dirac_weights[ beg ];
        if ( beg ) {
            const TF dp = sorted_dirac_positions[ beg - 1 ];
            const TF wp = sorted_dirac_weights[ beg - 1 ];
            i0 = ( dp + d0 + ( wp - w0 ) / ( d0 - dp ) ) / 2;
            o0 = d0 - dp;    
        }

        for( PI n = beg + 1; n < end; ++n ) {
            const TF d1 = sorted_dirac_positions[ n ];
            const TF w1 = sorted_dirac_weights[ n ];
    
            const TF i1 = ( d0 + d1 + ( w0 - w1 ) / ( d1 - d0 ) ) / 2;
            func( d0, w0, n - 1, o0, d1 - d0, i0, i1, num_thread );
    
            o0 = d1 - d0;
            i0 = i1;
            d0 = d1;
            w0 = w1;
        }

        TF i1 = numeric_limits<TF>::max();
        TF o1 = numeric_limits<TF>::max();
        if ( end < nb_diracs() ) {
            const TF dn = sorted_dirac_positions[ end ];
            const TF wn = sorted_dirac_weights[ end ];
            i1 = ( d0 + dn + ( w0 - wn ) / ( dn - d0 ) ) / 2;
            o1 = dn - d0;
        }

        func( d0, w0, end - 1, o0, o1, i0, i1, num_thread );
    } );
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

DTP TF UTP::normalized_density_x_integral_ap( TF x0, TF x1, PI nb_steps ) const {
    TF res = 0;
    for( TF x : Vec<TF>::linspace( x0, x1, nb_steps ) )
        res += x * normalized_density_value( x );
    return res * ( x1 - x0 ) / nb_steps;
}

DTP TF UTP::normalized_density_x_primitive( TF x ) const {
    PI s = normalized_density_values.size();
    using namespace std;

    if ( x >= s ) {
        const TF dx = ( s - x ) / normalized_density_gaussian_width;
        const TF ex = normalized_density_gaussian_width * ( 1 - exp( - pow( dx, 2 ) ) )
                    -  sqrt( numbers::pi_v<TF> ) * s * erf( dx );
        return normalized_density_x_primitives.back() + normalized_density_gaussian_width / 2 * min_density_value * ex;
    }

    if ( x < 0 ) { 
        const TF dx = x / normalized_density_gaussian_width;
        return min_density_value * pow( normalized_density_gaussian_width, 2 ) * ( 1 - exp( - pow( dx, 2 ) ) ) / 2;
    }

    // integral( ( ix + x ) * density_values[ ix ], x, 0, fx )
    // density_values[ ix ] * ( ix * fx + fx^2/2 )
    PI ix( x );
    TF fx = x - ix;
    return normalized_density_x_primitives[ ix ] + normalized_density_values[ ix ] * fx * ( ix + fx / 2 );
}

DTP TF UTP::normalized_density_x_integral( TF x0, TF x1 ) const {
    return normalized_density_x_primitive( x1 ) - normalized_density_x_primitive( x0 );
}

DTP TF UTP::normalized_density_primitive( TF x ) const {
    PI s = normalized_density_values.size();
    using namespace std;

    if ( x >= s )
        return normalized_density_primitives.back() + min_density_value * erf( ( x - s ) / normalized_density_gaussian_width ) * normalized_density_gaussian_width * sqrt( numbers::pi_v<TF> ) / 2;

    if ( x < 0 ) 
        return min_density_value * erf( x / normalized_density_gaussian_width ) * normalized_density_gaussian_width * sqrt( numbers::pi_v<TF> ) / 2;

    
    PI ix( x );
    TF fx = x - ix;
    return normalized_density_primitives[ ix ] + fx * normalized_density_values[ ix ];
}

DTP TF UTP::normalized_density_integral( TF x0, TF x1 ) const {
    return normalized_density_primitive( x1 ) - normalized_density_primitive( x0 );
}

DTP TF UTP::normalized_density_value( TF x ) const {
    PI s = normalized_density_values.size();
    using namespace std;

    if ( x >= s ) return min_density_value * exp( - pow( ( x - normalized_density_values.size() ) / normalized_density_gaussian_width, 2 ) );
    if ( x < 0 ) return min_density_value * exp( - pow( x / normalized_density_gaussian_width, 2 ) );
    return normalized_density_values[ x ];
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
                        ( normalized_density_value( b0 ) + normalized_density_value( b1 ) ) / sqrt( dirac_weight ), 
                        0, 
                        2 * ( v + normalized_density_integral( b0, b1 ) ),
                        false
                    );
                }

                // interface on the right
                return func( num_dirac, 
                    normalized_density_value( c0 ) / d0 - normalized_density_value( b0 ) / sqrt( dirac_weight ),
                    0,
                    2 * ( v + normalized_density_integral( b0, c0 ) ),
                    false
                );
            }
            
            // interface on the left
            if ( b1 < c0 ) { // ball cut on the right
                const TF cl = normalized_density_value( c1 ) / d1;
                return func( num_dirac, 
                    cl - normalized_density_value( b1 ) / sqrt( dirac_weight ),
                    - cl,
                    2 * ( v + normalized_density_integral( c1, b1 ) ),
                    false
                );
            }
            
            // interface on the right
            const TF cl = normalized_density_value( c1 ) / d1;
            const TF cr = normalized_density_value( c0 ) / d0;
            return func( num_dirac, 
                - ( cl + cr ),
                cr,
                2 * ( v + normalized_density_integral( c1, c0 ) ),
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
                    ( normalized_density_value( b0 ) + normalized_density_value( b1 ) ) / sqrt( dirac_weight ),
                    0,
                    2 * ( v - normalized_density_integral( b0, b1 ) ),
                    false
                );
            }

            // interface on the right
            const TF cr = normalized_density_value( c1 ) / d1;
            return func( num_dirac, 
                normalized_density_value( b0 ) / sqrt( dirac_weight ) + cr,
                - cr,
                2 * ( v - normalized_density_integral( b0, c1 ) ),
                false
            );
        } 
        
        // interface on the left
        if ( b1 < c1 ) { // ball cut on the right
            const TF cl = normalized_density_value( c0 ) / d0;
            return func( num_dirac, 
                cl + normalized_density_value( b1 ) / sqrt( dirac_weight ),
                0,
                2 * ( v - normalized_density_integral( c0, b1 ) ),
                false
            );
        }
        
        // // interface on the right
        // PI s = normalized_density_values.size();
        // TF vc0, pc0;
        // if ( c0 >= s ) {
        //     vc0 = min_density_value * exp( - pow( ( c0 - normalized_density_values.size() ) / normalized_density_gaussian_width, 2 ) );
        //     pc0 = normalized_density_primitives.back() + min_density_value * erf( ( c0 - s ) / normalized_density_gaussian_width ) * normalized_density_gaussian_width * sqrt( numbers::pi_v<TF> ) / 2;
        // }
        // else if ( c0 < 0 ) {
        //     vc0 = min_density_value * exp( - pow( c0 / normalized_density_gaussian_width, 2 ) );
        //     pc0 = min_density_value * erf( c0 / normalized_density_gaussian_width ) * normalized_density_gaussian_width * sqrt( numbers::pi_v<TF> ) / 2;
        // }
        // else { 
        //     PI ix( c0 );
        //     TF fx = c0 - ix;
        //     vc0 = normalized_density_values[ c0 ];
        //     pc0 = normalized_density_primitives[ ix ] + fx * normalized_density_values[ ix ];
        // }

        // TF vc1, pc1;
        // if ( c1 >= s ) {
        //     vc1 = min_density_value * exp( - pow( ( c1 - normalized_density_values.size() ) / normalized_density_gaussian_width, 2 ) );
        //     pc1 = normalized_density_primitives.back() + min_density_value * erf( ( c1 - s ) / normalized_density_gaussian_width ) * normalized_density_gaussian_width * sqrt( numbers::pi_v<TF> ) / 2;
        // }
        // else if ( c1 < 0 ) {
        //     vc1 = min_density_value * exp( - pow( c1 / normalized_density_gaussian_width, 2 ) );
        //     pc1 = min_density_value * erf( c1 / normalized_density_gaussian_width ) * normalized_density_gaussian_width * sqrt( numbers::pi_v<TF> ) / 2;
        // }
        // else { 
        //     PI ix( c1 );
        //     TF fx = c1 - ix;
        //     vc1 = normalized_density_values[ c1 ];
        //     pc1 = normalized_density_primitives[ ix ] + fx * normalized_density_values[ ix ];
        // }
        
        // const TF cl = vc0 / d0;
        // const TF cr = vc1 / d1;
        // return func( num_dirac, 
        //     cl + cr,
        //     - cr,
        //     2 * ( v - ( pc1 - pc0 ) ),
        //     false
        // );
        // 6763 us

        const TF cl = normalized_density_value( c0 ) / d0;
        const TF cr = normalized_density_value( c1 ) / d1;
        return func( num_dirac, 
            cl + cr,
            - cr,
            2 * ( v - normalized_density_integral( c0, c1 ) ),
            false
        );
        // 6868
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
            return func( num_dirac, normalized_density_integral( a0, a1 ), a1 > a0 );
        }

        const TF a0 = max( b0, c0 );
        const TF a1 = min( b1, c1 );
        return func( num_dirac, normalized_density_integral( a0, a1 ), a0 > a1 );
    } );
}

DTP TF UTP::l2_error() const {
    using namespace std;

    TF res = 0;
    bool has_bad_cell = 0;
    for_each_normalized_cell_mass( [&]( PI index, TF v, bool bad_cell ) {
        res += pow( sorted_dirac_masses[ index ] - v, 2 );
        has_bad_cell |= bad_cell;
    } );

    return has_bad_cell ? numeric_limits<TF>::max() : sqrt( res );
}

DTP Vec<TF> UTP::newton_dir() {    
    // integrated LDL solver
    Vec<TF> x( FromSize(), nb_diracs() );
    Vec<TF> e( FromSize(), nb_diracs() );
    TF prev_ldds = 0;
    TF prev_v = 0;
    for_each_normalized_system_item( [&]( PI index, TF m0, TF m1, TF v, bool bad_cell ) {
        if ( index == 0 && global_mass_ratio == 1 && current_contrast_ratio == 0 )
            m0 += 1;
    
        const TF d = m0 - prev_ldds;
        const TF l = m1 / d;
        e[ index ] = l;

        prev_ldds = d * l * l;

        const TF v0 = v - prev_v;
        x[ index ] = v0 / d;
        prev_v = l * v0;
    } );

    for( PI i = nb_diracs() - 1; i--; )
        x[ i ] -= e[ i ] * x[ i + 1 ];
        
    return x;
}


DTP void UTP::plot( Str filename ) const {
    std::ofstream fs( filename );

    const PI s = normalized_density_values.size();
    Vec<TF> xs = Vec<TF>::linspace( - 0.1 * s, 1.1 * s, 20 * s );
    Vec<TF> ys = map_vec( xs, [&]( TF x ) { return normalized_density_value( x ); } );

    Vec<TF> bx = normalized_cell_barycenters();
    Vec<TF> by = bx * 0;

    fs << "from matplotlib import pyplot\n";
    fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
    fs << "pyplot.plot( " << to_string( bx ) << ", " << to_string( by ) << ", '+' )\n";
    fs << "pyplot.show()\n";
}

// DTP TF UTP::dirac_convolution( TF normalized_pos, TF convolution_width ) const {
//     using namespace std;

//     TF res = 0;
//     for( PI i = 0; i < nb_diracs(); ++i )
//         res += sorted_dirac_masses[ i ] * exp( - pow( ( normalized_pos - sorted_dirac_positions[ i ] ) / convolution_width, 2 ) );
//     return res;
// }

DTP Vec<TF> UTP::normalized_cell_barycenters() const {
    Vec<TF> res( FromSize(), nb_diracs() );
    for_each_normalized_cell( [&]( TF dirac_position, TF dirac_weight, PI num_dirac, TF, TF, TF c0, TF c1 ) {
        const TF rd = sqrt( dirac_weight );
        if ( c0 > c1 ) {
            const TF x0 = max( dirac_position - rd, c1 );
            const TF x1 = min( dirac_position + rd, c0 );
            const TF it = normalized_density_integral( x0, x1 );
            res[ num_dirac ] = it ? - normalized_density_x_integral( x0, x1 ) / it : 0;
        } else {
            const TF x0 = max( dirac_position - rd, c0 );
            const TF x1 = min( dirac_position + rd, c1 );
            const TF it = normalized_density_integral( x0, x1 );
            res[ num_dirac ] = it ? normalized_density_x_integral( x0, x1 ) / it : 0;
        }        
    } );
    return res;
}

DTP void UTP::set_density_contrast( TF max_ratio ) {    
    // first take
    min_density_value = max_of_original_density_values * max_ratio;
    current_contrast_ratio = max_ratio;
    need_lower_contrast_ratio = false;
    
    const PI total_size = original_density_values.size();
    normalized_density_x_primitives.resize( total_size + 1 );
    normalized_density_primitives.resize( total_size + 1 );
    normalized_density_values.resize( total_size );
    
    TF integral_of_density_x_values = 0;
    TF integral_of_density_values = 0;
    for( PI i = 0; i < original_density_values.size(); ++i ) {
        TF dv = original_density_values[ i ];
        if ( dv < min_density_value ) {
            need_lower_contrast_ratio = true;
            dv = min_density_value;
        }

        normalized_density_x_primitives[ i ] = integral_of_density_x_values;
        normalized_density_primitives[ i ] = integral_of_density_values;
        normalized_density_values[ i ] = dv;

        integral_of_density_x_values += dv * ( i + 0.5 );
        integral_of_density_values += dv;
    }

    normalized_density_x_primitives.back() = integral_of_density_x_values;
    normalized_density_primitives.back() = integral_of_density_values;

    // erf on the left and on the right
    // sum_of_density_values += 2 * min_density_value * density_values.size();

    // normalization
    const TF coeff = sum_of_dirac_masses / ( integral_of_density_values * global_mass_ratio );
    min_density_value *= coeff;
    for( TF &v : normalized_density_x_primitives )
        v *= coeff;
    for( TF &v : normalized_density_primitives )
        v *= coeff;
    for( TF &v : normalized_density_values )
        v *= coeff;
}

// DTP bool UTP::converged( const Errors &er ) const {
//     return er.ni_residual < max_mass_ratio_error_target;
// }

DTP bool UTP::update_weights() {
    using namespace std;
    TF old_error = l2_error();
    for( PI num_iter = 0; num_iter < 50; ++num_iter ) {
        // get a new direction
        Vec<TF> dir = newton_dir();
        
        // find a first relaxation coeff
        Vec<TF> old_dirac_weights = sorted_dirac_weights;
        for( TF a = 1; ; a /= 2 ) {
            if ( a <= 1e-20 )
                throw runtime_error( "bad direction" );
            
            sorted_dirac_weights = old_dirac_weights + a * dir;
            TF new_error = l2_error();
            
            if ( new_error < old_error ) {
                P( num_iter, a, new_error );
                if ( new_error < target_l2_error )
                    return true;

                old_error = new_error;
                break;
            }
        }
    }

    return false;
}

DTP void UTP::solve() {
    update_weights();

    if ( target_contrast_ratio != current_contrast_ratio ) {
        set_density_contrast( target_contrast_ratio );
        update_weights();
    }
}

#undef DTP
#undef UTP


} // namespace usdot
