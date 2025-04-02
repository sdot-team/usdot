#pragma once

#include <tl/support/operators/self_max.h>
#include <tl/support/operators/norm_2.h>
#include <tl/support/operators/argmin.h>
#include <tl/support/operators/max.h>
#include <tl/support/operators/sum.h>
#include <tl/support/ASSERT.h>
#include <tl/support/P.h>
#include <numbers>
#include <fstream>
#include <limits>

#include "LogGridSolver.h"
#include "glot.h"

namespace usdot {

    
#define DTP template<class TF>
#define UTP LogGridSolver<TF>

DTP UTP::LogGridSolver( LogGridSolverInput<TF> &&input ) {
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
                        ( normalized_density_value( b0 ) + normalized_density_value( b1 ) ) / ( 2 * sqrt( dirac_weight ) ),
                        0, 
                        - normalized_density_integral( b0, b1 ),
                        false
                    );
                }

                // interface on the right
                return func( num_dirac, 
                    ( normalized_density_value( c0 ) / d0 - normalized_density_value( b0 ) / sqrt( dirac_weight ) ) / 2,
                    0,
                    - normalized_density_integral( b0, c0 ),
                    false
                );
            }
            
            // interface on the left
            if ( b1 < c0 ) { // ball cut on the right
                const TF cl = normalized_density_value( c1 ) / ( 2 * d1 );
                return func( num_dirac, 
                    cl - normalized_density_value( b1 ) / ( 2 * sqrt( dirac_weight ) ),
                    - cl,
                    - normalized_density_integral( c1, b1 ),
                    false
                );
            }
            
            // interface on the right
            const TF cl = normalized_density_value( c1 ) / ( 2 * d1 );
            const TF cr = normalized_density_value( c0 ) / ( 2 * d0 );
            return func( num_dirac, 
                - ( cl + cr ),
                cr,
                - normalized_density_integral( c1, c0 ),
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
                    ( normalized_density_value( b0 ) + normalized_density_value( b1 ) ) / ( 2 * sqrt( dirac_weight ) ),
                    0,
                    normalized_density_integral( b0, b1 ),
                    false
                );
            }

            // interface on the right
            const TF cr = normalized_density_value( c1 ) / ( 2 * d1 );
            return func( num_dirac, 
                normalized_density_value( b0 ) / ( 2 * sqrt( dirac_weight ) ) + cr,
                - cr,
                normalized_density_integral( b0, c1 ),
                false
            );
        } 
        
        // interface on the left
        if ( b1 < c1 ) { // ball cut on the right
            const TF cl = normalized_density_value( c0 ) / ( 2 * d0 );
            return func( num_dirac, 
                cl + normalized_density_value( b1 ) / ( 2 * sqrt( dirac_weight ) ),
                0,
                normalized_density_integral( c0, b1 ),
                false
            );
        }
        
        const TF cl = normalized_density_value( c0 ) / ( 2 * d0 );
        const TF cr = normalized_density_value( c1 ) / ( 2 * d1 );
        const TF in = normalized_density_integral( c0, c1 );
        return func( num_dirac, cl + cr, - cr, in, false );
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

DTP Vec<Vec<TF,2>> UTP::normalized_cell_boundaries() const {
    Vec<Vec<TF,2>> res( FromSize(), nb_diracs() );
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

DTP Vec<TF> UTP::cell_barycenters() const {
    const TF mul = ( end_x_density - beg_x_density ) / original_density_values.size();
    Vec<TF> bar = normalized_cell_barycenters();

    Vec<TF> res( FromSize(), nb_diracs() );
    for( PI i = 0; i < nb_diracs(); ++i )
        res[ sorted_dirac_nums[ i ] ] = beg_x_density + mul * bar[ i ];
    return res;
}

DTP Vec<TF> UTP::dirac_positions() const {
    const TF mul = ( end_x_density - beg_x_density ) / original_density_values.size();
    Vec<TF> res( FromSize(), nb_diracs() );
    for( PI i = 0; i < nb_diracs(); ++i )
        res[ sorted_dirac_nums[ i ] ] = beg_x_density + mul * sorted_dirac_positions[ i ];
    return res;
}

DTP void UTP::set_density_contrast( TF max_ratio ) {    
    if ( verbosity > 0 )
        P( max_ratio );
    
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

DTP TF UTP::normalized_error() const {
    using namespace std;

    TF res = 0;
    for_each_normalized_cell_mass( [&]( PI i, TF v, bool bad_cell ) {
        res += pow( log( v / sorted_dirac_masses[ i ] ), 2 );
    } );

    return sqrt( res );
}

DTP void UTP::plot( Str filename ) const {
    std::ofstream fs( filename );

    Vec<TF> xs;
    Vec<TF> ys;
    xs << beg_x_density;
    ys << original_density_values[ 0 ];
    for( PI i = 1; i < original_density_values.size(); ++i ) {
        const TF x = beg_x_density + ( end_x_density - beg_x_density ) * ( i + 0 ) / original_density_values.size();
        ys << original_density_values[ i - 1 ] << original_density_values[ i - 0 ];
        xs << x << x;
    }
    ys << original_density_values.back();
    xs << end_x_density;


    Vec<TF> bx = cell_barycenters();
    Vec<TF> by = bx * 0 - 0.1;

    Vec<TF> dx = dirac_positions();
    Vec<TF> dy = dx * 0 - 0.2;

    fs << "from matplotlib import pyplot\n";
    fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
    fs << "pyplot.plot( " << to_string( bx ) << ", " << to_string( by ) << ", '+' )\n";
    fs << "pyplot.plot( " << to_string( dx ) << ", " << to_string( dy ) << ", '+' )\n";
    fs << "pyplot.show()\n";
}

DTP Vec<TF> UTP::newton_dir() const {
    // using TM = Eigen::Matrix<TF,Eigen::Dynamic,Eigen::Dynamic>;
    // using TV = Eigen::Matrix<TF,Eigen::Dynamic,1>;

    // TM M( nb_diracs(), nb_diracs() );    
    // TV V( nb_diracs() );
    // M.fill( 0 );

    // TF mp = 0;
    // for_each_normalized_system_item( [&]( PI i, TF m0, TF m1, TF v, bool bad_cell ) {
    //     M.coeffRef( i, i + 0 ) = m0 / v;
    //     if ( m1 )
    //         M.coeffRef( i, i + 1 ) = m1 / v;
    //     if ( mp )
    //         M.coeffRef( i, i - 1 ) = mp / v;
    //     mp = m1;

    //     V.coeffRef( i ) = log( sorted_dirac_masses[ i ] / v );    
    //     // M.coeffRef( i + 0, i + 0 ) = m0;
    //     // if ( m1 )
    //     //     M.coeffRef( i + 0, i + 1 ) = m1;
    //     // if ( mp )
    //     //     M.coeffRef( i + 0, i - 1 ) = mp;
    //     // mp = m1;

    //     // V.coeffRef( i + 0 ) = sorted_dirac_masses[ i ] - v;
    // } );

    // // std::cout << M << std::endl;
    // // std::cout << V << std::endl;
    // Eigen::PartialPivLU<TM> lu( M );
    // return lu.solve( V );

    // int n = nb_diracs();
    // TM upper( n, n ); upper.fill( 0 );
    // TM lower( n, n ); lower.fill( 0 );
    // for (int i = 0; i < n; i++) {
    //     // diagonal
    //     TF sum = 0;
    //     if ( i )
    //         sum += lower(i,i-1) * upper(i-1,i);
    //     upper.coeffRef(i,i) = M(i,i) - sum;

    //     // Upper Triangular
    //     if ( i + 1 < n )
    //         upper.coeffRef(i,i+1) = M(i,i+1);

    //     // Lower Triangular
    //     lower.coeffRef(i,i) = 1;
    //     if ( i + 1 < n )
    //         lower.coeffRef(i+1,i) = M(i+1,i) / upper(i,i);
    // }
    
    // std::cout << upper << std::endl;
    // std::cout << lower << std::endl;

    /*
        upper.coeffRef(i,i) = M(i,i) - lower(i,i-1) * upper(i-1,i);
        upper.coeffRef(i,i+1) = M(i,i+1);
        lower.coeffRef(i+1,i) = M(i+1,i) / upper(i,i);
    */

    // direction at 0
    TF prev_m1 = 0;
    TF prev_x = 0;
    TF prev_u = 0;
    TF prev_d = 0;
    Vec<TF> us( FromSize(), nb_diracs() );
    Vec<TF> ds( FromSize(), nb_diracs() );
    Vec<TF> x0( FromSize(), nb_diracs() );
    for_each_normalized_system_item( [&]( PI index, TF m0, TF m1, TF v, bool bad_cell ) {
        const TF y = log( sorted_dirac_masses[ index ] / v );
        const TF l = prev_m1 ? prev_m1 / ( v * prev_d ) : 0;
        const TF d = m0 / v - prev_u * l;
        const TF u = m1 / v;
        ds[ index ] = d;
        us[ index ] = u;

        const TF x = y - l * prev_x;
        x0[ index ] = x;
        
        prev_m1 = m1;
        prev_x = x;
        prev_u = u;
        prev_d = d;
    } );

    x0.back() /= ds.back();
    for( PI i = nb_diracs() - 1; i--; )
        x0[ i ] = ( x0[ i ] - us[ i ] * x0[ i + 1 ] ) / ds[ i ];
    
    // first derivative
    prev_m1 = 0;
    prev_x = 0;
    prev_u = 0;
    prev_d = 0;
    for_each_normalized_system_item( [&]( PI index, TF m0, TF m1, TF v, bool bad_cell ) {
        // const TF y = - log( ( v + x0[ index ] * a ) / sorted_dirac_masses[ index ] );
        const TF y = - x0[ index ] * sorted_dirac_masses[ index ] / v;

        const TF l = prev_m1 ? prev_m1 / ( v * prev_d ) : 0;
        const TF d = m0 / v - prev_u * l;
        const TF u = m1 / v;
        ds[ index ] = d;
        us[ index ] = u;

        const TF x = y - l * prev_x;
        x0[ index ] = x;
        
        prev_m1 = m1;
        prev_x = x;
        prev_u = u;
        prev_d = d;
    } );

    x0.back() /= ds.back();
    for( PI i = nb_diracs() - 1; i--; )
        x0[ i ] = ( x0[ i ] - us[ i ] * x0[ i + 1 ] ) / ds[ i ];

    return x0;
}

DTP void UTP::update_weights() {
    TF last_error = normalized_error();
    for( PI num_iter = 0; ; ++num_iter ) {
        if ( num_iter == 15 )
            throw std::runtime_error( "too many iterations" );
        
        Vec<TF> old_dirac_weights = sorted_dirac_weights;
        Vec<TF> dir = newton_dir();
        
        TF end_a = 1, end_error = 0;
        for( ;; end_a /= 2 ) {
            // if ( end_a == 0 )
            //     throw std::runtime_error( "bad direction" );
            sorted_dirac_weights = old_dirac_weights + end_a * dir;
            end_error = normalized_error();
            if ( end_error < last_error )
                break;
        }
    
        Vec<TF> as, errors;
        for( TF a : Vec<TF>::linspace( 0, end_a, 1000, false ) ) {
            sorted_dirac_weights = old_dirac_weights + a * dir;
             TF error = normalized_error();
            if ( std::isnan( error ) )
                throw std::runtime_error( "unexpected nan" );
            errors << error;
            as << a;
        }
        errors << end_error;
        as << end_a;
        // glot_vec( xs, errors );
    
        PI best_i = argmin( errors );
        TF best_a = as[ best_i ];
        if ( best_a == 0 )
            throw std::runtime_error( "weirdox" );
        if ( verbosity > 0 )
            P( best_a, errors[ best_i ] );
        
        sorted_dirac_weights = old_dirac_weights + best_a * dir;
        last_error = errors[ best_i ];

        if ( errors[ best_i ] < target_l2_error )
            break;
    }
}

DTP void UTP::solve() {
    update_weights();

    while ( current_contrast_ratio > target_contrast_ratio ) {
        set_density_contrast( std::max( current_contrast_ratio / 2, target_contrast_ratio ) );
        update_weights();
    }
}

#undef DTP
#undef UTP


} // namespace usdot
