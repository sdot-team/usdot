#pragma once

// #include <tl/support/operators/self_max.h>
// #include <tl/support/operators/norm_2.h>
#include <tl/support/operators/argmin.h>
// #include <tl/support/operators/max.h>
// #include <tl/support/operators/sum.h>
#include <algorithm>
#include <stdexcept>
#include <tl/support/ASSERT.h>
#include <tl/support/P.h>
#include <numeric>
#include <fstream>
#include <limits>

#include "ConvGridSolver.h"
#include "glot.h"
#include "partial1D/Extrapolation.h"

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

    // density_gaussian_width
    initialize_filter_value( input.starting_filter_value );

    // weights
    sorted_dirac_weights.resize( nb_diracs );
    for( PI i = 0; i < nb_diracs; ++i ) {
        sorted_dirac_weights[ i ] = sorted_dirac_nums[ i ] < input.dirac_weights.size() ? 
            input.dirac_weights[ sorted_dirac_nums[ i ] ] * pow( mul_dirac, 2 ) : 
            pow( 0.5 * sorted_dirac_masses[ i ] / density->value( sorted_dirac_positions[ i ] ), 2 );
    }
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
    TV res( nb_diracs() );
    for_each_normalized_cell( [&]( TF dirac_position, TF dirac_weight, PI num_dirac, TF, TF, TF c0, TF c1 ) {
        const TF rd = sqrt( dirac_weight );
        if ( c0 > c1 ) {
            const TF x0 = max( dirac_position - rd, c1 );
            const TF x1 = min( dirac_position + rd, c0 );
            res[ num_dirac ] = - density->integral( x0, x1 );
        } else {
            const TF x0 = max( dirac_position - rd, c0 );
            const TF x1 = min( dirac_position + rd, c1 );
            res[ num_dirac ] = density->integral( x0, x1 );
        }        
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
    TV by( FromSizeAndItemValue(), bx.size(), -0.1 );

    TV dx = dirac_positions();
    TV dy( FromSizeAndItemValue(), dx.size(), -0.2 );

    fs << "from matplotlib import pyplot\n";
    fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
    fs << "pyplot.plot( " << to_string( bx ) << ", " << to_string( by ) << ", '+' )\n";
    fs << "pyplot.plot( " << to_string( dx ) << ", " << to_string( dy ) << ", '+' )\n";
    fs << "pyplot.show()\n";
}

// DTP std::function<typename UTP::TV( TF )> UTP::newton_dirs() {
//     using TM = Eigen::Matrix<TF,Eigen::Dynamic,Eigen::Dynamic>;
//     using EV = Eigen::Matrix<TF,Eigen::Dynamic,1>;
//     using namespace std;

//     Vec<TF> V = Vec<TF>( sorted_dirac_masses ) - Vec<TF>( normalized_cell_masses() );
//     TM M( nb_diracs(), nb_diracs() );
//     const TF eps = 1e-8;
//     for( PI n = 0; n < nb_diracs(); ++n ) {
//         // en 0, w = u^2 -> dw = 2 * sqrt( w ) * du
//         TF &ref_weight = sorted_dirac_weights[ n ];
//         TF old_weight = ref_weight;
//         if ( n == 0 || n + 1 == nb_diracs() )
//             ref_weight += 2 * sqrt( old_weight ) * eps;
//         else
//             ref_weight += eps;

//         Vec<TF> A = Vec<TF>( sorted_dirac_masses ) - Vec<TF>( normalized_cell_masses() );
//         for( PI m = 0; m < nb_diracs(); ++m )
//             M( m, n ) = ( V[ m ] - A[ m ] ) / eps;

//         ref_weight = old_weight;
//     }

//     EV Y( nb_diracs() );
//     for( PI m = 0; m < nb_diracs(); ++m )
//         Y[ m ] = V[ m ];

//     Eigen::FullPivLU<TM> lu( M );
//     EV X = lu.solve( Y );

//     return [X,sorted_dirac_weights=this->sorted_dirac_weights,n=nb_diracs()]( TF a ) {
//         TV res( n );
//         for( PI i = 0; i < n; ++i ) {
//             // en 0, w = u^2 -> dw = 2 * sqrt( w ) * du
//             if ( i == 0 || i + 1 == n )
//                 res[ i ] = pow( sqrt( sorted_dirac_weights[ i ] ) + a * X[ i ], 2 );
//             else
//                 res[ i ] = sorted_dirac_weights[ i ] + a * X[ i ];
//         }
//         return res;
//     };

    
//     // TV d = newton_dir().first;

//     // for( PI i = 0; i < nb_diracs(); ++i )
//     //     sorted_dirac_weights[ i ] += 1e-6 * d[ i ];
//     // TV e = newton_dir().first;

//     // for( PI i = 0; i < nb_diracs(); ++i )
//     //    sorted_dirac_weights[ i ] += 1e-6 * d[ i ];
//     // TV f = newton_dir().first;

//     // for( PI i = 0; i < nb_diracs(); ++i )
//     //     sorted_dirac_weights[ i ] -= 2e-6 * d[ i ];

//     // for( PI i = 0; i < nb_diracs(); ++i )
//     //    f[ i ] = ( d[ i ] + f[ i ] - 2 * e[ i ] ) / 40e-12;

//     // for( PI i = 0; i < nb_diracs(); ++i )
//     //    e[ i ] = ( e[ i ] - d[ i ] ) / 1e-6;


//     // return [d,e,f,sorted_dirac_weights=this->sorted_dirac_weights,n=nb_diracs()]( TF a ) {
//     //     TV res( n );
//     //     for( PI i = 0; i < n; ++i )
//     //         res[ i ] = sorted_dirac_weights[ i ]
//     //                  + a * d[ i ]
//     //                  + a * a / 3.4 * e[ i ]
//     //                  + a * a * a / 3 * f[ i ] * 0
//     //                  ;
//     //     return res;
//     // };
// }


DTP void UTP::for_each_normalized_system_item( auto &&bad_cell, auto &&bb, auto &&bi, auto &&ib, auto &&ii ) const {
    using namespace std;

    for_each_normalized_cell( [&]( TF dirac_position, TF dirac_weight, PI num_dirac, TF d0, TF d1, TF c0, TF c1 ) {
        const TF v = sorted_dirac_masses[ num_dirac ];
        if ( dirac_weight <= 0 )
            return bad_cell( num_dirac );

        const TF rd = sqrt( dirac_weight );
        const TF b0 = dirac_position - rd;
        const TF b1 = dirac_position + rd;

        // "wrong" bounds ordering
        if ( c0 > c1 )
            return bad_cell( num_dirac );

        // if ( c0 > c1 ) {
        //     // void cell
        //     if ( c0 < b0 || c1 > b1 )
        //         return bad_cell( num_dirac );

        //     // ball cut on the left
        //     if ( b0 > c1 ) {
        //         // ball cut on the right
        //         if ( b1 < c0 )
        //             return bb( num_dirac );

        //         // interface on the right
        //         return bc( num_dirac, 
        //             ( density->value( c0 ) / d0 - density->value( b0 ) / sqrt( dirac_weight ) ) / 2,
        //             0,
        //             - density->integral( b0, c0 ),
        //             false
        //         );
        //     }
            
        //     // interface on the left
        //     if ( b1 < c0 ) { // ball cut on the right
        //         const TF cl = density->value( c1 ) / ( 2 * d1 );
        //         return func( num_dirac, 
        //             cl - density->value( b1 ) / ( 2 * sqrt( dirac_weight ) ),
        //             - cl,
        //             - density->integral( c1, b1 ),
        //             false
        //         );
        //     }
            
        //     // interface on the right
        //     const TF cl = density->value( c1 ) / ( 2 * d1 );
        //     const TF cr = density->value( c0 ) / ( 2 * d0 );
        //     return func( num_dirac, 
        //         - ( cl + cr ),
        //         cr,
        //         - density->integral( c1, c0 ),
        //         false
        //     );
        // }

        // void cell
        if ( c1 < b0 || c0 > b1 )
            return bad_cell( num_dirac );

        if ( b0 > c0 ) { // ball cut on the left
            // ball cut on the right
            if ( b1 < c1 )
                return bb( num_dirac, dirac_position, rd );

            // interface on the right
            return bi( num_dirac, - density->value( c1 ) / ( 2 * d1 ) );
        } 
        
        // interface on the left
        if ( b1 < c1 ) // ball cut on the right
            return ib( num_dirac );
        
        const TF cl = density->value( c0 ) / ( 2 * d0 );
        const TF cr = density->value( c1 ) / ( 2 * d1 );
        const TF in = density->integral( c0, c1 );
        return ii( num_dirac, cl + cr, - cr, in );
    } );
}

DTP typename UTP::TV UTP::newton_dir() const {
    TV res( FromSize(), nb_diracs() );
    TV x0( FromSize(), nb_diracs() );
    TV x1( FromSize(), nb_diracs() );
    TV e( FromSize(), nb_diracs() );
    bool has_bad_cell = false;
    TF prev_dll = 0;
    TF prev_lv = 0;
    TF error = 0;
    for_each_normalized_system_item( 
        // bad cell
        [&]( PI index ) {
            has_bad_cell = true;
        },

        // bb
        [&]( PI index, TF center, TF prev_rad ) {
            const TF rad = density->radius_for( sorted_dirac_masses[ index ], 1e-4, center, prev_rad * 2 );
            res[ index ] = pow( rad, 2 ) - sorted_dirac_weights[ index ];
            P( pow( rad, 2 ), sorted_dirac_weights[ index ] );
        },

        // bc
        [&]( PI index, TF m1 ) {
            TODO;
        },

        // cb
        [&]( PI index ) {
            TODO;
        },

        // cc
        [&]( PI index, TF m0, TF m1, TF v ) {
            TODO;

            v = sorted_dirac_masses[ index ] - v;

            TF d = m0 - prev_dll;
            if ( d <= 0 )
                d = 1;

            const TF l = m1 / d;
            e[ index ] = l;

            prev_dll = d * l * l;

            const TF v0 = v - prev_lv;
            x0[ index ] = v0 / d;
            prev_lv = l * v0;
        }
    );

    // for( PI i = nb_diracs() - 1; i--; ) {
    //     x0[ i ] -= e[ i ] * x0[ i + 1 ];
    // }

    return res;
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

DTP int UTP::update_weights() {
    Vec<Vec<Vec<TF,2>>> bnds;
    for( PI num_iter = 0; ; ++num_iter ) {
        bnds << normalized_cell_boundaries(); 
        if ( num_iter == 500 )
            throw std::runtime_error( "too many iterations" );

        // line search
        Vec<TF> base_dirac_weights = sorted_dirac_weights;
        TF base_error = normalized_error();
        Vec<TF> dir = newton_dir();

        Vec<TF> as = Vec<TF>::linspace( 0, 1.5, 1000, false );
        Vec<TF> errors;
        for( TF a : as ) {
            sorted_dirac_weights = base_dirac_weights + a * dir;
            errors << normalized_error();            
        }

        PI best_i = argmin( errors );
        TF best_a = as[ best_i ];

        sorted_dirac_weights = base_dirac_weights + best_a * dir;

        if ( errors[ best_i ] < target_l2_error ) {
            plot_bnds_evolution( bnds );
            P( num_iter );
            return 0;
        }

        if ( best_a == 0 )
            throw std::runtime_error( "bad direction" );
        
        P( best_a, errors[ best_i ] );

    }

    return 0;
}

DTP void UTP::solve() {
    if ( nb_diracs() == 0 )
        return;

    // solve
    initialize_filter_value( pow( 0.5, 1.0 / original_density_values.size() ) );
    int not_solved = update_weights();
    if ( not_solved )
        throw std::runtime_error( "not solved with the first filter value" );

    // // increase the filter width to get a first solution
    // std::cout << "increase the filter width to get a first solution" << std::endl;
    // const TF initial_filter_value = current_filter_value;
    // while ( true ) {
    //     if ( update_weights() == 0 )
    //         break;
    //     std::fill( sorted_dirac_weights.begin(), sorted_dirac_weights.end(), init_w );
    //     initialize_filter_value( 0.5 * ( 1 + current_filter_value ) );
    //     P( current_filter_value );

    //     if ( current_filter_value > 0.95 ) {
    //         glot( ::Vec<TF>::linspace( -2, 5, 1000 ), 
    //             [&]( TF x ) { return density_value( x ); }
    //         );
    //         throw std::runtime_error( "unable to start" );
    //         return;
    //     }
    // }

    // // decrease the filter width to go to the target value
    // std::cout << "decrease the filter width to go to the target value" << std::endl;
    // while ( current_filter_value > target_filter_value ) {
    //     const TF prev_filter_value = current_filter_value;
    //     Vec old_weights = sorted_dirac_weights;
    //     for( TF next_filter_value = target_filter_value; ; next_filter_value = 0.5 * ( prev_filter_value + next_filter_value ) ) {
    //         go_to_filter_value( next_filter_value );
    //         P( current_filter_value );
    //         if ( update_weights() == 0 )
    //             break;
    //         sorted_dirac_weights = old_weights;
    //     }
    // }
}

#undef DTP
#undef UTP


} // namespace usdot
