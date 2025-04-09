#pragma once

// #include <tl/support/operators/self_max.h>
// #include <tl/support/operators/norm_2.h>
// #include <tl/support/operators/argmin.h>
// #include <tl/support/operators/max.h>
// #include <tl/support/operators/sum.h>
#include <tl/support/ASSERT.h>
#include <tl/support/P.h>
#include <algorithm>
#include <stdexcept>
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

    // values that have to be copied
    global_mass_ratio = input.global_mass_ratio;
    target_l2_error = input.target_l2_error;
    throw_if_error = input.throw_if_error;
    end_x_density = input.end_x_density;
    beg_x_density = input.beg_x_density;
    multithread = input.multithread;
    verbosity = input.verbosity;

    // values that have to be moved
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

    // densities
    TF cut_ratio = 1e-6, target_mass = sum_of_dirac_masses / global_mass_ratio;
    densities.push_back( new ControledDensity<TF>( original_density_values, 0, target_mass ) );
    for( TF ratio : input.filter_values )
        densities.push_back( new ControledDensity<TF>( densities[ 0 ]->values, ratio ) );
}


DTP PI UTP::nb_cells() const {
    return sorted_dirac_nums.size();
}

DTP void UTP::for_each_normalized_cell( auto &&get_weight, auto &&func ) const {
    using namespace std;

    auto i0 = numeric_limits<TF>::lowest();
    auto d0 = sorted_dirac_positions[ 0 ];
    auto w0 = get_weight( 0 );
    for( PI n = 0; n < nb_cells(); ++n ) {
        const auto d1 = n + 1 < nb_cells() ? sorted_dirac_positions[ n + 1 ] : numeric_limits<TF>::max();
        const auto w1 = n + 1 < nb_cells() ? get_weight( n + 1 ) : 0;
        if ( w0 < 0 ) {
            func( n, 0, 0, true );
            continue;
        }

        const auto i1 = ( d0 + d1 + ( w0 - w1 ) / ( d1 - d0 ) ) / 2;
        const auto r0 = sqrt( w0 );
        const auto bl = d0 - r0;
        const auto br = d0 + r0;

        if ( i0 <= i1 ) {
            const auto il = max( i0, bl );
            const auto ir = min( i1, br );
            func( n, il, ir, il > ir );
        } else {
            const auto il = max( i1, bl );
            const auto ir = min( i0, br );
            func( n, ir, il, il > ir );
        }

        i0 = i1;
        d0 = d1;
        w0 = w1;
    }
}

DTP void UTP::for_each_normalized_cell( auto &&func ) const {
    for_each_normalized_cell( [&]( PI n ) { return sorted_dirac_weights[ n ]; }, func );
}

DTP void UTP::for_each_normalized_system_item( auto &&func ) const {
    using namespace std;

    Vec<TF> mds( FromSize(), densities.size() );
    Vec<TF> mes( FromSize(), densities.size() );

    auto c0 = numeric_limits<TF>::lowest();
    auto d0 = numeric_limits<TF>::max();
    auto p0 = sorted_dirac_positions[ 0 ];
    auto w0 = sorted_dirac_weights[ 0 ];
    for( PI num_dirac = 0; num_dirac < nb_cells(); ++num_dirac ) {
        const auto p1 = num_dirac + 1 < nb_cells() ? sorted_dirac_positions[ num_dirac + 1 ] : numeric_limits<TF>::max();
        const auto w1 = num_dirac + 1 < nb_cells() ? sorted_dirac_weights[ num_dirac + 1 ] : 0;
        const auto d1 = p1 - p0;
        if ( w0 < 0 ) {
            func( num_dirac, mds, mes, 0, true );
            continue;
        }

        const auto c1 = ( p0 + p1 + ( w0 - w1 ) / d1 ) / 2;
        const auto r0 = sqrt( w0 );
        const auto b0 = p0 - r0;
        const auto b1 = p0 + r0;

        if ( c0 > c1 ) { // "wrong" ordering
            if ( c0 < b0 || c1 > b1 ) { // void cell
                func( num_dirac, mds, mes, 0, true );
            } else { // ball cut on the left
                if ( b0 > c1 ) { // ball cut on the right
                    const TF rsq = 0.5 / sqrt( w0 );
                    if ( b1 < c0 ) {
                        for( PI n = 0; n < densities.size(); ++n ) {
                            mds[ n ] = rsq * ( densities[ n ]->value( b0 ) + densities[ n ]->value( b1 ) );
                            mes[ n ] = 0;                        
                        }                        
                        func( num_dirac, mds, mes, densities[ 0 ]->integral( b1, b0 ), false );
                    } else { // interface on the right
                        const TF rd0 = 0.5 / d0;
                        for( PI n = 0; n < densities.size(); ++n ) {
                            mds[ n ] = rd0 * densities[ n ]->value( c0 ) - rsq * densities[ n ]->value( b0 );
                            mes[ n ] = 0;                        
                        }
                        func( num_dirac, mds, mes, densities[ 0 ]->integral( c0, b0 ), false );
                    }
                } else { // interface on the left
                    if ( b1 < c0 ) { // ball cut on the right
                        const TF rsq = 0.5 / sqrt( w0 );
                        const TF rd1 = 0.5 / d1;
                        for( PI n = 0; n < densities.size(); ++n ) {
                            const TF cl = rd1 * densities[ n ]->value( c1 );
                            mds[ n ] = cl - rsq * densities[ n ]->value( b1 );
                            mes[ n ] = - cl;                        
                        }                        
                        func( num_dirac, mds, mes, densities[ 0 ]->integral( b1, c1 ), false );
                    } else { // interface on the right
                        const TF rd0 = 0.5 / d0;
                        const TF rd1 = 0.5 / d1;
                        for( PI n = 0; n < densities.size(); ++n ) {
                            const TF cl = rd1 * densities[ n ]->value( c1 );
                            const TF cr = rd0 * densities[ n ]->value( c0 );
                            mds[ n ] = - ( cl + cr );
                            mes[ n ] = cr;
                        }
                        func( num_dirac, mds, mes, densities[ 0 ]->integral( c0, c1 ), false );
                    }
                }
            }
        } else {
            if ( c1 < b0 || c0 > b1 ) { // void cell
                func( num_dirac, mds, mes, 0, true );
            } else { // non empty cell
                if ( b0 > c0 ) { // ball cut on the left
                    const TF rsq = 0.5 / sqrt( w0 );
                    if ( b1 < c1 ) { // ball cut on the right
                        for( PI n = 0; n < densities.size(); ++n ) {
                            mds[ n ] = rsq * ( densities[ n ]->value( b0 ) + densities[ n ]->value( b1 ) );
                            mes[ n ] = 0;
                        }
                        func( num_dirac, mds, mes, densities[ 0 ]->integral( b0, b1 ), false );
                    } else { // interface on the right
                        const TF rd1 = 0.5 / d1;
                        for( PI n = 0; n < densities.size(); ++n ) {
                            const TF cr = rd1 * densities[ n ]->value( c1 );
                            mds[ n ] = rsq * densities[ n ]->value( b0 ) + cr;
                            mes[ n ] = - cr;
                        }
                        func( num_dirac, mds, mes, densities[ 0 ]->integral( b0, c1 ), false );
                    }
                } else { // interface on the left
                    const TF rd0 = 0.5 / d0;
                    if ( b1 < c1 ) { // ball cut on the right
                        const TF rsq = 0.5 / sqrt( w0 );
                        for( PI n = 0; n < densities.size(); ++n ) {
                            mds[ n ] = rd0 * densities[ n ]->value( c0 ) + rsq * densities[ n ]->value( b1 );
                            mes[ n ] = 0;
                        }
                        func( num_dirac, mds, mes, densities[ 0 ]->integral( c0, b1 ), false );
                    } else {
                        const TF rd1 = 0.5 / d1;
                        for( PI n = 0; n < densities.size(); ++n ) {
                            const TF cl = rd0 * densities[ n ]->value( c0 );
                            const TF cr = rd1 * densities[ n ]->value( c1 );
                            mds[ n ] = cl + cr;
                            mes[ n ] = - cr;
                        }
                        func( num_dirac, mds, mes, densities[ 0 ]->integral( c0, c1 ), false );
                    }
                }
            }
        }
    
        c0 = c1;
        d0 = d1;
        p0 = p1;
        w0 = w1;
    }
}

DTP Vec<Vec<TF,2>> UTP::normalized_cell_boundaries() const {
    Vec<Vec<TF,2>> res( FromSize(), nb_cells() );
    for_each_normalized_cell( [&]( PI n ) { return sorted_dirac_weights[ n ]; }, [&]( PI num_dirac, TF i0, TF i1, bool ) {
        res[ num_dirac ] = { i0, i1 };
    } );
    return res;
}

DTP UTP::TV UTP::normalized_cell_barycenters() const {
    TV res( FromSize(), nb_cells() );
    for_each_normalized_cell( [&]( PI num_dirac, TF cl, TF cr, bool ) {
        const TF it = densities[ 0 ]->integral( cl, cr );
        res[ num_dirac ] = it ? densities[ 0 ]->x_integral( cl, cr ) / it : 0;
    } );
    return res;
}

DTP UTP::TV UTP::cell_barycenters() const {
    const TF mul = ( end_x_density - beg_x_density ) / ( original_density_values.size() - 1 );
    TV bar = normalized_cell_barycenters();

    TV res( FromSize(), nb_cells() );
    for( PI i = 0; i < nb_cells(); ++i )
        res[ sorted_dirac_nums[ i ] ] = beg_x_density + mul * bar[ i ];
    return res;
}

DTP TF UTP::density_value( TF pos, PI n ) const {
    const TF mul = ( original_density_values.size() - 1 ) / ( end_x_density - beg_x_density );
    return densities[ n ]->value( mul * ( pos - beg_x_density ) );
}

DTP UTP::TV UTP::dirac_positions() const {
    const TF mul = ( end_x_density - beg_x_density ) / ( original_density_values.size() - 1 );
    TV res( FromSize(), nb_cells() );
    for( PI i = 0; i < nb_cells(); ++i )
        res[ sorted_dirac_nums[ i ] ] = beg_x_density + mul * sorted_dirac_positions[ i ];
    return res;
}

DTP TF UTP::error() const {
    using namespace std;

    TF res = 0;
    bool has_bad_cell = false;
    for_each_normalized_cell( [&]( PI n ) { return sorted_dirac_weights[ n ]; }, [&]( PI i, TF cl, TF cr, bool bad_cell ) {
        res += pow( sorted_dirac_masses[ i ] - densities[ 0 ]->integral( cl, cr ), 2 );
        if ( cr - cl < 1e-3 )
            has_bad_cell = true;   
        has_bad_cell |= bad_cell;
    } );

    return has_bad_cell ? numeric_limits<TF>::max() : sqrt( res );
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

DTP Vec<Vec<TF>> UTP::newton_dirs() const {
    bool has_bad_cell = false;
    Vec<TV> x( FromSizeAndItemValue(), densities.size(), FromSize(), nb_cells() );
    Vec<TV> e( FromSizeAndItemValue(), densities.size(), FromSize(), nb_cells() );
    Vec<TF> prev_dll( FromSizeAndItemValue(), densities.size(), 0 );
    Vec<TF> prev_lv( FromSizeAndItemValue(), densities.size(), 0 );
    Vec<TF> bad_dir( FromSizeAndItemValue(), densities.size(), 0 );
    TF error = 0;
    for_each_normalized_system_item( [&]( PI index, const auto &m0s, const auto &m1s, TF v, bool bad_cell ) {
        if ( bad_cell ) {
            has_bad_cell = true;
            return;
        }

        v -= sorted_dirac_masses[ index ];
        error += pow( v, 2 );

        for( PI nd = 0; nd < densities.size(); ++nd ) {            
            TF d = m0s[ nd ] - prev_dll[ nd ];
            if ( d == 0 ) {
                bad_dir[ nd ] = true;
                d = 1;
            }
            const TF l = m1s[ nd ] / d;
            e[ nd ][ index ] = l;

            prev_dll[ nd ] = d * l * l;

            const TF v0 = v - prev_lv[ nd ];
            x[ nd ][ index ] = v0 / d;
            prev_lv[ nd ] = l * v0;
        }
    } );

    for( PI nd = 0; nd < densities.size(); ++nd ) {
        if ( bad_dir[ nd ] ) {
            for( PI i = nb_cells() - 1; i--; )
                x[ nd ][ i ] = 0;
            continue;
        }
        for( PI i = nb_cells() - 1; i--; )
            x[ nd ][ i ] -= e[ nd ][ i ] * x[ nd ][ i + 1 ];
    }

    return x; // { x, has_bad_cell ? std::numeric_limits<TF>::max() : error };
}

DTP Vec<TF> UTP::cell_masses() const {
    TV res( FromSize(), nb_cells() );
    for_each_normalized_cell( [&]( PI num_dirac, TF i0, TF i1, bool wrong_cell ) {
        res[ num_dirac ] = densities[ 0 ]->integral( i0, i1 );
    } );
    return res;
}

DTP int UTP::update_weights() {
    using namespace std;

    //
    Vec<Vec<TF>> dirs = newton_dirs();
    // P( dirs );
    
    auto old_weights = sorted_dirac_weights;
    TF best_error = error();
    PI best_i = -1;
    for( TF relax = 1; ; relax *= 0.5 ) {
        if ( relax < 1e-9 )
            throw std::runtime_error( "bad relaxation" );

        for( PI i = 0; i < dirs.size(); ++i ) {
            sorted_dirac_weights = old_weights - relax * dirs[ i ];
            const TF e = error();
            // P( relax, i, e );
            if ( e < best_error ) {
                best_error = e;
                best_i = i;
            }        
        }

        if ( best_i != -1 ) {
            // if ( best_i != dirs.size() - 1 )
            sorted_dirac_weights = old_weights - relax * dirs[ best_i ];
            P( best_i, relax, best_error );
            break;
        }
    }
    
    return 0;
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

DTP void UTP::solve() {
    if ( nb_cells() == 0 )
        return;

    // solve
    Vec<Vec<Vec<TF,2>>> bnds;
    bnds << normalized_cell_boundaries(); 
    for( PI num_iter = 0; num_iter < 17; ++num_iter ) {
        update_weights();
        bnds << normalized_cell_boundaries(); 
    }

    //
    plot_bnds_evolution( bnds );
}

#undef DTP
#undef UTP


} // namespace usdot
