#pragma once

#include "utility/dichotomy.h"
#include "UsdotSystem.h"
#include <stdexcept>
#include <iostream>
#include <numeric>
#include <fstream>
#include <limits>
#include <cmath>

namespace usdot {

#define DTP template<class TF,class Density>
#define UTP UsdotSystem<TF,Density>

DTP UTP::UsdotSystem( Density *density, const VF &dirac_positions, TF global_mass_ratio, TF sep_equ_coords ) : density( density ), global_mass_ratio( global_mass_ratio ) {
    initialize_dirac_positions( dirac_positions, sep_equ_coords );
    initialize_dirac_masses( {}, global_mass_ratio );
    initialize_dirac_weights();
}

DTP void UTP::initialize_dirac_positions( const VF &dirac_positions, TF sep_equ_coords ) {
    using namespace std;

    // sorted_dirac_nums
    VI sorted_dirac_nums( dirac_positions.size() );
    iota( sorted_dirac_nums.begin(), sorted_dirac_nums.end(), PI( 0 ) );
    sort( sorted_dirac_nums.begin(), sorted_dirac_nums.end(), [&]( PI a, PI b ) {
        return dirac_positions[ a ] < dirac_positions[ b ];
    } );

    // sorted_dirac_positions
    sorted_dirac_num_offsets.reserve( dirac_positions.size() );
    sorted_dirac_num_values.reserve( dirac_positions.size() );
    sorted_dirac_positions.reserve( dirac_positions.size() );
    for( PI i = 0; i < sorted_dirac_nums.size(); ++i ) {
        const PI s = sorted_dirac_num_values.size();
        const PI n = sorted_dirac_nums[ i ];
        const TF p = dirac_positions[ n ];
        if ( s && sorted_dirac_positions.back() + sep_equ_coords >= p ) {
            sorted_dirac_num_values.push_back( n );
            continue;
        }

        sorted_dirac_num_offsets.push_back( s );
        sorted_dirac_num_values.push_back( n );
        sorted_dirac_positions.push_back( p );
    }
    sorted_dirac_num_offsets.push_back( sorted_dirac_num_values.size() );
}

DTP void UTP::initialize_dirac_weights() {
    using namespace std;

    // init weights as if cell were not overlapping. We use the same mass to ensure that no cell can "eat" another one.
    sorted_dirac_weights.resize( nb_sorted_diracs() );
    const TF mass_ref = global_mass_ratio / nb_sorted_diracs();
    for( PI i = 0; i < nb_sorted_diracs(); ++i ) {
        const TF dx = sorted_dirac_positions[ i ];
        auto err = [&]( TF r ) {
            return density->integral( dx - r, dx + r ) - mass_ref;
        };
        const TF r = dichotomy_growing_from_zero( err, mass_ref * 1e-2, TF( 1 ) );
        sorted_dirac_weights[ i ] = pow( r, 2 );
    }    

    //
    beg_unknowns = 0;
    TF last_b = density->min_x();
    VF bnds;
    for( ; beg_unknowns < nb_sorted_diracs(); ++beg_unknowns ) {
        const TF dx = sorted_dirac_positions[ beg_unknowns ];
        const TF we = sorted_dirac_weights[ beg_unknowns ];
        const TF dm = sorted_dirac_masses[ beg_unknowns ];
        const TF b0 = dx - sqrt( we );
        if ( b0 > last_b )
            break;

        auto err = [&]( TF add_b ) {
            // P( last_b, add_b, density->integral( last_b, last_b + add_b ) );
            return density->integral( last_b, last_b + add_b ) - dm;
        };

        TF add_b = dichotomy_growing_from_zero( err, target_max_error_ratio * dm, TF( 1 ) );
        last_b += add_b;
        bnds.push_back( last_b );
    }

    if ( beg_unknowns ) {
        const TF dx = sorted_dirac_positions[ beg_unknowns - 1 ];
        TF w = pow( bnds.back() - dx, 2 );
        for( PI i = beg_unknowns; --i; ) {
            const TF dx0 = sorted_dirac_positions[ i - 1 ];
            const TF dx1 = sorted_dirac_positions[ i - 0 ];
            const TF bnd = bnds[ i - 1 ];

            sorted_dirac_weights[ i ] = w;

            // ( w1 - w0 ) / ( d1 - d0 ) = d1 + d0 - 2 * bnd;
            w -= ( dx1 + dx0 - 2 * bnd ) * ( dx1 - dx0 );
        }

        sorted_dirac_weights[ 0 ] = w;            
    }

    //
    end_unknowns = nb_sorted_diracs();
    bnds.clear();
    last_b = density->max_x();
    for( ; end_unknowns > beg_unknowns; --end_unknowns ) {
        const TF dx = sorted_dirac_positions[ end_unknowns - 1 ];
        const TF we = sorted_dirac_weights[ end_unknowns - 1 ];
        const TF dm = sorted_dirac_masses[ end_unknowns - 1 ];
        const TF b1 = dx + sqrt( we );
        if ( b1 < last_b )
            break;

        auto err = [&]( TF sub_b ) {
            return density->integral( last_b - sub_b, last_b ) - dm;
        };

        last_b -= dichotomy_growing_from_zero( err, target_max_error_ratio * dm, TF( 1 ) );
        bnds.push_back( last_b );
    }

    if ( end_unknowns < nb_sorted_diracs() ) {
        const TF dx = sorted_dirac_positions[ end_unknowns ];
        TF w = pow( bnds.back() - dx, 2 );
        bnds.pop_back();
        for( PI i = end_unknowns; i + 1 < nb_sorted_diracs(); ++i ) {
            const TF dx0 = sorted_dirac_positions[ i + 0 ];
            const TF dx1 = sorted_dirac_positions[ i + 1 ];
            const TF bnd = bnds.back();
            bnds.pop_back();

            sorted_dirac_weights[ i ] = w;

            // ( w1 - w0 ) / ( d1 - d0 ) = d1 + d0 - 2 * bnd;
            w += ( dx1 + dx0 - 2 * bnd ) * ( dx1 - dx0 );
        }

        sorted_dirac_weights.back() = w;            
    }
}

DTP void UTP::initialize_dirac_masses( const VF &relative_dirac_masses, TF global_mass_ratio ) {
    sorted_dirac_masses.resize( nb_sorted_diracs() );
    if ( relative_dirac_masses.empty() ) {
        const TF bm = global_mass_ratio / nb_original_diracs();
        for( PI i = 0; i < nb_sorted_diracs(); ++i )
            sorted_dirac_masses[ i ] = bm * ( sorted_dirac_num_offsets[ i + 1 ] - sorted_dirac_num_offsets[ i + 0 ] );
    } else {
        TF sum = 0;
        for( TF v : relative_dirac_masses )
            sum += v;
        
        const TF c = global_mass_ratio / sum;
        for( PI i = 0; i < nb_sorted_diracs(); ++i ) {
            TF acc = 0;
            for( PI j = sorted_dirac_num_offsets[ i + 0 ]; j < sorted_dirac_num_offsets[ i + 1 ]; ++j )                
                acc += relative_dirac_masses[ j ];
            sorted_dirac_masses[ i ] = c * acc;
        }
    }
}

DTP typename UTP::VF UTP::original_cell_barycenters() const {
    VF res( nb_original_diracs() );
    for_each_cell( [&]( PI n, TF b, TF e ) {
        const TF it = density->integral( b, e );
        for( PI o = sorted_dirac_num_offsets[ n + 0 ]; o < sorted_dirac_num_offsets[ n + 1 ]; ++o )
            res[ sorted_dirac_num_values[ o ] ] = it ? density->x_integral( b, e ) / it : 0;
    } );
    return res;
}


DTP void UTP::for_each_cell( auto &&func ) const {
    using namespace std;

    // no diracs ?
    const PI nc = sorted_dirac_positions.size();
    if ( nc == 0 )
        return;

    //
    TF c1 = std::numeric_limits<TF>::lowest();
    TF d1 = sorted_dirac_positions[ 0 ];
    TF w1 = sorted_dirac_weights[ 0 ];
    for( PI i = 1; i <= nc; ++i ) {
        // cuts
        const TF c0 = c1;
        const TF d0 = d1;
        const TF w0 = w1;
        if ( i < nc ) {
            d1 = sorted_dirac_positions[ i ];
            w1 = sorted_dirac_weights[ i ];
            c1 = ( d1 + d0 - ( w1 - w0 ) / ( d1 - d0 ) ) / 2;
        } else
            c1 = std::numeric_limits<TF>::max();

        // ball
        const TF we = max( TF( 0 ), sorted_dirac_weights[ i - 1 ] );
        const TF rb = sqrt( we );
        const TF b0 = d0 - rb;
        const TF b1 = d0 + rb;

        // call
        const TF e0 = max( c0, b0 );
        const TF e1 = min( c1, b1 );
        if ( e0 <= e1 )
            func( i - 1, e0, e1 );
    }
}

#ifdef TL_DISPLAYER_IS_DEFINED
DTP void UTP::display( Displayer &ds ) const {
    ds.start_array();
    for_each_cell( [&]( PI num_sorted_cell, TF x0, TF x1 ) {
        ds.start_array();
        ds << x0 << x1;
        ds.end_array();
    } );
    ds.end_array();
}
#endif

DTP typename UTP::VF UTP::sorted_cell_masses() const {
    using namespace std;
    VF res( nb_sorted_diracs() );
    for_each_cell( [&]( PI num_sorted_cell, TF x0, TF x1 ) {
        res[ num_sorted_cell ] = density->integral( x0, x1 );
    } );
    return res;
}

DTP void UTP::plot( Str filename ) const {
    using namespace std;
    
    ofstream fs( filename );

    fs << "from matplotlib import pyplot\n";

    for_each_cell( [&]( PI num_sorted_cell, TF x0, TF x1 ) {
        // x0 = max( x0, density->min_x() );
        // x1 = min( x1, density->max_x() );
        const TF y0 = ( num_sorted_cell - TF( 1 ) ) / nb_sorted_diracs();
        const TF y1 = ( num_sorted_cell + TF( 1 ) ) / nb_sorted_diracs();
        fs << "pyplot.plot( [ "
           << x0 << ", " << x1 << ", " << x1 << ", " << x0 << ", " << x0 << " ], [ "
           << y0 << ", " << y0 << ", " << y1 << ", " << y1 << ", " << y0 << " ] )\n";
    } );

    // density
    density->plot( fs );

    // // diracs
    // fs << "pyplot.plot( [ ";
    // for( auto c : dirac_positions() )
    //     fs << c << ", ";
    // fs << " ], [";
    // for( PI i = 0; i < nb_original_diracs(); ++i )
    //     fs << -0.5 << ", ";
    // fs << " ], '+' )\n";

    fs << "pyplot.show()\n";
}

DTP int UTP::newton_dir( VF &dir, TF &max_error_ratio ) const {
    using namespace std;

    // {
    //     TridiagonalSymmetricMatrix<TF> M;
    //     // TF max_error_ratio;
    //     VF V;
    //     int err = get_sorted_newton_system( M, V, max_error_ratio );
    //     if ( err )
    //         return err;
    //     M.inplace_ldlt_decomposition();
    //     dir = M.solve_using_ldlt( V );
    //     return 0;
    // }

    // init
    dir.resize( nb_sorted_diracs() );
    max_error_ratio = 0;
    
    // no diracs ?
    const PI nc = nb_sorted_diracs();
    if ( nc == 0 )
        return 0;

    // system intermediate values
    VF e( nb_sorted_diracs() );
    TF prev_ldds = 0;
    TF prev_v = 0;

    //
    TF c1 = std::numeric_limits<TF>::lowest();
    TF ld = std::numeric_limits<TF>::max();
    TF d1 = sorted_dirac_positions[ 0 ];
    TF w1 = sorted_dirac_weights[ 0 ];
    for( PI i = 1; i <= nc; ++i ) {
        // cuts
        const TF c0 = c1;
        const TF d0 = d1;
        const TF w0 = w1;
        if ( i < nc ) {
            d1 = sorted_dirac_positions[ i ];
            w1 = sorted_dirac_weights[ i ];
            c1 = ( d1 + d0 - ( w1 - w0 ) / ( d1 - d0 ) ) / 2;
        } else {
            d1 = std::numeric_limits<TF>::max();
            c1 = std::numeric_limits<TF>::max();
        }

        if ( c0 >= c1 )
            return 1;

        // ball
        const TF we = sorted_dirac_weights[ i - 1 ];
        if ( we <= 0 )
            return 1; // negative or zero weight
        const TF rb = sqrt( we );
        const TF b0 = d0 - rb;
        const TF b1 = d0 + rb;

        if ( b0 >= c1 || b1 <= c0 )
            return 2;

        // values
        TF m0 = 0;
        if ( b0 > c0 ) {
            m0 += density->value( b0 ) / rb / 2;
        } else {
            const TF d = density->value( c0 ) / ld / 2;
            m0 += d;
        }

        TF m1 = 0;
        const TF rd = d1 - d0;
        if ( b1 < c1 ) {
            m0 += density->value( b1 ) / rb / 2;
        } else {
            const TF d = density->value( c1 ) / rd / 2;
            m0 += d;
            m1 -= d;
        }

        const TF mass = density->integral( max( b0, c0 ), min( b1, c1 ) );
        const TF tgt = sorted_dirac_masses[ i - 1 ];
        const TF y = tgt - mass;

        // storage
        // if ( dmass_dL ) M( i - 1, i - 2 ) += dmass_dL;
        // M( i - 1, i - 1 ) += dmass_dM;

        max_error_ratio = max( max_error_ratio, abs( y ) / tgt );

        // system
        const TF d = m0 - prev_ldds;
        if ( d <= 0 )
            return 3;
        const TF l = m1 / d;
        e[ i - 1 ] = l;

        prev_ldds = d * l * l;

        const TF v0 = y - prev_v;
        dir[ i - 1 ] = v0 / d;
        prev_v = l * v0;

        // swap values
        ld = rd;
    }

    for( PI i = nb_sorted_diracs() - 1; i--; )
        dir[ i ] -= e[ i ] * dir[ i + 1 ];

    return 0;
}

DTP int UTP::get_sorted_newton_system( TridiagonalSymmetricMatrix<TF> &M, VF &V, TF &max_error_ratio ) const {
    using namespace std;
    
    // no diracs ?
    const PI nc = nb_sorted_diracs();
    if ( nc == 0 )
        return 0;

    // init
    M.resize( nc );
    V.resize( nc );
    M.clear_values();
    for( PI i = 0; i < nc; ++i )
        V[ i ] = sorted_dirac_masses[ i ];

    //
    TF c1 = std::numeric_limits<TF>::lowest();
    TF ld = std::numeric_limits<TF>::max();
    TF d1 = sorted_dirac_positions[ 0 ];
    TF w1 = sorted_dirac_weights[ 0 ];
    for( PI i = 1; i <= nc; ++i ) {
        // cuts
        const TF c0 = c1;
        const TF d0 = d1;
        const TF w0 = w1;
        if ( i < nc ) {
            d1 = sorted_dirac_positions[ i ];
            w1 = sorted_dirac_weights[ i ];
            c1 = ( d1 + d0 - ( w1 - w0 ) / ( d1 - d0 ) ) / 2;
        } else
            c1 = std::numeric_limits<TF>::max();

        if ( c0 >= c1 )
            return 1;

        // ball
        const TF we = max( TF( 0 ), sorted_dirac_weights[ i - 1 ] );
        const TF rb = sqrt( we );
        const TF b0 = d0 - rb;
        const TF b1 = d0 + rb;

        if ( b0 >= c1 || b1 <= c0 )
            return 2;

        // values
        TF dmass_dL = 0, dmass_dM = 0;
        if ( b0 > c0 ) {
            dmass_dM += density->value( b0 ) / rb / 2;
        } else {
            const TF d = density->value( c0 ) / ld / 2;
            dmass_dM += d;
            dmass_dL -= d;
        }

        const TF rd = d1 - d0;
        if ( b1 < c1 ) {
            dmass_dM += density->value( b1 ) / rb / 2;
        } else {
            const TF d = density->value( c1 ) / rd / 2;
            dmass_dM += d;
        }

        const TF mass = density->integral( max( b0, c0 ), min( b1, c1 ) );

        // storage
        if ( dmass_dL ) M( i - 1, i - 2 ) += dmass_dL;
        M( i - 1, i - 1 ) += dmass_dM;
        V[ i - 1 ] -= mass;

        // swap values
        ld = rd;
    }

    // 
    max_error_ratio = 0;
    for( PI i = 0; i < V.size(); ++i ) {
        const TF err = abs( V[ i ] ) / sorted_dirac_masses[ i ];
        max_error_ratio = max( max_error_ratio, err );

        // void cell
        if ( V[ i ] == sorted_dirac_masses[ i ] )
            return 4;
    }

    return 0;
}

DTP bool UTP::no_nan( const VF &v ) {
    for( TF x : v )
        if ( std::isnan( x ) || std::isinf( x ) )
            return false;
    return true;
}

DTP void UTP::plot_bnds_evolution( const std::vector<VB> &bnds ) {
    auto ys = linspace<TF>( 0, 1, bnds.size() );

    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    for( PI i = 0; i < bnds[ 0 ].size(); ++i ) {
        for( PI j = 0; j < bnds[ 0 ][ 0 ].size(); ++j ) {
            VF xs( bnds.size() );
            for( PI n = 0; n < bnds.size(); ++n )
                xs[ n ] = bnds[ n ][ i ][ j ];
            fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
        }
    }
    fs << "pyplot.legend()\n";
    fs << "pyplot.show()\n";
}

DTP typename UTP::VB UTP::sorted_cell_boundaries() const {
    VB res( nb_sorted_diracs() );
    for_each_cell( [&]( PI num_cell, TF x0, TF x1 ) {
        res[ num_cell ][ 0 ] = x0;
        res[ num_cell ][ 1 ] = x1;
    } );
    return res;
}

DTP int UTP::newton_iterations( const TF min_relax ) {
    // first newton system
    TF max_error_ratio;
    VF new_dir;
    VF dir;
    int err = newton_dir( dir, max_error_ratio );
    if ( err )
        return 1;
    
    // {
    //     P( dir, max_error_ratio );

    //     TridiagonalSymmetricMatrix<TF> M;

    //     VF V;
    //     get_sorted_newton_system( M, V, max_error_ratio );
    //     M.inplace_ldlt_decomposition();
    //     P( M.solve_using_ldlt( V ) );
    //     P( max_error_ratio );
    // }
    // assert( 0 );

    for( nb_newton_iterations = 0; ; ++nb_newton_iterations ) {
        if ( nb_newton_iterations == 500000 )
            throw std::runtime_error( "too many newton iterations" );
        
        // on a solution ?
        if ( max_error_ratio < target_max_error_ratio )
            return 0;    

        // find the relaxation coefficient
        VF w0 = sorted_dirac_weights;
        for( TF relax = 1; ; relax /= 4 ) {
            for( PI i = 0; i < w0.size(); ++i )
                sorted_dirac_weights[ i ] = w0[ i ] + relax * dir[ i ];

            err = newton_dir( new_dir, max_error_ratio );
            // err = get_sorted_newton_system( M, V, max_error_ratio );
            if ( err == 0 ) {
                if ( verbosity >= 2 )   
                    std::cout << "  relax: " << relax << " nb_newton_iterations: " << nb_newton_iterations << " max_error_ratio: " << max_error_ratio << "\n";
                std::swap( dir, new_dir );

                // {
                //     P( dir[ 10 ], max_error_ratio );
                //     if ( std::isnan( dir[ 10 ] ) || std::isinf( dir[ 10 ] ) )
                //         throw std::runtime_error( "NaN or Inf in dir" );

                //     TridiagonalSymmetricMatrix<TF> M;
                //     VF V;
                //     get_sorted_newton_system( M, V, max_error_ratio );
                //     M.inplace_ldlt_decomposition();
                //     P( M.solve_using_ldlt( V )[ 10 ], max_error_ratio );
                // }


                break;
            }

            if ( relax <= min_relax ) {
                sorted_dirac_weights = w0;
                return 2;
            }
        }
    }

    return 0;
}

DTP void UTP::solve() {
    newton_iterations();
}

#undef DTP
#undef UTP

} // namespace usdot
