#pragma once


#include <eigen3/Eigen/SparseLU>
#include <eigen3/Eigen/Dense>
#include "utility/linspace.h"
#include "TestSystem.h"

#include <stdexcept>
#include <numeric>
#include <fstream>
#include <limits>
#include <utility>
#include <vector>

namespace usdot {
    
#define DTP template<class TF,class Density>
#define UTP TestSystem<TF,Density>

DTP UTP::TestSystem() {
    global_mass_ratio = 1;
    density = nullptr;
}

DTP void UTP::initialize_with_flat_density() {
    _update_system( /*need_weights*/ false );
    using namespace std;

    density->set_flattening_ratio( 1 );

    // make the agglomerates
    struct Agglomerate {
        TF beg_x;
        TF end_x;
        TF len_x;
        PI beg_n;
        PI end_n;
    };
    std::vector<Agglomerate> aggs;
    aggs.reserve( nb_sorted_diracs() / 2 );
    for( PI i = 0; i < nb_sorted_diracs(); ++i ) {
        const TF l = sorted_dirac_masses[ i ] * density->ptp_x();
        const TF b = min( density->max_x() - l, max( density->min_x(), sorted_dirac_positions[ i ] - l / 2 ) );
        if ( aggs.empty() || aggs.back().end_x - b < 0 ) {
            aggs.push_back( { b, b + l, l, i, i + 1 } );
            continue;
        }
        const TF d = aggs.back().end_x - b;

        Agglomerate &agg = aggs.back();

        agg.len_x += l;
        agg.beg_x -= d * l / agg.len_x;
        agg.end_x = agg.beg_x + agg.len_x;
        agg.end_n = i + 1;
        if ( agg.beg_x < density->min_x() ) {
            agg.beg_x = density->min_x();
            agg.end_x = density->min_x() + agg.len_x;
        }
        if ( agg.end_x > density->max_x() ) {
            agg.beg_x = density->max_x() - agg.len_x;
            agg.end_x = density->max_x();
        }

        while ( true ) {
            if ( aggs.size() < 2 )
                break;
            
            Agglomerate &a0 = aggs[ aggs.size() - 2 ];
            Agglomerate &a1 = aggs[ aggs.size() - 1 ];
            const TF d = a0.end_x - a1.beg_x;
            if ( d < 0 )
                break;

            a0.len_x += a1.len_x;
            a0.beg_x -= d * a1.len_x / a0.len_x;
            a0.end_x = a0.beg_x + a0.len_x;
            a0.end_n = a1.end_n;
            if ( a0.beg_x < density->min_x() ) {
                a0.beg_x = density->min_x();
                a0.end_x = density->min_x() + a0.len_x;
            }
            if ( a0.end_x > density->max_x() ) {
                a0.beg_x = density->max_x() - a0.len_x;
                a0.end_x = density->max_x();
            }
        
            aggs.pop_back();
        }
    }

    // set the weights
    sorted_dirac_weights.resize( nb_sorted_diracs() );
    for( const Agglomerate &agg : aggs ) {
        P( agg.beg_x, agg.end_x, agg.len_x );
        if ( agg.beg_x > density->min_x() ) {
            PI i = agg.beg_n;

            TF d0 = sorted_dirac_positions[ i ];
            TF x0 = agg.beg_x;
            TF w0 = pow( d0 - x0, 2 );

            sorted_dirac_weights[ i ] = w0;

            while( ++i < agg.end_n ) {
                const TF d1 = sorted_dirac_positions[ i ];
                const TF x1 = x0 + sorted_dirac_masses[ i - 1 ] * density->ptp_x();
                const TF w1 = w0 + ( d1 - d0 ) * ( d0 + d1 - 2 * x1 );

                sorted_dirac_weights[ i ] = w1;

                d0 = d1;
                x0 = x1;
                w0 = w1;
            }

            continue;
        }

        if ( agg.end_x < density->max_x() ) {
            PI i = agg.end_n - 1;

            TF d1 = sorted_dirac_positions[ i ];
            TF x1 = agg.end_x;
            TF w1 = pow( x1 - d1, 2 );
            
            sorted_dirac_weights[ i ] = w1;

            while( i-- > agg.beg_n ) {
                const TF d0 = sorted_dirac_positions[ i ];
                const TF x0 = x1 - sorted_dirac_masses[ i + 1 ] * density->ptp_x();
                const TF w0 = w1 - ( d1 - d0 ) * ( d0 + d1 - 2 * x0 );

                sorted_dirac_weights[ i ] = w0;

                d1 = d0;
                x1 = x0;
                w1 = w0;
            }

            continue;
        }

        //
        PI i = agg.beg_n;

        TF d0 = sorted_dirac_positions[ i ];
        TF x0 = density->min_x();
        TF w0 = 2 * max(
            pow( sorted_dirac_positions[ agg.end_n - 1 ] - density->min_x(), 2 ),
            pow( density->max_x() - sorted_dirac_positions[ agg.beg_n ], 2 )
        );
        sorted_dirac_weights[ i ] = w0;
        
        while( ++i < agg.end_n ) {
            const TF d1 = sorted_dirac_positions[ i ];
            const TF x1 = x0 + sorted_dirac_masses[ i - 1 ] * density->ptp_x();
            const TF w1 = w0 + ( d1 - d0 ) * ( d0 + d1 - 2 * x1 );

            sorted_dirac_weights[ i ] = w1;
            
            d0 = d1;
            x0 = x1;
            w0 = w1;
        }
    }
}

DTP T_T void UTP::_for_each_unintersected_cell( const T &func ) const {
    using namespace std;

    TF i0 = numeric_limits<TF>::lowest();
    TF d0 = sorted_dirac_positions[ 0 ];
    TF w0 = sorted_dirac_weights[ 0 ];
    TF od = numeric_limits<TF>::max();
    for( PI n = 1; n < nb_sorted_diracs(); ++n ) {
        const TF d1 = sorted_dirac_positions[ n ];
        const TF w1 = sorted_dirac_weights[ n ];

        const TF i1 = ( d0 + d1 + ( w0 - w1 ) / ( d1 - d0 ) ) / 2;
        func( d0, w0, n - 1, od, d1 - d0, i0, i1 );

        od = d1 - d0;
        i0 = i1;
        d0 = d1;
        w0 = w1;
    }

    func( d0, w0, nb_sorted_diracs() - 1, od, numeric_limits<TF>::max(), i0, numeric_limits<TF>::max() );
}

DTP typename UTP::VF UTP::newton_dir_ap() {    
    using namespace std;

    // M
    // using TM = Eigen::Matrix<TF,Eigen::Dynamic,Eigen::Dynamic>;
    using TV = Eigen::Matrix<TF,Eigen::Dynamic,1>;
    using TM = Eigen::SparseMatrix<TF>;
    TM M( nb_sorted_diracs(), nb_sorted_diracs() );    
    VF V = mass_errors();

    // M.fill( 0 );
    for( PI r = 0; r < nb_sorted_diracs(); ++r ) {
        TF &ref_weight = sorted_dirac_weights[ r ];
        TF old_weight = ref_weight;
        ref_weight += eps;

        VF U = mass_errors();

        for( PI c = max( r, 1ul ) - 1; c < min( r + 2, nb_sorted_diracs() ); ++c )
        // for( PI c = 0; c < nb_sorted_diracs(); ++c )
            M.coeffRef( r, c ) = ( U[ c ] - V[ c ] ) / eps;

        ref_weight = old_weight;
    }

    TV Y = Eigen::Map<TV,Eigen::Unaligned>( V.data(), V.size() );

    // Eigen::FullPivLU<TM> lu( M );
    // auto X = lu.solve( Y );
    Eigen::SparseLU<TM,Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern( M ); 
    solver.factorize( M ); 
    auto X = solver.solve( Y ); 

    // Eigen::EigenSolver<TM> es( M );
    // auto ev = es.eigenvalues();
    // std::cout << "M: " << M << std::endl;
    // std::cout << "ev: " << ev << std::endl;

    // VF eigs = VF{ ev.begin(), ev.end() };
    // std::sort( eigs.begin(), eigs.end() );
    // P( eigs );

    return VF{ X.begin(), X.end() };
}

DTP int UTP::newton_iterations( TF min_relax ) {
    // 
    for( nb_newton_iterations = 0; ; ++nb_newton_iterations ) {
        if ( nb_newton_iterations == 50 )
            throw std::runtime_error( "max iter" );

        auto dir = newton_dir_ap();
        VF w0 = sorted_dirac_weights;
        for( TF a = 1; ; a /= 2 ) {
            if ( a < min_relax ) {
                sorted_dirac_weights = w0;
                // plot_bnds_evolution( bnds_evolution );
                // P( cell_boundaries() );
                // plot();
                return 2;
            }

            for( PI i = 0; i < w0.size(); ++i )
                sorted_dirac_weights[ i ] = w0[ i ] - a * dir[ i ];

            if ( max_relative_mass_error() < std::numeric_limits<TF>::max() ) {
                P( nb_newton_iterations, a, max_relative_mass_error() );
                bnds_evolution.push_back( cell_boundaries() );
                break;
            }
        }
            
        if ( max_relative_mass_error() < 1e-3 )
            return 0;
    }
}

DTP void UTP::solve( bool use_approx_for_ders ) {
    using namespace std;

    initialize_with_flat_density();
    
    for( TF d : { 0.999, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0 } ) {
        density->set_flattening_ratio( d );

        // glot( linspace<TF>( density->min_x() - 0, density->max_x() + 0, 1000 ), 
        //     [&]( TF x ) { return density->value( x ); }
        // );
        plot();

        int err = newton_iterations( 1e-20 );
        if ( err ) {
            P( err );
            throw runtime_error( "err 1" );
        }
    }
}

DTP PI UTP::nb_original_diracs() const {
    return sorted_dirac_num_values.size();
}

DTP PI UTP::nb_sorted_diracs() const {
    return sorted_dirac_positions.size();
}

DTP TF UTP::density_value( TF pos ) const {
    return density->value( pos );
}

DTP typename UTP::VF UTP::dirac_positions() const {
    VF res( nb_original_diracs() );
    for( PI i = 0; i < nb_sorted_diracs(); ++i )
        for( PI o = sorted_dirac_num_offsets[ i + 0 ]; o < sorted_dirac_num_offsets[ i + 1 ]; ++o )
            res[ sorted_dirac_num_values[ o ] ] = sorted_dirac_positions[ i ];
    return res;
}

DTP typename UTP::VF UTP::dirac_weights() const {
    VF res( nb_original_diracs() );
    for( PI i = 0; i < nb_sorted_diracs(); ++i )
        for( PI o = sorted_dirac_num_offsets[ i + 0 ]; o < sorted_dirac_num_offsets[ i + 1 ]; ++o )
            res[ sorted_dirac_num_values[ o ] ] = sorted_dirac_weights[ i ];
    return res;
}

DTP TF UTP::x_tol() const {
    return target_max_mass_error * global_mass_ratio * density->ptp_x() / nb_sorted_diracs();
}

DTP void UTP::plot( Str filename ) const {
    std::ofstream fs( filename );

    fs << "from matplotlib import pyplot\n";

    // density
    density->plot( fs );

    // boundaries
    TF y = 0;
    for( auto c : cell_boundaries() ) {
        const TF x0 = c[ 0 ];
        const TF x1 = c[ 1 ];
        const TF y0 = ( y++ ) / nb_sorted_diracs() / 20;
        const TF y1 = y0 + 0.5;
        fs << "pyplot.plot( [ " << x0 << ", " << x1 << ", " << x1 << ", " << x0 << ", " << x0 << " ], ";
        fs << "[ " << y0 << ", " << y0 << ", " << y1 << ", " << y1 << ", " << y0 << " ] )\n";
    }


    // diracs
    fs << "pyplot.plot( [ ";
    for( auto c : dirac_positions() )
        fs << c << ", ";
    fs << " ], [";
    for( auto c : dirac_positions() )
        fs << -0.5 << ", ";
    fs << " ], '+' )\n";

    fs << "pyplot.show()\n";
}

DTP void UTP::plot_bnds_evolution( const std::vector<VB> &bnds ) {

    auto ys = linspace<TF>( 0, 1, bnds.size() );

    std::ofstream fs( "glot.py" );
    fs << "from matplotlib import pyplot\n";
    if ( bnds.size() ) {
        for( PI i = 0; i < bnds[ 0 ].size(); ++i ) {
            for( PI j = 0; j < bnds[ 0 ][ 0 ].size(); ++j ) {
                VF xs( bnds.size() );
                for( PI n = 0; n < bnds.size(); ++n )
                    xs[ n ] = bnds[ n ][ i ][ j ];
                fs << "pyplot.plot( " << to_string( xs ) << ", " << to_string( ys ) << " )\n";
            }
        }
    }
    fs << "pyplot.legend()\n";
    fs << "pyplot.show()\n";
}

DTP T_T void UTP::for_each_cell( const T &func ) const {
    using namespace std;

    if ( nb_sorted_diracs() == 0 )
        return;

    TF i0 = numeric_limits<TF>::lowest();
    TF d0 = sorted_dirac_positions[ 0 ];
    TF w0 = sorted_dirac_weights[ 0 ];

    for( PI n = 1; n < nb_sorted_diracs(); ++n ) {
        const TF d1 = sorted_dirac_positions[ n ];
        const TF w1 = sorted_dirac_weights[ n ];
        const TF i1 = ( d0 + d1 + ( w0 - w1 ) / ( d1 - d0 ) ) / 2;
        const TF r0 = sqrt( max( w0, TF( 0 ) ) );
        const TF b0 = d0 - r0;
        const TF e0 = d0 + r0;

        func( n - 1, max( i0, b0 ), min( i1, e0 ) );

        i0 = i1;
        d0 = d1;
        w0 = w1;
    }

    const TF r0 = sqrt( max( w0, TF( 0 ) ) );
    const TF b0 = d0 - r0;
    const TF e0 = d0 + r0;

    func( nb_sorted_diracs() - 1, max( i0, b0 ), e0 );
}

DTP typename UTP::VF UTP::cell_barycenters() const {
    _update_system();

    VF res( nb_original_diracs() );
    for_each_cell( [&]( PI n, TF b, TF e ) {
        const TF it = density->integral( b, e );
        for( PI o = sorted_dirac_num_offsets[ n + 0 ]; o < sorted_dirac_num_offsets[ n + 1 ]; ++o )
            res[ sorted_dirac_num_values[ o ] ] = it ? density->x_integral( b, e ) / it : 0;
    } );
    return res;
}

DTP typename UTP::VB UTP::cell_boundaries() const {
    using namespace std;

    _update_system();

    VB res( nb_original_diracs() );
    TF mi = density->min_x();
    TF ma = density->max_x();
    for_each_cell( [&]( PI n, TF b, TF e ) {
        b = max( mi, min( ma, b ) );
        e = max( mi, min( ma, e ) );
        for( PI o = sorted_dirac_num_offsets[ n + 0 ]; o < sorted_dirac_num_offsets[ n + 1 ]; ++o )
            res[ sorted_dirac_num_values[ o ] ] = { b, e }; 
    } );
    return res;
}

DTP TF UTP::cost() const {
    using namespace std;
    _update_system();
    TF res = 0;
    for_each_cell( [&]( PI n, TF b, TF e ) {
        res += density->x2_integral( b, e, sorted_dirac_positions[ n ] );
    } );
    return sqrt( res );
}

DTP typename UTP::VF UTP::cell_masses() const {
    _update_system();
    VF res( nb_original_diracs() );
    for_each_cell( [&]( PI n, TF b, TF e ) {
        PI bo = sorted_dirac_num_offsets[ n + 0 ];
        PI eo = sorted_dirac_num_offsets[ n + 1 ];
        for( PI o = bo; o < eo; ++o )
            res[ sorted_dirac_num_values[ o ] ] = density->integral( b, e ) / ( eo - bo );
    } );
    return res;
}

// DTP TF UTP::mass_error( bool max_if_bad_cell ) const {
//     using namespace std;
//     _update_system();

//     TF res = 0;
//     bool has_bad_cell = false;
//     for_each_cell( [&]( PI n, TF b, TF e ) {
//         if ( b >= e )
//             has_bad_cell = true;
//         const TF in = density->integral( b, e );
//         const TF ma = sorted_dirac_masses[ n ];
//         auto err = in / ma - inv_coeff * ma / in;
//         res += pow( err, 2 );
//     } );
//     if ( has_bad_cell )    
//         return numeric_limits<TF>::max();
//     return res;
// }

DTP typename UTP::VF UTP::mass_errors() const {
    using namespace std;
    _update_system();

    VF res( nb_sorted_diracs() );
    for_each_cell( [&]( PI n, TF b, TF e ) {
        const TF in = density->integral( b, e );
        const TF ma = sorted_dirac_masses[ n ];
        // res[ n ] = in / ma - inv_coeff * ma / in - ( 1 - inv_coeff );
        res[ n ] = in - ma;
    } );
    return res;
    // VF res( nb_sorted_diracs() );
    // TF V = mass_error();
    // for( PI r = 0; r < nb_sorted_diracs(); ++r ) {
    //     TF &ref_weight = sorted_dirac_weights[ r ];
    //     TF old_weight = ref_weight;
    //     ref_weight += eps;

    //     TF U = mass_error();
    //     res[ r ] = ( U - V ) / eps;

    //     ref_weight = old_weight;
    // }
    // return res;
}

DTP TF UTP::max_relative_mass_error() const {
    using namespace std;
    _update_system();
    TF res = 0;
    bool has_bad_cell = 0;
    for_each_cell( [&]( PI n, TF b, TF e ) {
        const TF m = density->integral( b, e );
        if ( m <= 0 )
            has_bad_cell = 1;
        res = max( res, abs( sorted_dirac_masses[ n ] - m ) / sorted_dirac_masses[ n ] );
    } );
    P( has_bad_cell );
    if ( has_bad_cell )
        return numeric_limits<TF>::max();
    return res;
}

DTP void UTP::set_dirac_positions( const VF &dirac_positions, TF min_dirac_separation ) {
    using namespace std;

    // sorted_dirac_nums
    VI sorted_dirac_nums( dirac_positions.size() );
    iota( sorted_dirac_nums.begin(), sorted_dirac_nums.end(), PI( 0 ) );
    sort( sorted_dirac_nums.begin(), sorted_dirac_nums.end(), [&]( PI a, PI b ) {
        return dirac_positions[ a ] < dirac_positions[ b ];
    } );

    // sorted_dirac_positions
    // const TF eps = numeric_limits<TF>::epsilon() * 100 * density->width();
    sorted_dirac_num_offsets.reserve( dirac_positions.size() );
    sorted_dirac_num_values.reserve( dirac_positions.size() );
    sorted_dirac_positions.reserve( dirac_positions.size() );
    for( PI i = 0; i < sorted_dirac_nums.size(); ++i ) {
        const PI s = sorted_dirac_num_values.size();
        const PI n = sorted_dirac_nums[ i ];
        const TF p = dirac_positions[ n ];
        if ( s && sorted_dirac_positions.back() + min_dirac_separation >= p ) {
            sorted_dirac_num_values.push_back( n );
            continue;
        }

        sorted_dirac_num_offsets.push_back( s );
        sorted_dirac_num_values.push_back( n );
        sorted_dirac_positions.push_back( p );
    }
    sorted_dirac_num_offsets.push_back( sorted_dirac_num_values.size() );
}

DTP void UTP::set_relative_dirac_masses( const VF &relative_mass_ratios ) {
    this->relative_mass_ratios = relative_mass_ratios;
}

DTP void UTP::set_global_mass_ratio( const TF &global_mass_ratio ) {
    this->global_mass_ratio = global_mass_ratio;
}

DTP void UTP::set_density( Density *density ) {
    this->density = density;
}

DTP void UTP::_update_system( bool need_weights ) const {
    if ( ! density )
        throw std::runtime_error( "density is not defined" );
    
    if ( sorted_dirac_masses.empty() ) {
        const TF bm = global_mass_ratio / nb_original_diracs();
        sorted_dirac_masses.resize( nb_sorted_diracs() );
        for( PI i = 0; i < nb_sorted_diracs(); ++i )
            sorted_dirac_masses[ i ] = bm * ( sorted_dirac_num_offsets[ i + 1 ] - sorted_dirac_num_offsets[ i + 0 ] );
    }

    if ( sorted_dirac_weights.empty() ) {
        const TF b = sorted_dirac_positions.front();
        const TF e = sorted_dirac_positions.back();
        const TF w = pow( global_mass_ratio * ( e - b ) / nb_sorted_diracs(), 2 );
        sorted_dirac_weights = fill<TF>( nb_sorted_diracs(), w );
    }
}

#undef DTP
#undef UTP

} // namespace usdot
