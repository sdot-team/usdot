#pragma once


#include "WeightInitializer.h"
#include "utility/linspace.h"
//#include "WeightUpdater.h"
//#include "utility/glot.h"
#include "System.h"

#include <stdexcept>
#include <numeric>
#include <fstream>
#include <chrono>
#include <limits>

namespace usdot {
    
#define DTP template<class TF,class Density>
#define UTP System<TF,Density>

DTP UTP::System() {
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
        const TF d = aggs.back().end_x - b;
        if ( aggs.empty() || d < 0 ) {
            aggs.push_back( { b, b + l, l, i, i + 1 } );
            continue;
        }

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
    for( const Agglomerate &agg : aggs ) {
        if ( agg.beg_x > density->min_x() ) {
            PI i = agg.beg_n;

            TF d0 = sorted_dirac_positions[ i ];
            TF x0 = agg.beg_x;
            TF w0 = pow( d0 - x0, 2 );

            sorted_dirac_weights[ i ] = w0;

            while( ++i < agg.end_n ) {
                const TF d1 = sorted_dirac_positions[ i ];
                const TF x1 = x0 + sorted_dirac_masses[ i ] * density->ptp_x();
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
                const TF x0 = x1 - sorted_dirac_masses[ i ] * density->ptp_x();
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
            const TF x1 = x0 + sorted_dirac_masses[ i ] * density->ptp_x();
            const TF w1 = w0 + ( d1 - d0 ) * ( d0 + d1 - 2 * x1 );

            sorted_dirac_weights[ i ] = w1;
            
            d0 = d1;
            x0 = x1;
            w0 = w1;
        }
    }
}

DTP T_T void UTP::_for_each_newton_item( const T &func ) const {
    _for_each_unintersected_cell( [&]( TF dirac_position, TF dirac_weight, PI num_dirac, TF d0, TF d1, TF c0, TF c1 ) {
        const TF v = sorted_dirac_masses[ num_dirac ];
        if ( dirac_weight <= 0 )
            return func( num_dirac, 0, 0, v, 1, false );

        const TF rd = sqrt( dirac_weight );
        const TF b0 = dirac_position - rd;
        const TF b1 = dirac_position + rd;

        // "wrong" bounds ordering
        if ( c0 > c1 ) {
            // void cell
            if ( c0 < b0 || c1 > b1 ) {
                // P( c0, c1 );
                // plot();
                return func( num_dirac, 0, 0, v, 2, false );
            }

            // ball cut on the left
            if ( b0 > c1 ) {
                // ball cut on the right
                if ( b1 < c0 ) {
                    const TF db = ( density->value( b0 ) + density->value( b1 ) ) / ( 2 * sqrt( dirac_weight ) );
                    return func( num_dirac, 
                        db,
                        0,
                        - density->integral( b0, b1 ),
                        0,
                        db != 0
                    );
                }

                // interface on the right
                const TF db = density->value( b0 ) / sqrt( dirac_weight );
                return func( num_dirac, 
                    ( density->value( c0 ) / d0 - db ) / 2,
                    0,
                    - density->integral( b0, c0 ),
                    0,
                    db != 0
                );
            }
            
            // interface on the left
            if ( b1 < c0 ) { // ball cut on the right
                const TF cl = density->value( c1 ) / ( 2 * d1 );
                const TF db = density->value( b1 ) / ( 2 * sqrt( dirac_weight ) );
                return func( num_dirac,
                    cl - db,
                    - cl,
                    - density->integral( c1, b1 ),
                    0,
                    db != 0
                );
            }
            
            // interface on the right
            const TF cl = density->value( c1 ) / ( 2 * d1 );
            const TF cr = density->value( c0 ) / ( 2 * d0 );
            return func( num_dirac, 
                - ( cl + cr ),
                cr,
                - density->integral( c1, c0 ),
                0,
                false
            );
        }

        // void cell
        if ( c1 < b0 || c0 > b1 )
            return func( num_dirac, 0, 0, v, 3, false );

        if ( b0 > c0 ) { // ball cut on the left
            // ball cut on the right
            if ( b1 < c1 ) {
                const TF db = ( density->value( b0 ) + density->value( b1 ) ) / ( 2 * sqrt( dirac_weight ) );
                return func( num_dirac, 
                    db,
                    0,
                    density->integral( b0, b1 ),
                    0,
                    db != 0
                );
            }

            // interface on the right
            const TF cr = density->value( c1 ) / ( 2 * d1 );
            const TF db = density->value( b0 ) / ( 2 * sqrt( dirac_weight ) );
            return func( num_dirac, 
                db + cr,
                - cr,
                density->integral( b0, c1 ),
                0,
                db != 0
            );
        } 
        
        // interface on the left
        if ( b1 < c1 ) { // ball cut on the right
            const TF cl = density->value( c0 ) / ( 2 * d0 );
            const TF db = density->value( b1 ) / ( 2 * sqrt( dirac_weight ) );
            return func( num_dirac, 
                cl + db,
                0,
                density->integral( c0, b1 ),
                0,
                db != 0
            );
        }
        
        const TF cl = density->value( c0 ) / ( 2 * d0 );
        const TF cr = density->value( c1 ) / ( 2 * d1 );
        const TF in = density->integral( c0, c1 );
        return func( num_dirac, cl + cr, - cr, in, 0, false );
    } );
}

DTP T_T void UTP::_for_each_newton_item( PI num_der, const T &func ) {
    _for_each_unintersected_cell( [&]( TF dirac_position, TF dirac_weight, PI num_dirac, TF d0, TF d1, TF c0, TF c1 ) {
        const TF v = sorted_dirac_masses[ num_dirac ];
        if ( dirac_weight <= 0 )
            return func( num_dirac, 0, 0, 0, 1 );

        const TF rd = sqrt( dirac_weight );
        const TF b0 = dirac_position - rd;
        const TF b1 = dirac_position + rd;

        // "wrong" bounds ordering
        if ( c0 > c1 ) {
            // void cell
            if ( c0 < b0 || c1 > b1 )
                return func( num_dirac, 0, 0, 0, 2 );

            // ball cut on the left
            if ( b0 > c1 ) {
                // ball cut on the right
                if ( b1 < c0 ) {
                    return func( num_dirac, 
                        ( density->value( b0, num_der ) + density->value( b1, num_der ) ) / ( 2 * sqrt( dirac_weight ) ),
                        0, 
                        - density->integral( b0, b1, num_der ),
                        0
                    );
                }

                // interface on the right
                return func( num_dirac, 
                    ( density->value( c0, num_der ) / d0 - density->value( b0, num_der ) / sqrt( dirac_weight ) ) / 2,
                    0,
                    - density->integral( b0, c0, num_der ),
                    0
                );
            }
            
            // interface on the left
            if ( b1 < c0 ) { // ball cut on the right
                const TF cl = density->value( c1, num_der ) / ( 2 * d1 );
                return func( num_dirac, 
                    cl - density->value( b1, num_der ) / ( 2 * sqrt( dirac_weight ) ),
                    - cl,
                    - density->integral( c1, b1, num_der ),
                    0
                );
            }
            
            // interface on the right
            const TF cl = density->value( c1, num_der ) / ( 2 * d1 );
            const TF cr = density->value( c0, num_der ) / ( 2 * d0 );
            return func( num_dirac, 
                - ( cl + cr ),
                cr,
                - density->integral( c1, c0, num_der ),
                0
            );
        }

        // void cell
        if ( c1 < b0 || c0 > b1 )
            return func( num_dirac, 0, 0, 0, 3 );

        if ( b0 > c0 ) { // ball cut on the left
            // ball cut on the right
            if ( b1 < c1 ) { 
                return func( num_dirac, 
                    ( density->value( b0, num_der ) + density->value( b1, num_der ) ) / ( 2 * sqrt( dirac_weight ) ),
                    0,
                    density->integral( b0, b1, num_der ),
                    0
                );
            }

            // interface on the right
            const TF cr = density->value( c1, num_der ) / ( 2 * d1 );
            return func( num_dirac, 
                density->value( b0, num_der ) / ( 2 * sqrt( dirac_weight ) ) + cr,
                - cr,
                density->integral( b0, c1, num_der ),
                0
            );
        } 
        
        // interface on the left
        if ( b1 < c1 ) { // ball cut on the right
            const TF cl = density->value( c0, num_der ) / ( 2 * d0 );
            return func( num_dirac, 
                cl + density->value( b1, num_der ) / ( 2 * sqrt( dirac_weight ) ),
                0,
                density->integral( c0, b1, num_der ),
                0
            );
        }
        
        const TF cl = density->value( c0, num_der ) / ( 2 * d0 );
        const TF cr = density->value( c1, num_der ) / ( 2 * d1 );
        const TF in = density->integral( c0, c1, num_der );
        return func( num_dirac, cl + cr, - cr, in, 0 );
    } );
}

DTP UTP::MF UTP::der_weights_wrt_lap_ratio( PI nb_ders ) {
    const TF eps = 1e-6;
    
    VF V( nb_sorted_diracs() );
    MF res;

    TF base_fr = density->current_flattening_ratio;
    VF w0 = sorted_dirac_weights;

    density->set_flattening_ratio( base_fr - 1 * eps );
    int err = newton_iterations();
    if ( err )
        throw std::runtime_error( "err newton in der w" );
    VF w1 = sorted_dirac_weights;

    density->set_flattening_ratio( base_fr - 2 * eps );
    err = newton_iterations();
    if ( err )
        throw std::runtime_error( "err newton in der w" );
    VF w2 = sorted_dirac_weights;

    VF d1( w0.size() );
    for( PI i = 0; i < w0.size(); ++i )
        d1[ i ] = ( w0[ i ] - w1[ i ] ) / eps; 
    res.push_back( d1 );

    VF d2( w0.size() );
    for( PI i = 0; i < w0.size(); ++i )
        d1[ i ] = ( w0[ i ] + w2[ i ] - 2 * w1[ i ] ) / pow( eps, 2 ); 
    res.push_back( d2 );

    // // helpers
    // auto set_vec = [&]( PI nd ) {
    //     _for_each_newton_item( nd, [&]( PI index, TF m0, TF m1, TF v, int ) {
    //         V[ index ] = - v;
    //     } );
    // };
    // auto add_mat = [&]( PI nd, PI nr, TF coeff ) {
    //     _for_each_newton_item( nd, [&]( PI index, TF m0, TF m1, TF v, int ) {
    //         V[ index ] -= 2 * m0 * res[ nr ][ index ];
    //         if ( m1 ) {
    //             V[ index ] -= 2 * m1 * res[ nr ][ index + 1 ];
    //             V[ index + 1 ] -= 2 * m1 * res[ nr ][ index ];
    //         }
    //     } );
    // };

    // // X'
    // if ( nb_ders >= 1 ) {
    //     set_vec( 1 );
    //     res.push_back( newton_matrix_ldlt.solve_using_ldlt( V ) );
    // }

    // // X''
    // if ( nb_ders >= 2 ) {
    //     set_vec( 2 );
    //     add_mat( 1, 0, 2 );
    //     res.push_back( newton_matrix_ldlt.solve_using_ldlt( V ) );
    // }

    // // X''', ...
    // if ( nb_ders >= 3 ) {
    //     set_vec( 3 );
    //     add_mat( 2, 0, 3 );
    //     add_mat( 1, 1, 3 );
    //     res.push_back( newton_matrix_ldlt.solve_using_ldlt( V ) );
    // }

    // // X'''', ...
    // if ( nb_ders >= 4 ) {
    //     set_vec( 4 );
    //     add_mat( 3, 0, 4 );
    //     add_mat( 2, 1, 6 );
    //     add_mat( 1, 2, 4 );
    //     res.push_back( newton_matrix_ldlt.solve_using_ldlt( V ) );
    // }

    // // ...
    // if ( nb_ders >= 5 ) {
    //     throw std::runtime_error( "TODO: nb_ders >= 4" );
    // }

    // //
    return res;
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

DTP int UTP::_make_newton_system() {
    using namespace std;

    newton_matrix_ldlt.resize( nb_sorted_diracs() );
    newton_vector.resize( nb_sorted_diracs() );
    newton_error = 0;

    bool has_bad_cell = false;
    bool has_arc = false;
    _for_each_newton_item( [&]( PI index, TF m0, TF m1, TF v, int bad_cell, bool arc ) {
        if ( bad_cell )
            has_bad_cell = bad_cell;
        if ( arc )
            has_arc = arc;

        newton_matrix_ldlt( index, index ) = m0;
        if ( m1 )
            newton_matrix_ldlt( index, index + 1 ) = m1;

        v -= sorted_dirac_masses[ index ];
        newton_vector[ index ] = v;
        newton_error += v * v;
    } );

    if ( has_bad_cell )
        return 1;

    if ( ! has_arc )
        newton_matrix_ldlt( 0, 0 ) *= 2;

    newton_matrix_ldlt.inplace_ldlt_decomposition();

    return 0;
}

DTP int UTP::newton_iterations() {
    TF prev_newton_error = std::numeric_limits<TF>::max();
    for( PI num_iter = 0; ; ++num_iter ) {
        if ( num_iter == 30 )
            throw std::runtime_error( "max iter" );

        int err = _make_newton_system();
        // P( num_iter, err, newton_error );
        if ( err || newton_error > prev_newton_error )
            return 1;
        prev_newton_error = newton_error;
        
        VF v = newton_matrix_ldlt.solve_using_ldlt( newton_vector );
        for( PI i = 0; i < v.size(); ++i )
            sorted_dirac_weights[ i ] -= v[ i ];
    
        if ( newton_error < 1e-20 ) {
            P( num_iter );
            return 0;
        }
    }
}

DTP void UTP::solve_using_cdf() {
    // _update_system( true );

    // auto t0 = std::chrono::high_resolution_clock::now();

    // WeightInitializer<TF,Density> wi( *this );
    // wi.max_nb_iterations = 800000;
    // wi.run();

    // nb_iterations_init = wi.nb_iterations;
    // auto t1 = std::chrono::high_resolution_clock::now();
    // time_in_init = std::chrono::duration<double>{ t1 - t0 }.count();
}

DTP void UTP::solve() {
    initialize_with_flat_density();
    newton_iterations();

    //     
    for( TF prev_flattening_ratio = 1; prev_flattening_ratio; ) {
        VF w0 = sorted_dirac_weights;
        MF ders = der_weights_wrt_lap_ratio( 3 );
        TF trat = 0.75;
        for( TF trial_flattening_ratio = 0; ; trial_flattening_ratio = ( 1 - trat ) * prev_flattening_ratio + trat * trial_flattening_ratio ) {
            const TF a = trial_flattening_ratio - prev_flattening_ratio;
            
            P( prev_flattening_ratio, trial_flattening_ratio, a );
            
            density->set_flattening_ratio( trial_flattening_ratio );
            for( PI i = 0; i < nb_sorted_diracs(); ++i )
                sorted_dirac_weights[ i ] = w0[ i ] + ders[ 0 ][ i ] * a + ders[ 1 ][ i ] * pow( a, 2 ) / 2; // + ders[ 2 ][ i ] * pow( a, 3 ) / 6;

            if ( abs( a ) < 1e-6 ) {
                throw std::runtime_error( "TODO: low a" );
            }

            // test 
            int err = newton_iterations();
            if ( err == 0 ) {
                prev_flattening_ratio = trial_flattening_ratio;
                break;
            }
        }
    }

    // if ( verbosity >= 2 && stream )
    //     *stream << "nb iteration init: " << nb_iterations_init << " update: " << nb_iterations_update << "\n";
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
        const TF y0 = ( y++ ) / nb_sorted_diracs();
        const TF y1 = y0 + 1;
        fs << "pyplot.plot( [ " << x0 << ", " << x1 << ", " << x1 << ", " << x0 << ", " << x0 << " ], ";
        fs << "[ " << y0 << ", " << y0 << ", " << y1 << ", " << y1 << ", " << y0 << " ] )\n";
    }


    // diracs
    fs << "pyplot.plot( [ ";
    for( auto c : dirac_positions() )
        fs << c << ", ";
    fs << " ], [";
    for( auto c : dirac_positions() )
        fs << 0 << ", ";
    fs << " ], '+' )\n";

    fs << "pyplot.show()\n";
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

DTP TF UTP::l2_mass_error( bool max_if_bad_cell ) const {
    using namespace std;
    _update_system();
    TF res = 0;
    bool has_bad_cell = false;
    for_each_cell( [&]( PI n, TF b, TF e ) {
        if ( b >= e )
            has_bad_cell = true;
        res += pow( sorted_dirac_masses[ n ] - density->integral( b, e ), 2 );
    } );
    if ( has_bad_cell )    
        return numeric_limits<TF>::max();
    return sqrt( res );
}

DTP TF UTP::max_mass_error() const {
    using namespace std;
    _update_system();
    TF res = 0;
    for_each_cell( [&]( PI n, TF b, TF e ) {
        res = max( res, abs( sorted_dirac_masses[ n ] - density->integral( b, e ) ) );
    } );
    return sqrt( res );
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
