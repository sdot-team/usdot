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
#include <vector>

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
        // const TF v = sorted_dirac_masses[ num_dirac ];
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

DTP typename UTP::MF UTP::der_weights_wrt_flat_ratio( PI nb_ders, bool use_approx_for_ders ) {
    MF res;
    if ( use_approx_for_ders ) {
        TF base_fr = density->current_flattening_ratio;
        VF w0 = sorted_dirac_weights;
        const TF eps = 1e-20;

        density->set_flattening_ratio( base_fr - 1 * eps );
        int err = newton_iterations( 1e-6 );
        if ( err )
            throw std::runtime_error( "err newton in der w" );
        VF w1 = sorted_dirac_weights;

        density->set_flattening_ratio( base_fr - 2 * eps );
        err = newton_iterations( 1e-6 );
        if ( err )
            throw std::runtime_error( "err newton in der w" );
        VF w2 = sorted_dirac_weights;

        if ( nb_ders >= 3 ) {
            density->set_flattening_ratio( base_fr - 3 * eps );
            err = newton_iterations( 1e-6 );
            if ( err )
                throw std::runtime_error( "err newton in der w" );
        }
        VF w3 = sorted_dirac_weights;


        // ----------------------
        VF d1( w0.size() );
        for( PI i = 0; i < w0.size(); ++i )
            d1[ i ] = ( w0[ i ] - w1[ i ] ) / eps; 
        res.push_back( d1 );

        VF d2( w0.size() );
        for( PI i = 0; i < w0.size(); ++i )
            d2[ i ] = ( w0[ i ] + w2[ i ] - 2 * w1[ i ] ) / pow( eps, 2 ); 
        res.push_back( d2 );

        if ( nb_ders >= 3 ) {
            VF d3( w0.size() );
            for( PI i = 0; i < w0.size(); ++i )
                d3[ i ] = ( w0[ i ] - 3 * w1[ i ] + 3 * w2[ i ] - w3[ i ] ) / pow( eps, 3 ); 
            res.push_back( d3 );
        }

        return res;
    }

    //
    density->compute_derivatives( nb_ders );
    VF V( nb_sorted_diracs() );

    // helpers
    auto set_vec = [&]( PI nd ) {
        _for_each_newton_item( nd, [&]( PI index, TF m0, TF m1, TF v, int ) {
            V[ index ] = - v;
        } );
    };
    auto add_mat = [&]( PI nd, PI nr, TF coeff ) {
        _for_each_newton_item( nd, [&]( PI index, TF m0, TF m1, TF v, int ) {
            V[ index ] -= 2 * m0 * res[ nr ][ index ];
            if ( m1 ) {
                V[ index ] -= 2 * m1 * res[ nr ][ index + 1 ];
                V[ index + 1 ] -= 2 * m1 * res[ nr ][ index ];
            }
        } );
    };

    // X'
    if ( nb_ders >= 1 ) {
        set_vec( 1 );
        res.push_back( newton_matrix_ldlt.solve_using_ldlt( V ) );
    }

    // X''
    if ( nb_ders >= 2 ) {
        set_vec( 2 );
        add_mat( 1, 0, 2 );
        res.push_back( newton_matrix_ldlt.solve_using_ldlt( V ) );
    }

    // X''', ...
    if ( nb_ders >= 3 ) {
        set_vec( 3 );
        add_mat( 2, 0, 3 );
        add_mat( 1, 1, 3 );
        res.push_back( newton_matrix_ldlt.solve_using_ldlt( V ) );
    }
    
    // X'''', ...
    if ( nb_ders >= 4 ) {
        set_vec( 4 );
        add_mat( 3, 0, 4 );
        add_mat( 2, 1, 6 );
        add_mat( 1, 2, 4 );
        res.push_back( newton_matrix_ldlt.solve_using_ldlt( V ) );
    }

    // ...
    if ( nb_ders >= 5 ) {
        throw std::runtime_error( "TODO: nb_ders >= 4" );
    }

    //
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

    newton_matrix_ldlt.clear_values();

    bool has_bad_cell = false;
    bool has_arc = false;
    TF mm0 = numeric_limits<TF>::max();
    _for_each_newton_item( [&]( PI index, TF m0, TF m1, TF v, int bad_cell, bool arc ) {
        if ( bad_cell )
            has_bad_cell = bad_cell;
        if ( arc )
            has_arc = arc;

        mm0 = min( mm0, m0 );

        newton_matrix_ldlt( index, index ) = m0;
        if ( m1 )
            newton_matrix_ldlt( index, index + 1 ) = m1;

        if ( m0 == 0 )
            has_bad_cell = 3;

        v -= sorted_dirac_masses[ index ];
        newton_vector[ index ] = v;
        newton_error += v * v;
    } );

    for( auto &v : newton_matrix_ldlt.values )
        if ( isnan( v ) )
            throw std::runtime_error( "nan mat (before fact)" );

    if ( has_bad_cell )
        return 1;

    if ( ! has_arc )
        newton_matrix_ldlt( 0, 0 ) *= 2;

    if ( int err = newton_matrix_ldlt.inplace_ldlt_decomposition() )
        return 1 + err;

    for( auto &v : newton_matrix_ldlt.values )
        if ( isnan( v ) )
            return 1;
        // throw std::runtime_error( "nan mat (after fact)" );

    if ( mm0 < 0 )
        return 2;
    // P( mm0 );

    return 0;
}

DTP void UTP::get_system_con( ConnectedCells &connected_cells, std::vector<Poly> &opt_weights, TF &max_a, TF &cell_error, int &has_bad_cell, const VF &sorted_dirac_weights ) const {
    using namespace std;

    //
    opt_weights.resize( nb_sorted_diracs() );
    connected_cells.clear();
    has_bad_cell = 0;
    cell_error = 0;
    max_a = 1;
    
    TF i0 = numeric_limits<TF>::lowest(); // prev cut position
    TF d0 = sorted_dirac_positions[ 0 ]; // prev dirac position
    TF w0 = sorted_dirac_weights[ 0 ]; // prev weight
    TF dp = 0; // prev prev dirac position
    TF wp = 0; // prev prev weight

    // PI bi_ind; // index of the bi cell
    for( PI n = 0; n < nb_sorted_diracs(); ++n ) {
        // negative weight ?
        if ( w0 < 0 ) {
            has_bad_cell = 1;
            return;
        }

        // next intersection
        const TF d1 = n + 1 < nb_sorted_diracs() ? sorted_dirac_positions[ n + 1 ] : numeric_limits<TF>::max();
        const TF w1 = n + 1 < nb_sorted_diracs() ? sorted_dirac_weights[ n + 1 ] : 0;
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
                connected_cells.push_back( { n, n } );
                // if ( sva == 0 ) {
                //     has_bad_cell = 4;
                //     return;
                // }
                const TF mde = 0; // sorted_dirac_masses[ n ] / original_density_values.size();

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
                connected_cells.push_back( { n, n } );
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

                const TF mde = 0; //sorted_dirac_masses[ n ] / original_density_values.size();
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

                const TF mde = 0; // sorted_dirac_masses[ n ] / original_density_values.size();
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

DTP TF UTP::best_r_for_ib( const std::vector<Poly> &polys, PI n, TF a ) const {
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

    const TF ter = target_max_mass_error * sorted_dirac_masses[ n ] * 1e-5;

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
                if ( sgn( err( x0 ) ) == sgn( err( x1 ) ) )
                    return -2;
                return dichotomy( err, ter, max( TF( 0 ), x0 ), x1 );
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
                if ( max_r > density->ptp_x() )
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
            if ( max_r > density->ptp_x() )
                return -7;
            max_r *= 2;
        }
        // P( __LINE__, max_r );
        return dichotomy( err, ter, TF( 0 ), max_r );
    }

    throw std::runtime_error( "TODO" );
}

DTP int UTP::get_weights_for( VF &new_dirac_weights, const ConnectedCells &connected_cells, const std::vector<Poly> &polys, TF a ) {
    new_dirac_weights.resize( nb_sorted_diracs() );
    for( auto inds : connected_cells ) {
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

DTP int UTP::newton_con( TF min_a ) {
    if ( nb_sorted_diracs() == 0 )
        return 0;

    // system content
    ConnectedCells new_connected_cells;
    ConnectedCells connected_cells;
    std::vector<Poly> new_opt_weights;
    std::vector<Poly> opt_weights;
    int has_bad_cell;
    TF new_cell_error;
    TF cell_error;
    TF max_a;

    // history
    // Vec<Vec<Vec<TF,2>>> bnds;

    // get the first system
    get_system_con( connected_cells, opt_weights, max_a, cell_error, has_bad_cell, sorted_dirac_weights );
    if ( has_bad_cell )
        return 100 + has_bad_cell; // throw std::runtime_error( "bad initialization" );
    // P( cell_error );
    if ( cell_error < target_max_mass_error / nb_sorted_diracs() )
        return 0;

    // iterate
    VF new_sorted_dirac_weights;
    for( PI num_iter = 0; ; ++num_iter ) {
        if ( num_iter == 1000 )
            return 2; // throw std::runtime_error( "too much iterations" );
        // bnds << normalized_cell_boundaries();

        for( TF a = max_a; ; a *= 0.8 ) {
            if ( a < min_a ) {
                // plot_bnds_evolution( bnds );
                return 3;
            }

            int esw = get_weights_for( new_sorted_dirac_weights, connected_cells, opt_weights, a );
            if ( esw == 0 ) {
                get_system_con( new_connected_cells, new_opt_weights, max_a, new_cell_error, has_bad_cell, new_sorted_dirac_weights );
                if ( has_bad_cell == 0 && new_cell_error < cell_error ) {
                    std::swap( new_sorted_dirac_weights, sorted_dirac_weights );
                    std::swap( new_connected_cells, connected_cells );
                    std::swap( new_opt_weights, opt_weights );
                    std::swap( new_cell_error, cell_error );
    
                    // P( cell_error, a );
                    P( sqrt( cell_error ) );
                    if ( sqrt( cell_error ) < target_max_mass_error / nb_sorted_diracs() ) {
                        if ( verbosity >= 2 && stream )
                            *stream << "    nb iteration con: " << num_iter << "\n";
                        return 0;
                    }
                    break;
                }
            }
        }
    }
}

DTP int UTP::newton_iterations( TF min_relax ) {
    using namespace std;

    int err = _make_newton_system();
    if ( err )
        return 1;

    TF prev_newton_error = newton_error;
    nb_newton_iterations = 0;
    for( nb_newton_iterations = 0; ; ++nb_newton_iterations ) {
        if ( nb_newton_iterations == 1000 )
            throw std::runtime_error( "max iter" );

        for( auto &v : newton_matrix_ldlt.values )
            if ( isnan( v ) )
                throw std::runtime_error( "nan mat" );

        for( auto &v : newton_vector )
            if ( isnan( v ) )
                throw runtime_error( "nan vec" );

        VF dir = newton_matrix_ldlt.solve_using_ldlt( newton_vector );
        for( auto &v : dir )
            if ( isnan( v ) )
                throw runtime_error( "nan dir" );

        VF w0 = sorted_dirac_weights;
        for( TF a = 1; ; a /= 2 ) {
            if ( a < min_relax ) {
                sorted_dirac_weights = w0;
                return 2;
            }

            for( PI i = 0; i < w0.size(); ++i )
                sorted_dirac_weights[ i ] = w0[ i ] - a * dir[ i ];

            int err = _make_newton_system();
            if ( err ) // || newton_error > prev_newton_error
                continue;

            prev_newton_error = newton_error;
            break;
        }
            
        if ( newton_error < std::numeric_limits<TF>::epsilon() )
            return 0;
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

DTP void UTP::solve( bool use_approx_for_ders ) {
    using namespace std;

    initialize_with_flat_density();

    density->set_flattening_ratio( 1 - 1e-6 );
    int err = newton_iterations( 1e-3 );
    if ( err )
        throw std::runtime_error( "bad init" );

    // for( auto b : cell_boundaries() ) {
    //     assert( b[ 0 ] < b[ 1 ] );
    //     assert( b[ 0 ] >= density->min_x() );
    //     assert( b[ 1 ] <= density->max_x() );
    // }
    // P( cell_boundaries()[ 0 ] );
    // P( cell_boundaries()[ 1 ] );
    // P( density->min_x() );

    density->set_flattening_ratio( 1e-3 );
    err = newton_iterations( 1e-6 );
    if ( err ) {
        // glot( linspace<TF>( density->min_x() - 0, density->max_x() + 0, 1000 ), 
        //     [&]( TF x ) { return density->value( x ); }
        // );
        plot();

        // P( cell_boundaries()[ 0 ] );
        // P( cell_boundaries()[ 1 ] );
        // P( density->min_x() );
        // P( err );
        throw std::runtime_error( "bad newton" );
    }
    return;

    //
    for( TF prev_flattening_ratio = density->current_flattening_ratio; prev_flattening_ratio; ) {
        // compute the derivatives
        VF w0 = sorted_dirac_weights;
        MF ders = der_weights_wrt_flat_ratio( 2, use_approx_for_ders );
        // for( PI i = 0; i < nb_sorted_diracs(); ++i )
        //     if ( isnan( ders[ 0 ][ i ] ) )
        //         throw std::runtime_error( "nan der 0" );
        // for( PI i = 0; i < nb_sorted_diracs(); ++i )
        //     if ( isnan( ders[ 1 ][ i ] ) )
        //         throw std::runtime_error( "nan der 1" );

        // go to the next `t`
        TF trat = 0.5;
        for( TF trial_flattening_ratio = 0; ; trial_flattening_ratio = ( 1 - trat ) * prev_flattening_ratio + trat * trial_flattening_ratio ) {
            //
            const TF a = trial_flattening_ratio - prev_flattening_ratio;
            density->set_flattening_ratio( trial_flattening_ratio );
            if ( abs( a ) < 1e-6 ) {
                throw std::runtime_error( "low a" );
            }


            // test with the derivatives
            for( PI i = 0; i < nb_sorted_diracs(); ++i )
                sorted_dirac_weights[ i ] = w0[ i ] + ders[ 0 ][ i ] * a + ders[ 1 ][ i ] * pow( a, 2 ) / 2; // + ders[ 2 ][ i ] * pow( a, 3 ) / 6;
            err = newton_iterations( 1e-2 );
            if ( err == 0 ) {
                prev_flattening_ratio = trial_flattening_ratio;
                if ( verbosity >= 2 && stream )
                    *stream << "  flattening_ratio: " << prev_flattening_ratio << " nb_newton_iterations: " << nb_newton_iterations << " (2)\n";
                break;
            }

            // // test with the derivatives
            // for( PI i = 0; i < nb_sorted_diracs(); ++i )
            //     sorted_dirac_weights[ i ] = w0[ i ] + ders[ 0 ][ i ] * a;
            // err = newton_iterations( 1e-2 );
            // if ( err == 0 ) {
            //     prev_flattening_ratio = trial_flattening_ratio;
            //     if ( verbosity >= 2 && stream )
            //         *stream << "  flattening_ratio: " << prev_flattening_ratio << " nb_newton_iterations: " << nb_newton_iterations << " (1)\n";
            //     break;
            // }

            // test without the derivatives
            for( PI i = 0; i < nb_sorted_diracs(); ++i )
                sorted_dirac_weights[ i ] = w0[ i ];
            err = newton_iterations( 13e-3 );
            if ( err == 0 ) {
                prev_flattening_ratio = trial_flattening_ratio;
                if ( verbosity >= 2 && stream )
                    *stream << "  flattening_ratio: " << prev_flattening_ratio << " nb_newton_iterations: " << nb_newton_iterations << " (0)\n";
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
