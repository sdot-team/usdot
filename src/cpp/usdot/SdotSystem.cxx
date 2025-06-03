#pragma once

#include "SdotSystem.h"
#include <stdexcept>
#include <iostream>
#include <numeric>
#include <fstream>
#include <limits>
#include <cmath>

namespace usdot {

#define DTP template<class TF,class Density>
#define UTP SdotSystem<TF,Density>

DTP UTP::SdotSystem( Density *density, const VF &dirac_positions, TF global_mass_ratio, TF sep_equ_coords ) : density( density ), global_mass_ratio( global_mass_ratio ) {
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
    sorted_dirac_weights.resize( nb_sorted_diracs() );
    std::fill( sorted_dirac_weights.begin(), sorted_dirac_weights.end(), 0.1 );
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

DTP int UTP::for_each_interval( auto &&func ) const {
    using namespace std;

    const PI nc = sorted_dirac_positions.size();
    if ( nc == 0 )
        return 0;

    //
    int error = 0;
    auto call_with_ball_cut = [&]( PI num_sorted_dirac, TF c0, TF c1 ) {
        if ( c0 >= c1 ) {
            error |= 2;
            return;
        }

        const TF we = sorted_dirac_weights[ num_sorted_dirac ];
        if ( we <= 0 ) {
            error |= 1;
            return;
        }

        const TF xc = sorted_dirac_positions[ num_sorted_dirac ];
        const TF rb = sqrt( we );
        const TF b0 = xc - rb;
        const TF b1 = xc + rb;

        // void intersection ?
        if ( c0 >= b1 || c1 <= b0 ) {
            error |= 4;
            return;
        }

        // epsilon cut ?
        if ( epsilon < rb ) {
            const TF re = sqrt( we - epsilon * epsilon );
            const TF e0 = xc - re;
            const TF e1 = xc + re;

            // start with the left disk
            if ( c0 > e1 ) {
                // cut between eps and end disk ?
                if ( c1 < b1 )
                    return func( num_sorted_dirac, Interval{ BndType::DiskToDisk, BndType::DiskToDisk, c0, c1 } );

                // -> end with the disk
                return func( num_sorted_dirac, Interval{ BndType::DiskToDisk, BndType::DiskToVoid, c0, b1 } );
            }

            // start with the cut ?
            if ( c0 > e0 ) {            
                // end with the cut ?
                if ( c1 < e1 )
                    return func( num_sorted_dirac, Interval{ BndType::EpsToEps, BndType::EpsToEps, c0, c1 } );

                // cut between eps and end disk ?
                if ( c1 < b1 ) {
                    func( num_sorted_dirac, Interval{ BndType::EpsToEps, BndType::EpsToDisk, c0, e1 } );
                    func( num_sorted_dirac, Interval{ BndType::EpsToDisk, BndType::DiskToDisk, e1, c1 } );
                    return;
                }

                // -> end with the disk
                func( num_sorted_dirac, Interval{ BndType::EpsToEps, BndType::EpsToDisk, c0, e1 } );
                func( num_sorted_dirac, Interval{ BndType::EpsToDisk, BndType::DiskToVoid, e1, b1 } );
                return;
            }

            // in the middle of b0 and e0 ?
            if ( c0 > b0 ) {
                //
                if ( c1 < e0 )
                    return func( num_sorted_dirac, Interval{ BndType::DiskToDisk, BndType::DiskToEps, c0, c1 } );

                // end with the cut ?
                if ( c1 < e1 ) {
                    func( num_sorted_dirac, Interval{ BndType::DiskToDisk, BndType::DiskToEps, c0, e0 } );
                    func( num_sorted_dirac, Interval{ BndType::DiskToEps, BndType::EpsToEps, e0, c1 } );
                    return;
                }

                // cut between eps and end disk ?
                if ( c1 < b1 ) {
                    func( num_sorted_dirac, Interval{ BndType::DiskToDisk, BndType::DiskToEps, c0, e0 } );
                    func( num_sorted_dirac, Interval{ BndType::DiskToEps, BndType::EpsToDisk, e0, e1 } );
                    func( num_sorted_dirac, Interval{ BndType::EpsToDisk, BndType::DiskToDisk, e1, c1 } );
                    return;
                }

                // -> end with the disk
                func( num_sorted_dirac, Interval{ BndType::DiskToDisk, BndType::DiskToEps, c0, e0 } );
                func( num_sorted_dirac, Interval{ BndType::DiskToEps, BndType::EpsToDisk, e0, e1 } );
                func( num_sorted_dirac, Interval{ BndType::EpsToDisk, BndType::DiskToVoid, e1, b1 } );
                return;
            }

            // -> start with the disk
            // end with the left disk ?
            if ( c1 < e0 ) {
                func( num_sorted_dirac, Interval{ BndType::VoidToDisk, BndType::DiskToDisk, b0, c1 } );
                return;
            }

            // end with the cut ?
            if ( c1 < e1 ) {
                func( num_sorted_dirac, Interval{ BndType::VoidToDisk, BndType::DiskToEps, b0, e0 } );
                func( num_sorted_dirac, Interval{ BndType::DiskToEps, BndType::EpsToEps, e0, c1 } );
                return;
            }

            // cut between eps and end disk ?
            if ( c1 < b1 ) {
                func( num_sorted_dirac, Interval{ BndType::VoidToDisk, BndType::DiskToEps, b0, e0 } );
                func( num_sorted_dirac, Interval{ BndType::DiskToEps, BndType::EpsToDisk, e0, e1 } );
                func( num_sorted_dirac, Interval{ BndType::EpsToDisk, BndType::DiskToDisk, e1, c1 } );
                return;
            }

            // -> end with the disk
            func( num_sorted_dirac, Interval{ BndType::VoidToDisk, BndType::DiskToEps, b0, e0 } );
            func( num_sorted_dirac, Interval{ BndType::DiskToEps, BndType::EpsToDisk, e0, e1 } );
            func( num_sorted_dirac, Interval{ BndType::EpsToDisk, BndType::DiskToVoid, e1, b1 } );
            return;
        } 

        // start with the disk ?
        if ( c0 < b0 ) {
            // end with the disk ?
            if ( c1 > b1 )
                return func( num_sorted_dirac, Interval{ BndType::VoidToDisk, BndType::DiskToVoid, b0, b1 } );
            return func( num_sorted_dirac, Interval{ BndType::VoidToDisk, BndType::DiskToDisk, b0, c1 } );
        }

        // end with the disk ?
        if ( c1 > b1 )
            return func( num_sorted_dirac, Interval{ BndType::DiskToDisk, BndType::DiskToVoid, c0, b1 } );
        return func( num_sorted_dirac, Interval{ BndType::DiskToDisk, BndType::DiskToDisk, c0, c1 } );
    };

    //
    TF c0 = std::numeric_limits<TF>::lowest();
    TF d0 = sorted_dirac_positions[ 0 ];
    TF w0 = sorted_dirac_weights[ 0 ];
    for( PI i = 1; i < nc; ++i ) {
        const TF d1 = sorted_dirac_positions[ i ];
        const TF w1 = sorted_dirac_weights[ i ];

        const TF c1 = ( d1 + d0 - ( w1 - w0 ) / ( d1 - d0 ) ) / 2;
        if ( std::isnan( c1 ) )
            return 8;

        call_with_ball_cut( i - 1, c0, c1 );

        c0 = c1;
        d0 = d1;
        w0 = w1;
    }

    call_with_ball_cut( nc - 1, c0, std::numeric_limits<TF>::max() );

    return error;
}

DTP T_T int UTP::for_each_cell_polyline( T &&func, PI nb_divs ) const {
    using namespace std;

    PI prev_num_sorted_cell = -1;
    VF xs;
    VF ys;
    int res = for_each_interval( [&]( PI num_sorted_cell, Interval<TF> interval ) {
        if ( num_sorted_cell != prev_num_sorted_cell ) {
            if ( xs.size() ) {
                func( prev_num_sorted_cell, xs, ys );
                xs.clear();
                ys.clear();
            }
            prev_num_sorted_cell = num_sorted_cell;
        }

        const TF d = sorted_dirac_positions[ num_sorted_cell ];
        const TF w = sorted_dirac_weights[ num_sorted_cell ];
        const TF r = sqrt( w );

        auto app = [&]( TF x, TF y ) {
            if ( xs.empty() || xs.back() != x || ys.back() != y ) {
                xs.push_back( x );
                ys.push_back( y );
            }
        };

        if ( interval.is_disk() ) {
            const TF a0 = acos( clamp( TF( ( interval.x0 - d ) / r ), TF( -1 ), TF( +1 ) ) );
            const TF a1 = acos( clamp( TF( ( interval.x1 - d ) / r ), TF( -1 ), TF( +1 ) ) );
            const PI na = PI( nb_divs * abs( a1 - a0 ) / acos( TF( -1 ) ) + 1 );
            for( PI ia = 0; ia <= na; ++ia ) {
                app( 
                    d + r * cos( a0 + ( a1 - a0 ) * ia / na ),
                    r * sin( a0 + ( a1 - a0 ) * ia / na )
                );
            }
        } else {
            app( interval.x0, epsilon );
            app( interval.x1, epsilon );
        }
    } );

    if ( prev_num_sorted_cell != -1 && xs.size() )
        func( prev_num_sorted_cell, xs, ys );

    return res;
}

#ifdef TL_DISPLAYER_IS_DEFINED
DTP void UTP::display( Displayer &ds ) const {
    ds.start_array();
    for_each_interval( [&]( PI num_sorted_cell, Interval<TF> interval ) {
        ds << interval;
    } );
    ds.end_array();
}
#endif

DTP typename UTP::VF UTP::sorted_cell_masses() const {
    using namespace std;
    VF res( nb_sorted_diracs(), 0 );
    for_each_interval( [&]( PI num_sorted_cell, Interval<TF> interval ) {
        if ( interval.is_disk() ) {
            const TF c = sorted_dirac_positions[ num_sorted_cell ];
            const TF w = sorted_dirac_weights[ num_sorted_cell ];
            const TF r = sqrt( w );
            res[ num_sorted_cell ] += density->disk_integral( interval.x0, interval.x1, c, r );
        } else
            res[ num_sorted_cell ] += 2 * epsilon * density->integral( interval.x0, interval.x1 );
    } );
    for( TF &v : res )
        v /= 2 * epsilon;
    return res;
}

DTP void UTP::plot( Str filename ) const {
    std::ofstream fs( filename );

    fs << "from matplotlib import pyplot\n";

    for_each_cell_polyline( [&]( PI num_sorted_cell, const VF &xs, const VF &ys ) {
        if ( xs.empty() )
            return;

        fs << "pyplot.plot( [ ";

        for( PI i = 0; i < xs.size(); ++i )
            fs << xs[ i ] << ",";
        for( PI i = xs.size(); i--; )
            fs << xs[ i ] << ",";
        fs << xs[ 0 ] << ",";
        
        fs << "], [\n";
        for( PI i = 0; i < ys.size(); ++i )
            fs << ys[ i ] << ",";
        for( PI i = ys.size(); i--; )
            fs << - ys[ i ] << ",";
        fs << ys[ 0 ] << ",";
        
        fs << "] )\n";
    } );

    // // density
    // density->plot( fs );

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

DTP typename UTP::VF UTP::sorted_cell_masses_ap( PI ni ) const {
    VF res( nb_sorted_diracs() );

    for_each_cell_polyline( [&]( PI num_sorted_cell, const VF &xs, const VF &ys ) {
        TF tot = 0;
        for( PI i = 1; i < xs.size(); ++i ) {
            const TF x0 = xs[ i - 1 ];
            const TF x1 = xs[ i - 0 ];
            const TF y0 = ys[ i - 1 ];
            const TF y1 = ys[ i - 0 ];

            tot += density->integral( x0, x1 ) * ( y1 + y0 );
        }
        res[ num_sorted_cell ] = tot;
    }, ni );

    return res;
}

DTP void UTP::check_newton_system( TF eps ) {
    TridiagonalSymmetricMatrix<TF> M;
    TF max_error_ratio;
    VF V;
    int err = get_sorted_newton_system( M, V, max_error_ratio );

    TridiagonalSymmetricMatrix<TF> N;
    TF max_error_ratio_a;
    VF U;
    get_sorted_newton_system_ap( N, U, max_error_ratio_a, eps / 1e10 );

    if ( err )
        return;
    // P( V, U );
    // P( M, N );
    // P( max_error_ratio_a, max_error_ratio );
    // P( U, V );
    // P( M, N );
    assert( abs( max_error_ratio_a - max_error_ratio ) < eps );
    for( PI i = 0; i < V.size(); ++i )
        assert( abs( U[ i ] - V[ i ] ) < eps );
    for( PI i = 0; i < M.values.size(); ++i ) {
        if ( abs( N.values[ i ] - M.values[ i ] ) < eps )
            continue;
        std::cout << "check_newton_system: error at index " << i << ": " << N.values[ i ] - M.values[ i ] << std::endl;
        assert( 0 );
    }
    std::cout << "check_newton_system: ok" << std::endl;
}

DTP int UTP::get_sorted_newton_system_ap( TridiagonalSymmetricMatrix<TF> &M, VF &V, TF &max_error_ratio, TF eps ) {
    using namespace std;
    int error = 0;

    V = sorted_cell_masses();

    M.resize( nb_sorted_diracs() );
    for( PI n = 0; n < nb_sorted_diracs(); ++n ) {
        TF &ref_weight = sorted_dirac_weights[ n ];
        TF old_weight = ref_weight;
        ref_weight += eps;

        VF A = sorted_cell_masses();
        for( PI m = max( n, 1ul ) - 1; m <= n; ++m )
            M( m, n ) = ( A[ m ] - V[ m ] ) / eps;

        ref_weight = old_weight;
    }

    for( PI i = 0; i < V.size(); ++i )
        V[ i ] = sorted_dirac_masses[ i ] - V[ i ];

    max_error_ratio = 0;
    for( PI i = 0; i < V.size(); ++i ) {
        const TF err = abs( V[ i ] ) / sorted_dirac_masses[ i ];
        max_error_ratio = max( max_error_ratio, err );

        // void cell
        if ( V[ i ] == sorted_dirac_masses[ i ] )
            error = 4;
    }

    return error;
}

DTP int UTP::get_sorted_newton_system( TridiagonalSymmetricMatrix<TF> &M, VF &V, TF &max_error_ratio ) const {
    using namespace std;
    
    M.resize( nb_sorted_diracs() );
    V.resize( nb_sorted_diracs() );
    M.clear_values();
   
    for( PI i = 0; i < nb_sorted_diracs(); ++i )
        V[ i ] = sorted_dirac_masses[ i ];

    TF mass_arcs = 0;
    int error = for_each_interval( [&]( PI num_cell, Interval<TF> interval ) {
        const TF c = sorted_dirac_positions[ num_cell ];
        const TF w = sorted_dirac_weights[ num_cell ];
        const TF r = sqrt( max( TF( 0 ), w ) );
        assert( interval.x0 <= interval.x1 );

        const TF rdist = num_cell + 1 < nb_sorted_diracs() ? sorted_dirac_positions[ num_cell + 1 ] - c : numeric_limits<TF>::max();
        const TF ldist = num_cell ? c - sorted_dirac_positions[ num_cell - 1 ] : numeric_limits<TF>::max();
        
        TF mass = 0, dmass_dL = 0, dmass_dM = 0;

        if ( interval.is_disk() ) {
            TF loc_dmass_dM = 0;
            TF loc_mass = 0;
            density->disk_integral( interval.x0, interval.x1, c, r, loc_mass, loc_dmass_dM );
            loc_dmass_dM /= 2 * epsilon;
            loc_mass /= 2 * epsilon;

            // P( loc_dmass_dM );
            if ( loc_dmass_dM < 0 ) {
                loc_dmass_dM = 0;
                error |= 16;
            }

            dmass_dM += loc_dmass_dM;
            mass += loc_mass;

            mass_arcs += loc_mass;
        } else
            mass += density->integral( interval.x0, interval.x1 );

        const TF v0 = max( density->value( interval.x0 ), 1e-2 * density->max_value );
        if ( interval.t0 == BndType::DiskToDisk ) {
            const TF h = v0 * sqrt( pow( r, 2 ) - pow( interval.x0 - c, 2 ) );
            const TF d = h / ( ldist * 2 * epsilon );
            dmass_dM += d;
            dmass_dL -= d;
        } else if ( interval.t0 == BndType::EpsToEps ) {
            const TF h = v0;
            const TF d = h / ldist / 2;
            dmass_dM += d;
            dmass_dL -= d;
        }

        const TF v1 = max( density->value( interval.x1 ), 1e-2 * density->max_value );
        if ( interval.t1 == BndType::DiskToDisk ) {
            const TF h = v1 * sqrt( pow( r, 2 ) - pow( interval.x1 - c, 2 ) );
            const TF d = h / ( rdist * 2 * epsilon );
            dmass_dM += d;
        } else if ( interval.t1 == BndType::EpsToEps ) {
            const TF h = v1;
            const TF d = h / rdist / 2;
            dmass_dM += d;
        }

        if ( dmass_dL ) M( num_cell, num_cell - 1 ) += dmass_dL;
        M( num_cell, num_cell ) += dmass_dM;
        V[ num_cell ] -= mass;
    } );

    if ( mass_arcs == 0 )
        error = 3;

    // 
    max_error_ratio = 0;
    for( PI i = 0; i < V.size(); ++i ) {
        const TF err = abs( V[ i ] ) / sorted_dirac_masses[ i ];
        max_error_ratio = max( max_error_ratio, err );

        // void cell
        if ( V[ i ] == sorted_dirac_masses[ i ] )
            error = 4;
    }

    return error;
}

DTP void UTP::initialize_weights() {
    using namespace std;

    // agglomerate = cells glued together
    struct Agglomerate {
        TF beg_x;
        TF end_x;
        TF len_x;
        PI beg_n;
        PI end_n;
    };
    std::vector<Agglomerate> aggs;
    aggs.reserve( nb_sorted_diracs() / 2 );
    
    // make the agglomerates
    const TF base_len = density->ptp_x() * global_mass_ratio / nb_original_diracs();
    for( PI i = 0; i < nb_sorted_diracs(); ++i ) {
        const TF l = base_len * ( sorted_dirac_num_offsets[ i + 1 ] - sorted_dirac_num_offsets[ i ] );
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
    for( auto &v : res )
        v[ 0 ] = std::numeric_limits<TF>::max(), v[ 1 ] = std::numeric_limits<TF>::lowest();
    for_each_interval( [&]( PI num_cell, Interval<TF> interval ) {
        res[ num_cell ][ 0 ] = std::min( res[ num_cell ][ 0 ], interval.x0 );
        res[ num_cell ][ 1 ] = std::max( res[ num_cell ][ 1 ], interval.x1 );
    } );
    return res;
}

DTP int UTP::newton_iterations( const TF min_relax ) {
    // first newton system
    TridiagonalSymmetricMatrix<TF> M;
    TF max_error_ratio;
    VF V;
    int err = get_sorted_newton_system( M, V, max_error_ratio );
    if ( err )
        return 1;
    
    for( nb_newton_iterations = 0; ; ++nb_newton_iterations ) {
        if ( nb_newton_iterations == 500 )
            throw std::runtime_error( "too many newton iterations" );
        
        // on a solution ?
        if ( max_error_ratio < target_max_error_ratio )
            return 0;
        
        assert( no_nan( M.values ) );
        assert( no_nan( V ) );

        // get a direction
        M.inplace_ldlt_decomposition();
        VF dir = M.solve_using_ldlt( V );

        // plot();
        assert( no_nan( M.values ) );
        assert( no_nan( dir ) );

        // find the relaxation coefficient
        VF w0 = sorted_dirac_weights;
        for( TF relax = 1; ; relax /= 2 ) {
            for( PI i = 0; i < w0.size(); ++i )
                sorted_dirac_weights[ i ] = w0[ i ] + relax * dir[ i ];

            // check_newton_system( 1e-10 );

            err = get_sorted_newton_system( M, V, max_error_ratio );
            if ( err == 0 ) {
                bnds.push_back( sorted_cell_boundaries() );

                // if ( max_error_ratio > 20 ) {
                //     P( M.diag() );
                //     plot();
                //     TODO;
                // }
                // mass_history.push_back( sorted_cell_masses() );
                if ( nb_newton_iterations == 39 ) {
                    plot_bnds_evolution( bnds );
                    // P( sorted_cell_masses() );
                    // TF s = 0;
                    // for( TF v : sorted_cell_masses() )
                    //     s += v;
                    // P( s );
                    // auto &a = mass_history.back();
                    // auto &b = mass_history[ mass_history.size() - 2 ];
                    // for( PI i = 0; i < nb_sorted_diracs(); ++i )
                    //     P( a[ i ], b[ i ], a[ i ] - b[ i ] );
                    assert( 0 );
                }
                if ( verbosity >= 2 )   
                    std::cout << "  relax: " << relax << " nb_newton_iterations: " << nb_newton_iterations << " max_error_ratio: " << max_error_ratio << "\n";
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
    solve_with_diffusion();
    // initialize_weights();

    // int err = newton_iterations( 1e-5, 10 );
    // P( err );
}

DTP void UTP::solve_with_diffusion() {
    initialize_weights();
    int err = newton_iterations( 1e-5 );
    if ( err )
        throw std::runtime_error( "bad init" );

    // // flat density
    // density->set_flattening_ratio( 1 );

    // int err = newton_iterations( 1e-5 );
    // if ( err )
    //     throw std::runtime_error( "bad init" );

    // //
    // for( TF prev_flattening_ratio = density->current_flattening_ratio; prev_flattening_ratio; ) {
    //     // go to the next `t`
    //     VF w0 = sorted_dirac_weights;
    //     const TF t_ratio = 0.5;
    //     for( TF trial_flattening_ratio = 0; ; trial_flattening_ratio = ( 1 - t_ratio ) * prev_flattening_ratio + t_ratio * trial_flattening_ratio ) {
    //         // change the density
    //         if ( verbosity >= 2 )
    //             std::cout << "try flattening_ratio: " << trial_flattening_ratio << "\n";
    //         density->set_flattening_ratio( trial_flattening_ratio );

    //         // test delta t
    //         const TF a = trial_flattening_ratio - prev_flattening_ratio;
    //         if ( abs( a ) < 1e-10 ) {
    //             // plot();
    //             throw std::runtime_error( "low a" );
    //         }

    //         // test newton iterations
    //         for( PI i = 0; i < nb_sorted_diracs(); ++i )
    //             sorted_dirac_weights[ i ] = w0[ i ];
    //         int err = newton_iterations( 1e-4 );
    //         if ( err == 0 ) {
    //             prev_flattening_ratio = trial_flattening_ratio;
    //             if ( verbosity >= 2 )
    //                 std::cout << "flattening_ratio: " << prev_flattening_ratio << " nb_newton_iterations: " << nb_newton_iterations << "\n";
    //             break;
    //         }
    //     }
    // }
}

// DTP VF UTP::barycenters_ap( const Density<TF> &density, bool sorted_nums, PI ni ) const {
//     VF res( FromSize(), nb_cells() );

//     auto diter = density->iterator();
//     for_each_cell( [&]( PI num_cell, Interval<TF> interval ) {
//         TF pint = 0, mass = 0;
//         for_each_sub_interval( *diter, interval, num_cell, [&]( Interval<TF> interval ) {
//             // diter.barycenter( pint, area, interval.x0, interval.x1 );
//             for( PI i = 0; i < ni; ++i ) {
//                 const TF x = interval.x0 + ( interval.x1 - interval.x0 ) * ( i + 0.5 ) / ni;
//                 const TF a = ( interval.x1 - interval.x0 ) * diter->value( x );
//                 pint += a * x;
//                 mass += a;
//             }            
//         } );
//         res[ sorted_nums ? num_cell : sorted_dirac_nums[ num_cell ] ] = mass ? pint / mass : 0;
//     } );

//     return res;
// }

// DTP VF UTP::barycenters( const Density<TF> &density, bool sorted_nums ) const {
//     VF res( FromSize(), nb_cells() );

//     auto diter = density.iterator();
//     for_each_cell( [&]( PI num_cell, Interval<TF> interval ) {
//         TF pint = 0, area = 0;

//         for_each_sub_interval( *diter, interval, num_cell, [&]( Interval<TF> interval ) {
//             diter->barycenter( pint, area, interval.x0, interval.x1 );
//         } );

//         res[ sorted_nums ? num_cell : sorted_dirac_nums[ num_cell ] ] = area ? pint / area : 0;
//     } );

//     return res;
// }

// DTP VF UTP::masses( const Density<TF> &density, bool sorted_nums ) const {
//     VF res( FromSize(), nb_cells() );

//     auto diter = density.iterator();
//     for_each_cell( [&]( PI num_cell, Interval<TF> interval ) {
//         TF area = 0;

//         for_each_sub_interval( *diter, interval, num_cell, [&]( Interval<TF> interval ) {
//             area += diter->integral( interval.x0, interval.x1 );
//         } );

//         res[ sorted_nums ? num_cell : sorted_dirac_nums[ num_cell ] ] = area;
//     } );

//     return res;
// }

// DTP void UTP::set_weights( const VF &weights, bool sorted_nums ) const {
//     for( PI i = 0; i < nb_cells(); ++i )
//         sorted_dirac_weights[ i ] = weights[ sorted_nums ? i : sorted_dirac_nums[ i ] ];
// }

// DTP VF UTP::get_weights( bool sorted_nums ) const {
//     VF res( FromSize(), nb_cells() );
//     for( PI i = 0; i < nb_cells(); ++i )
//         res[ sorted_nums ? i : sorted_dirac_nums[ i ] ] = sorted_dirac_weights[ i ];
//     return res;
// }

#undef DTP
#undef UTP

} // namespace usdot
