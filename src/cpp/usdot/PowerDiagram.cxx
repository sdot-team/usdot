#pragma once

#include "PowerDiagram.h"
#include <limits>

namespace usdot {

#define DTP template<class TF>
#define UTP PowerDiagram<TF>

DTP UTP::PowerDiagram( const VF &seed_coords, TF sep_equ_coords ) {
    sorted_seed_nums = { FromSizeAndFunctionOnIndex(), seed_coords.size(), []( auto i ) { return i; } };
    std::sort( sorted_seed_nums.begin(), sorted_seed_nums.end(), [&]( PI a, PI b ) {
        return seed_coords[ a ] < seed_coords[ b ];
    } );

    sorted_seed_weights.resize( seed_coords.size() );
    sorted_seed_coords.resize( seed_coords.size() );
    for( PI i = 0; i < seed_coords.size(); ++i ) {
        sorted_seed_weights[ i ] = seed_weights[ sorted_seed_nums[ i ] ];
        sorted_seed_coords[ i ] = seed_coords[ sorted_seed_nums[ i ] ];
    }

    // for( PI i = 1; i < sorted_seed_coords.size(); ++i ) {
    //     PI b = i - 1;
    //     while ( i < sorted_seed_coords.size() && sorted_seed_coords[ i ] == sorted_seed_coords[ b ] )
    //         ++i;
    //     if ( i > b + 1 ) {
    //         const TF pe = i < sorted_seed_coords.size() ? sorted_seed_coords[ i ]  : sorted_seed_coords[ i - 1 ];
    //         const TF pb = b ? sorted_seed_coords[ b - 1 ]  : sorted_seed_coords[ b ];
    //         for( PI j = b; j < i; ++j ) {
    //             const TF t = 0.5 + sep_equ_coords / ( i - b ) * ( TF( j - b ) / ( i - 1 - b ) - 0.5 );
    //             sorted_seed_coords[ j ] = pb + ( pe - pb ) * t;
    //         }
    //     }
    // }

    // 
    if ( const PI s = sorted_seed_coords.size() ) {
        const TF de = ( sorted_seed_coords.back() - sorted_seed_coords.front() ) / ( 1e2 * s );
        TF px = sorted_seed_coords[ 0 ];
        for( PI i = 1; i < sorted_seed_coords.size(); ++i ) {
            if ( sorted_seed_coords[ i ] < px + de )
                sorted_seed_coords[ i ] = px + de;
            px = sorted_seed_coords[ i ];
        }
    }
}

DTP void UTP::for_each_cell( auto &&func ) const {
    using namespace std;
    using Bt = BndType;

    const PI nc = sorted_seed_coords.size();
    if ( nc == 0 )
        return;

    //
    auto call_with_ball_cut = [&]( PI num_seed, Bt t0, Bt t1, TF x0, TF x1 ) {
        if ( allow_ball_cut ) {
            const TF ra = sqrt( max( 0, sorted_seed_weights[ num_seed ] ) );
            const TF xc = sorted_seed_coords[ num_seed ];
            const TF bl = xc - ra;
            const TF br = xc + ra;

            if ( x0 < bl ) {
                t0 = Bt::Ball;
                x0 = bl;
            } else if ( x0 > br ) {
                t0 = Bt::Ball;
                x0 = br;
            }

            if ( x1 < bl ) {
                t1 = Bt::Ball;
                x1 = bl;
            } else if ( x1 > br ) {
                t1 = Bt::Ball;
                x1 = br;
            }
        }

        func( num_seed, Interval{ t0, t1, x0, x1 } );
    };

    //
    TF x0 = std::numeric_limits<TF>::lowest();
    Bt t0 = Bt::Inf;
    for( PI i = 1; i < nc; ++i ) {
        const TF w0 = sorted_seed_weights[ i - 1 ];
        const TF w1 = sorted_seed_weights[ i - 0 ];
        const TF d0 = sorted_seed_coords[ i - 1 ];
        const TF d1 = sorted_seed_coords[ i - 0 ];

        TF x1 = ( d1 + d0 - ( w1 - w0 ) / ( d1 - d0 ) ) / 2;
        Bt t1 = Bt::Cell;

        call_with_ball_cut( i - 1, t0, t1, x0, x1 );

        x0 = x1;
        t0 = t1;
    }

    call_with_ball_cut( nc - 1, t0, Bt::Inf, x0, std::numeric_limits<TF>::max() );
}

DTP void UTP::__for_each_sub_interval( DensityIterator<TF> &density_iterator, const Interval<TF> &cell_interval, PI cell_num, auto &&func ) const {
    const BndType c_t0 = cell_interval.t0;
    const BndType c_t1 = cell_interval.t1;
    const TF c_x0 = cell_interval.x0;
    const TF c_x1 = cell_interval.x1;

    // move density_iterator at the beginning
    while ( density_iterator.x0 > c_x0 )
        if ( ! density_iterator.move_backward() )
            break;

    while ( density_iterator.x1 <= c_x0 )
        if ( ! density_iterator.move_forward() )
            return;

    // if no possible intersection, early return 
    if ( density_iterator.x0 >= c_x1 )
        return;

    // if density_iterator contains c_x0
    if ( density_iterator.x0 <= c_x0 ) {
        // if density_iterator contains also c_x1, we only have to make one call
        if ( density_iterator.x1 >= c_x1 ) {
            func( Interval<TF>( c_t0, c_t1, c_x0, c_x1 ) );
            return;
        }

        // else, call with interval ending at density_iterator.x1
        func( Interval<TF>( c_t0, BndType::Density, c_x0, density_iterator.x1 ) );
        if ( ! density_iterator.move_forward() )
            return;
    }

    // for each density intervals fully contained in the cell interval
    while ( density_iterator.x1 < c_x1 ) {
        func( Interval<TF>( BndType::Density, BndType::Density, density_iterator.x0, density_iterator.x1 ) );
        if ( ! density_iterator.move_forward() )
            return;
    }

    // if still common data
    if ( density_iterator.x0 < c_x1 )
        func( Interval<TF>( BndType::Density, c_t1, density_iterator.x0, c_x1 ) );
}

DTP void UTP::_for_each_sub_interval( DensityIterator<TF> &density_iterator, const Interval<TF> &cell_interval, PI cell_num, auto &&func ) const {
    using namespace std;
    
    if ( allow_ball_cut && epsilon && ( cell_interval.t0 == BndType::Ball || cell_interval.t1 == BndType::Ball ) ) {
        // const TF radius = sqrt( max( 0, sorted_seed_weights[ cell_num ] ) );
        // const TF center = sorted_seed_coords[ cell_num ];
        TODO;
        return;
    }

    __for_each_sub_interval( density_iterator, cell_interval, cell_num, FORWARD( func ) );
}

DTP void UTP::for_each_sub_interval( DensityIterator<TF> &density_iterator, const Interval<TF> &cell_interval, PI cell_num, auto &&func ) const {
    // need to swap the cell boundaries ?
    if ( cell_interval.x0 > cell_interval.x1 ) {
        Interval<TF> si( cell_interval.t1, cell_interval.t0, cell_interval.x1, cell_interval.x0 );
        return _for_each_sub_interval( density_iterator, si, cell_num, [&]( Interval<TF> interval ) {
            std::swap( interval.t0, interval.t1 );
            std::swap( interval.x0, interval.x1 );
            func( interval );
        } );
    }

    //
    return _for_each_sub_interval( density_iterator, cell_interval, cell_num, FORWARD( func ) );
}

DTP void UTP::get_newton_system_ap( SymmetricBandMatrix<TF> &M, VF &V, PI &nb_arcs, const Density<TF> &density, TF eps ) {
    using std::max;
    using std::min;

    VF R = masses( density );

    for( PI n = 0; n < nb_cells(); ++n ) {
        TF &ref_weight = sorted_seed_weights[ n ];
        TF old_weight = ref_weight;
        ref_weight += eps;

        VF A = masses( density );
        for( PI m = max( n, 1ul ) - 1; m <= n; ++m )
            M( m, n ) += ( A[ m ] - R[ m ] ) / eps;

        ref_weight = old_weight;
    }

    V -= R;
}

DTP void UTP::get_newton_system( SymmetricBandMatrix<TF> &M, VF &V, PI &nb_arcs, const Density<TF> &density, TF coeff ) const {
    RcPtr<DensityIterator<TF>> diter = density.iterator();
    for_each_cell( [&]( PI num_cell, Interval<TF> interval ) {
        const TF rdist = num_cell + 1 < nb_cells() ? sorted_seed_coords[ num_cell + 1 ] - sorted_seed_coords[ num_cell + 0 ] : 0;
        const TF ldist = num_cell ? sorted_seed_coords[ num_cell - 0 ] - sorted_seed_coords[ num_cell - 1 ] : 0;
        
        TF mass = 0, dmass_dL = 0, dmass_dM = 0;
        for_each_sub_interval( *diter, interval, num_cell, [&]( Interval<TF> interval ) {
            switch ( interval.t0 ) {
                case BndType::Cell: {
                    const TF value = diter->value( interval.x0 ) / ldist / 2;
                    dmass_dL -= value;
                    dmass_dM += value;
                    break;
                }
                case BndType::Ball: {
                    const TF value = diter->value( interval.x0 ) / sqrt( sorted_seed_weights[ num_cell + 0 ] ) / 2;
                    dmass_dM += value;
                    ++nb_arcs;
                    break;
                }
                default:
                    break;
            }

            switch ( interval.t1 ) {
                case BndType::Cell: {
                    const TF value = diter->value( interval.x1 ) / rdist / 2;
                    dmass_dM += value;
                    break;
                }
                case BndType::Ball: {
                    const TF value = diter->value( interval.x0 ) / sqrt( sorted_seed_weights[ num_cell + 0 ] ) / 2;
                    dmass_dM += value;
                    ++nb_arcs;
                    break;
                }
                default:
                    break;
            }

            mass += diter->integral( interval.x0, interval.x1 );
        } );

        if ( dmass_dL ) M( num_cell, num_cell - 1 ) += coeff * dmass_dL;
        M( num_cell, num_cell ) += coeff * dmass_dM;
        V[ num_cell ] -= coeff * mass;
    } );
}

DTP Vec<Vec<TF,2>> UTP::cell_boundaries( bool sorted_nums ) const {
    Vec<Vec<TF,2>> res( FromSize(), nb_cells() );

    if ( sorted_nums ) {
        for_each_cell( [&]( PI num_cell, Interval<TF> interval ) {
            res[ num_cell ][ 0 ] = interval.x0;
            res[ num_cell ][ 1 ] = interval.x1;
        } );
    } else {
        for_each_cell( [&]( PI num_cell, Interval<TF> interval ) {
            PI n = sorted_seed_nums[ num_cell ];
            res[ n ][ 0 ] = interval.x0;
            res[ n ][ 1 ] = interval.x1;
        } );
    }

    return res;
}

DTP VF UTP::barycenters_ap( const Density<TF> &density, bool sorted_nums, PI ni ) const {
    VF res( FromSize(), nb_cells() );

    auto diter = density->iterator();
    for_each_cell( [&]( PI num_cell, Interval<TF> interval ) {
        TF pint = 0, mass = 0;
        for_each_sub_interval( *diter, interval, num_cell, [&]( Interval<TF> interval ) {
            // diter.barycenter( pint, area, interval.x0, interval.x1 );
            for( PI i = 0; i < ni; ++i ) {
                const TF x = interval.x0 + ( interval.x1 - interval.x0 ) * ( i + 0.5 ) / ni;
                const TF a = ( interval.x1 - interval.x0 ) * diter->value( x );
                pint += a * x;
                mass += a;
            }            
        } );
        res[ sorted_nums ? num_cell : sorted_seed_nums[ num_cell ] ] = mass ? pint / mass : 0;
    } );

    return res;
}

DTP VF UTP::barycenters( const Density<TF> &density, bool sorted_nums ) const {
    VF res( FromSize(), nb_cells() );

    auto diter = density.iterator();
    for_each_cell( [&]( PI num_cell, Interval<TF> interval ) {
        TF pint = 0, area = 0;

        for_each_sub_interval( *diter, interval, num_cell, [&]( Interval<TF> interval ) {
            diter->barycenter( pint, area, interval.x0, interval.x1 );
        } );

        res[ sorted_nums ? num_cell : sorted_seed_nums[ num_cell ] ] = area ? pint / area : 0;
    } );

    return res;
}

DTP VF UTP::masses_ap( const Density<TF> &density, bool sorted_nums, PI ni ) const {
    VF res( FromSize(), nb_cells() );

    auto diter = density->iterator();
    for_each_cell( [&]( PI num_cell, Interval<TF> interval ) {
        TF area = 0;

        for_each_sub_interval( *diter, interval, num_cell, [&]( Interval<TF> interval ) {
            TF loc = 0;
            for( PI i = 0; i < ni; ++i ) {
                const TF x = interval.x0 + ( interval.x1 - interval.x0 ) * ( i + 0.5 ) / ni;
                loc += diter->value( x );
            }

            area += loc * ( interval.x1 - interval.x0 ) / ni;
        } );

        res[ sorted_nums ? num_cell : sorted_seed_nums[ num_cell ] ] = area;
    } );

    return res;
}

DTP VF UTP::masses( const Density<TF> &density, bool sorted_nums ) const {
    VF res( FromSize(), nb_cells() );

    auto diter = density.iterator();
    for_each_cell( [&]( PI num_cell, Interval<TF> interval ) {
        TF area = 0;

        for_each_sub_interval( *diter, interval, num_cell, [&]( Interval<TF> interval ) {
            area += diter->integral( interval.x0, interval.x1 );
        } );

        res[ sorted_nums ? num_cell : sorted_seed_nums[ num_cell ] ] = area;
    } );

    return res;
}

DTP void UTP::set_weights( const VF &weights, bool sorted_nums ) const {
    for( PI i = 0; i < nb_cells(); ++i )
        sorted_seed_weights[ i ] = weights[ sorted_nums ? i : sorted_seed_nums[ i ] ];
}

DTP VF UTP::get_weights( bool sorted_nums ) const {
    VF res( FromSize(), nb_cells() );
    for( PI i = 0; i < nb_cells(); ++i )
        res[ sorted_nums ? i : sorted_seed_nums[ i ] ] = sorted_seed_weights[ i ];
    return res;
}

#undef DTP
#undef UTP

} // namespace usdot
