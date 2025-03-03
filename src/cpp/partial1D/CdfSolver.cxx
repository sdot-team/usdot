#pragma once

#include "CdfSolver.h"
#include <limits>

namespace usdot {
    
#define DTP template<class TF>
#define UTP CdfSolver<TF>

DTP UTP::CdfSolver( CdfApproximation<TF> &&cdf, const Vec<TF> &seed_coords, const Vec<TF> &mass_ratios ) : seed_coords( seed_coords ), mass_ratios( mass_ratios ), cx( std::move( cdf.xs ) ), cy( std::move( cdf.ys ) ), cz( std::move( cdf.zs ) ), last( nullptr ) {
    if ( cy.empty() )
        return;

    // normalization
    if ( const TF ey = cy.back() )
        for( auto &y : cy ) 
            y /= ey;    
}

DTP void UTP::solve( Vec<TF> &seed_weights ) {
    using namespace std;

    seed_weights.resize( nb_cells() );
    last = nullptr;

    // find the `y` position (position in the cumulative density function) of the center of the cells, assuming first they are not touching
    for( PI n = 0, i = 0; n < nb_cells(); ++n ) {
        TF x = max( cx.front(), min( cx.back(), seed_coords[ n ] ) );
        while( cx[ i + 1 ] < x )
            ++i;
        
        // single cell position
        TF y;
        if ( TF dx = cx[ i + 1 ] - cx[ i + 0 ] )
            y = cy[ i + 0 ] + ( cy[ i + 1 ] - cy[ i + 0 ] ) * ( x - cx[ i + 0 ] ) / dx;
        else
            y = ( cy[ i + 0 ] + cy[ i + 1 ] ) / 2;

        // boundaries
        const TF mr = mass_ratios[ n ];
        const TF m2 = mr / 2;
        y = min( y, 1 - m2 );
        y = max( y, m2 );

        // append and check if touching
        append_cells( n, n + 1, y - m2, mr);
    }

    // helper to get x (seed coord) from y (cdf axis)
    PI i = cy.size() - 2;
    auto get_x = [&]( TF y ) {
        while ( cy[ i + 0 ] > y )
            --i;
        return cx[ i + 0 ] + ( cx[ i + 1 ] - cx[ i + 0 ] ) * ( y - cy[ i + 0 ] ) / ( cy[ i + 1 ] - cy[ i + 0 ] );
    };

    // get weights
    for( Agglomeration *item = last; item; item = item->prev ) {
        if ( item->len_n() == 0 )
            continue;

        // data from the last cell
        const PI n = item->end_n - 1;
        TF c = seed_coords[ n ];

        TF end_y = item->end_y(), beg_y = end_y - mass_ratios[ n ];
        TF end_x = get_x( end_y ), beg_x = get_x( beg_y );

        // compute the weight from the required radius
        TF w = max( pow( end_x - c, 2 ), pow( c - beg_x, 2 ) );
        seed_weights[ item->end_n - 1 ] = w;

        //  
        TF w_to_add = 0;
        for( PI n = item->end_n; --n > item->beg_n; ) {
            // update bounds
            end_y = std::exchange( beg_y, max( 0, beg_y - mass_ratios[ n ] ) );
            end_x = std::exchange( beg_x, get_x( beg_y ) );
            
            // update w
            TF d = std::exchange( c, seed_coords[ n - 1 ] );
            w -= ( d + c - 2 * end_x ) * ( d - c );
            
            const TF min_w = max( pow( c - beg_x, 2 ), pow( end_x - c, 2 ) );
            if ( min_w > w )
                w_to_add = max( w_to_add, min_w - w );

            // store the result
            seed_weights[ n - 1 ] = w;
        }

        if ( w_to_add )
            for( PI n = item->beg_n; n < item->end_n; ++n )
                seed_weights[ n ] += w_to_add;
    }
}

DTP void UTP::append_cells( PI beg_n, PI end_n, TF beg_y, TF len_y ) {
    // non touching cell ?
    if ( last == nullptr || last->end_y() < beg_y ) {
        Agglomeration *item = pool.create<Agglomeration>();
        item->beg_y = beg_y;
        item->len_y = len_y;
        item->beg_n = beg_n;
        item->end_n = end_n;
        item->prev = last;
        last = item;
        return;
    }

    // append to the new cell(s)
    TF move_back = ( last->end_y() - beg_y ) * 1 / ( last->len_n() + 1 );
    last->len_y += len_y;
    last->end_n = end_n;

    // while we have an Agglomerate to move back
    while ( true ) {
        last->beg_y -= move_back;
        if ( last->beg_y < 0 )
            last->beg_y = 0;
        if ( last->end_y() > 1 )
            last->beg_y = max( 0, 1 - last->len_y );

        // break 
        Agglomeration *prev = last->prev;
        if ( prev == nullptr || prev->end_y() < last->beg_y )
            break;

        // append the last set
        move_back = ( prev->end_y() - last->beg_y ) * last->len_n() / ( prev->len_n() + last->len_n() );
        prev->len_y += last->len_y;
        prev->end_n = last->end_n;
        last->prev = prev->prev;
        last = prev;
    }
}

#undef DTP
#undef UTP

} // namespace usdot
