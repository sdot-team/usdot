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
    if ( const TF ey = cy.back() ) {
        for( auto &y : cy ) 
            y /= ey;
        for( auto &z : cz ) 
            z /= ey;
    }

    // Vec<TF> c;
    // for( TF x : Vec<TF>::linspace( 0, 1, 21 ) )
    //     c << cdf_value( x );
    // P( c );

    // Vec<TF> d;
    // for( TF x : Vec<TF>::linspace( 0, 1, 21 ) )
    //     d << ( cdf_value( x + 1e-4 ) - cdf_value( x ) ) / 1e-4;
    // P( d );
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
        const TF x0 = cx[ i + 0 ];
        const TF x1 = cx[ i + 1 ];
        if ( TF dx = x1 - x0 ) {
            const TF fr = ( x - x0 ) / dx;
            const TF y0 = cy[ i + 0 ];
            const TF y1 = cy[ i + 1 ];
            const TF zm = cz[ i ];

            y = y0 + fr * ( y1 - y0 ) + 4 * fr * ( 1 - fr ) * zm;        
        } else
            y = ( cy[ i + 0 ] + cy[ i + 1 ] ) / 2;

        // boundaries
        const TF mr = mass_ratios[ n ];
        const TF m2 = mr / 2;
        y = min( y, 1 - m2 );
        y = max( y, m2 );

        // append and check if touching
        append_cells( n, n + 1, y - m2, mr );
    }

    // helper to get x (seed coord) from y (cdf axis)
    PI i = cy.size() - 2;
    auto get_x = [&]( TF y ) {
        while ( cy[ i + 0 ] > y )
            --i;

        // y0 + x * ( y1 - y0 ) + 4 * x * ( 1 - x ) * zm = y
        // y - y0 - x * ( y1 - y0 + 4 * zm ) + 4 * x * x * zm
        const TF x0 = cx[ i + 0 ];
        const TF x1 = cx[ i + 1 ];
        const TF y0 = cy[ i + 0 ];
        const TF y1 = cy[ i + 1 ];
        
        if ( const TF zm = cz[ i ] ) {
            const TF b = y1 - y0 + 4 * zm;
            const TF d = max( 0, b * b - 16 * zm * ( y - y0 ) );
            const TF fr = ( b - sqrt( d ) ) / ( 8 * zm );
            return x0 + ( x1 - x0 ) * fr;
        }
        
        const TF fr = ( y - y0 ) / ( y1 - y0 );
        return x0 + ( x1 - x0 ) * fr;
    };

    // get weights
    for( Agglomeration *item = last; item; item = item->prev ) {
        if ( item->len_n() == 0 )
            continue;

        P( item->beg_y, item->end_y(), item->len_n() );

        // data from the last cell
        const PI n = item->end_n - 1;
        TF c = seed_coords[ n ];

        TF end_y = item->end_y(), beg_y = end_y - mass_ratios[ n ];
        TF end_x = get_x( end_y ), beg_x = get_x( beg_y );

        P( beg_x, end_x, c );

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
            
            P( beg_x, end_x, c, d );
            
            const TF min_w = max( pow( c - beg_x, 2 ), pow( end_x - c, 2 ) );
            if ( min_w > w )
                w_to_add = max( w_to_add, min_w - w );

            // store the result
            seed_weights[ n - 1 ] = w;
        }
        
        P( w_to_add );


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

DTP TF UTP::cdf_value( TF x ) const {
    if ( cx.empty() )
        return 0;

    if ( x <= cx.front() )
        return 0;
    if ( x >= cx.back() )
        return 1;

    PI i = 0;
    while ( cx[ i + 1 ] < x )
        ++i;

    const TF x0 = cx[ i + 0 ];
    const TF x1 = cx[ i + 1 ];
    if ( x0 == x1 )
        return cy[ i + 0 ];
    
    const TF fr = ( x - x0 ) / ( x1 - x0 );
    const TF y0 = cy[ i + 0 ];
    const TF y1 = cy[ i + 1 ];
    const TF zm = cz[ i ];

    return y0 + fr * ( y1 - y0 ) + 4 * fr * ( 1 - fr ) * zm;
}

#undef DTP
#undef UTP

} // namespace usdot
