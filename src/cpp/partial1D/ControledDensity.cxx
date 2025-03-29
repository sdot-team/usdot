#pragma once

#include "ControledDensity.h"
#include <tl/support/P.h>
#include <algorithm>

namespace usdot {

#define DTP template<class TF>
#define UTP ControledDensity<TF>

DTP UTP::ControledDensity( const Vec &source, TF ratio_at_end, TF target_mass ) {
    using namespace std;

    // (pixel_ratio)^nb_pixels = ratio_at_end
    TF pixel_ratio = ratio_at_end ? pow( ratio_at_end, 1. / source.size() ) : 0;

    // forward
    values.reserve( 2 * source.size() );
    TF prev = 0;
    for( TF v : source ) {
        const TF n = max( prev, v );
        prev = n * pixel_ratio;
        values.push_back( n );
    }
    while ( prev > 1e-9 ) {
        values.push_back( prev );
        prev *= pixel_ratio;
    }

    // backward
    std::reverse( values.begin(), values.end() );
    const PI old_size = values.size();
    prev = 0;
    for( TF &v : values ) {
        if ( v < prev )
            v = prev;
        prev = v * pixel_ratio;
    }
    while ( prev > 1e-9 ) {
        values.push_back( prev );
        prev *= pixel_ratio;
    }

    std::reverse( values.begin(), values.end() );
    offset = values.size() - old_size;

    // primitives
    x_primitives.resize( values.size() );
    primitives.resize( values.size() );
    TF o = 0, x_o = 0, der_o = 0;
    for( PI i = 0; i + 1 < values.size(); ++i ) {
        x_primitives[ i ] = x_o;
        primitives[ i ] = o;

        const TF v0 = values[ i + 0 ]; 
        const TF v1 = values[ i + 1 ]; 
        const TF m = ( v0 + v1 ) / 2;
        // primitive( ( x0 + f ) * ( v0 + f * ( v1 - v0 ) ) )
        // m + primitive( f * ( v0 + f * ( v1 - v0 ) ) )
        // m + v0 / 2 + primitive( f * f * ( v1 - v0 ) )
        x_o += ( ( TF( i ) - offset ) * m + v0 / 2 + ( v1 - v0 ) / 3 );
        o += m;
    }
    x_primitives.back() = x_o;
    primitives.back() = o;
    
    // normalize
    if ( target_mass > 0 ) {
        const TF mul = target_mass / integral( 0, source.size() - 1 );
        for( auto &v : x_primitives ) v *= mul;
        for( auto &v : primitives ) v *= mul;
        for( auto &v : values ) v *= mul;
    }
}

DTP TF UTP::x_primitive( TF x ) const {
    x += offset;
    if ( x < 0 )
        return x_primitives.front();
    
    const PI i( x );
    if ( i >= primitives.size() - 1 )
        return x_primitives.back();

    // px0 + primitive( ( i - offset + f ) * ( v0 + f * ( v1 - v0 ) ) )
    // px0 + ( i - offset ) * ( f * v0 + f^2/2 * ( v1 - v0 ) ) + primitive( f * ( v0 + f * ( v1 - v0 ) ) )
    const TF f = x - i;
    const TF f2 = f * f;
    const TF f3 = f2 * f;
    const TF v0 = values[ i + 0 ];
    const TF v1 = values[ i + 1 ];
    return x_primitives[ i ] 
        + ( TF( i ) - offset ) * ( f * v0 + f2 / 2 * ( v1 - v0 ) )
        + f2 / 2 * v0 + f3 / 3 * ( v1 - v0 );
}

DTP TF UTP::primitive( TF x ) const {
    x += offset;
    if ( x < 0 )
        return 0;
    
    const PI i( x );
    if ( i >= primitives.size() - 1 )
        return primitives.back();

    // p0 + primitive( v0 + f * ( v1 - v0 ) )
    const TF f = x - i;
    const TF v0 = values[ i + 0 ];
    const TF v1 = values[ i + 1 ];
    return primitives[ i ] + f * v0 + f * f * ( v1 - v0 ) / 2;
}

DTP TF UTP::value( TF x ) const {
    x += offset;
    if ( x < 0 || x >= values.size() - 1 )
        return 0;

    PI i( x );
    TF f = x - i;
    return values[ i ] * ( 1 - f ) + values[ i + 1 ] * f;
}

DTP TF UTP::x_integral( TF x0, TF x1 ) const {
    return x_primitive( x1 ) - x_primitive( x0 );
}

DTP TF UTP::integral( TF x0, TF x1 ) const {
    return primitive( x1 ) - primitive( x0 );
}

DTP TF UTP::contrast( PI beg, PI end ) const {
    //TF mi
    return 0;
}

#undef DTP
#undef UTP

} // namespace usdot
