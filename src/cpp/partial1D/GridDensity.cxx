#pragma once

#include <algorithm>
#include <tl/support/P.h>
#include "GridDensity.h"

namespace usdot {

    
#define DTP template<class TF>
#define UTP GridDensity<TF>

DTP UTP::GridDensity( const Vec &original_values, TF filter, PI mul_x, TF cut_ratio, TF target_mass ) {
    using namespace std;

    //
    auto interp = [&]( auto &&func ) {
        for( PI x = 0; x < mul_x; ++x ) {
            TF a = TF( x ) / mul_x;
            func( a * original_values.front() );
        }
        for( PI i = 0; i + 1 < original_values.size(); ++i ) {
            for( PI x = 0; x < mul_x; ++x ) {
                TF a = TF( x ) / mul_x;
                func( ( 1 - a ) * original_values[ i + 0 ] + a * original_values[ i + 1 ] );
            }
        }
        for( PI x = 0; x <= mul_x; ++x ) {
            TF a = TF( x ) / mul_x;
            func( ( 1 - a ) * original_values.back() );
        }
    };

    // 
    Vec der_forward; der_forward.reserve( 2 * mul_x * values.size() );
    Vec forward; forward.reserve( 2 * mul_x * values.size() );
    TF der_prev = 0, prev = 0, m = 0;
    interp( [&]( TF v ) {
        const TF n = filter * prev + ( 1 - filter ) * v;
        const TF der_n = filter * der_prev + prev - v;
        der_forward.push_back( der_n );
        forward.push_back( n );

        m = max( m, n );

        der_prev = der_n;
        prev = n;
    } );
    while ( prev > m * cut_ratio ) {
        const TF n = filter * prev;
        const TF der_n = prev + filter * der_prev;

        der_forward.push_back( der_n );
        forward.push_back( n );
        
        der_prev = der_n;
        prev = n;
    }

    //
    der_values.reserve( 2 * forward.size() - mul_x * values.size() );
    values.reserve( 2 * forward.size() - mul_x * values.size() );
    for( PI i = forward.size(); i--; ) {
        const TF n = filter * prev + ( 1 - filter ) * forward[ i ];
        const TF der_n = filter * der_prev + prev - forward[ i ] + ( 1 - filter ) * der_forward[ i ];

        der_values.push_back( der_n );
        values.push_back( n );

        der_prev = der_n;
        prev = n;
    }
    while ( prev > m * cut_ratio ) {
        const TF n = filter * prev;
        const TF der_n = prev + filter * der_prev;

        der_values.push_back( der_n );
        values.push_back( n );

        der_prev = der_n;
        prev = n;
    }

    //
    std::reverse( der_values.begin(), der_values.end() );
    std::reverse( values.begin(), values.end() );
    offset = values.size() - forward.size() + mul_x;

    // primitives
    der_primitives.resize( values.size() );
    x_primitives.resize( values.size() );
    primitives.resize( values.size() );
    TF o = 0, x_o = 0, der_o = 0;
    for( PI i = 0; i + 1 < values.size(); ++i ) {
        der_primitives[ i ] = der_o;
        x_primitives[ i ] = x_o;
        primitives[ i ] = o;

        const TF dv0 = der_values[ i + 0 ]; 
        const TF dv1 = der_values[ i + 1 ]; 
        const TF dm = ( dv0 + dv1 ) / 2;
        const TF v0 = values[ i + 0 ]; 
        const TF v1 = values[ i + 1 ]; 
        const TF m = ( v0 + v1 ) / 2;
        // primitive( ( x0 + f ) * ( v0 + f * ( v1 - v0 ) ) )
        // m + primitive( f * ( v0 + f * ( v1 - v0 ) ) )
        // m + v0 / 2 + primitive( f * f * ( v1 - v0 ) )
        x_o += ( ( TF( i ) - offset ) * m + v0 / 2 + ( v1 - v0 ) / 3 ) / mul_x;
        der_o += dm;
        o += m;
    }
    der_primitives.back() = der_o;
    x_primitives.back() = x_o;
    primitives.back() = o;

    // normalize
    const TF mul = target_mass / integral( 0, original_values.size() - 1 );
    for( auto &v : der_primitives ) v *= mul;
    for( auto &v : der_values ) v *= mul;

    for( auto &v : x_primitives ) v *= mul;
    for( auto &v : primitives ) v *= mul;
    for( auto &v : values ) v *= mul;
}

DTP TF UTP::der_primitive( TF x ) const {
    x += offset;
    if ( x < 0 )
        return 0;
    
    const PI i( x );
    if ( i >= primitives.size() - 1 )
        return der_primitives.back();

    // p0 + primitive( v0 + f * ( v1 - v0 ) )
    const TF f = x - i;
    const TF v0 = der_values[ i + 0 ];
    const TF v1 = der_values[ i + 1 ];
    return der_primitives[ i ] + f * v0 + f * f * ( v1 - v0 ) / 2;
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

DTP TF UTP::der_value( TF x ) const {
    x += offset;
    if ( x < 0 || x >= der_values.size() - 1 )
        return 0;

    PI i( x );
    TF f = x - i;
    return der_values[ i ] * ( 1 - f ) + der_values[ i + 1 ] * f;
}

DTP TF UTP::value( TF x ) const {
    x += offset;
    if ( x < 0 || x >= values.size() - 1 )
        return 0;

    PI i( x );
    TF f = x - i;
    return values[ i ] * ( 1 - f ) + values[ i + 1 ] * f;
}

DTP TF UTP::der_integral( TF x0, TF x1 ) const {
    return der_primitive( x1 ) - der_primitive( x0 );
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
