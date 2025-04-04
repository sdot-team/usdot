#pragma once

#include "GridDensity.h"
#include <algorithm>

namespace usdot {

    
#define DTP template<class TF>
#define UTP GridDensity<TF>

DTP UTP::GridDensity( VF &&values ) : values( std::move( values ) ) {
    using namespace std;

    // primitives
    x_primitives.resize( this->values.size() );
    primitives.resize( this->values.size() );
    TF o = 0, x_o = 0, der_o = 0;
    for( TI i = 0; i + 1 < this->values.size(); ++i ) {
        const TF v0 = this->values[ i + 0 ]; 
        const TF v1 = this->values[ i + 1 ]; 
        const TF m = ( v0 + v1 ) / 2;

        x_primitives[ i ] = x_o;
        primitives[ i ] = o;

        // primitive( ( x0 + f ) * ( v0 + f * ( v1 - v0 ) ) )
        // m + primitive( f * ( v0 + f * ( v1 - v0 ) ) )
        // m + v0 / 2 + primitive( f * f * ( v1 - v0 ) )
        x_o += ( i * m + v0 / 2 + ( v1 - v0 ) / 3 );
        o += m;
    }

    x_primitives.back() = x_o;
    primitives.back() = o;
}

DTP TF UTP::x_primitive( TF x ) const {
    if ( x < 0 )
        return x_primitives.front();
    
    const TI i( x );
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
        + TF( i ) * ( f * v0 + f2 / 2 * ( v1 - v0 ) )
        + f2 / 2 * v0 + f3 / 3 * ( v1 - v0 );
}

DTP TF UTP::primitive( TF x ) const {
    if ( x < 0 )
        return 0;
    
    const TI i( x );
    if ( i >= primitives.size() - 1 )
        return primitives.back();

    // p0 + primitive( v0 + f * ( v1 - v0 ) )
    const TF f = x - i;
    const TF v0 = values[ i + 0 ];
    const TF v1 = values[ i + 1 ];
    return primitives[ i ] + f * v0 + f * f * ( v1 - v0 ) / 2;
}

DTP TF UTP::value( TF x ) const {
    if ( x < 0 || x >= values.size() - 1 )
        return 0;

    TI i( x );
    TF f = x - i;
    return values[ i ] * ( 1 - f ) + values[ i + 1 ] * f;
}

DTP TF UTP::x_integral( TF x0, TF x1 ) const {
    return x_primitive( x1 ) - x_primitive( x0 );
}

DTP TF UTP::integral( TF x0, TF x1 ) const {
    return primitive( x1 ) - primitive( x0 );
}

DTP TF UTP::mass() const {
    return primitives.back();
}

#undef DTP
#undef UTP


} // namespace usdot
