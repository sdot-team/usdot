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
    for( PI i = 0; i + 1 < this->values.size(); ++i ) {
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
        + TF( i ) * ( f * v0 + f2 / 2 * ( v1 - v0 ) )
        + f2 / 2 * v0 + f3 / 3 * ( v1 - v0 );
}

DTP TF UTP::primitive( TF x ) const {
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

DTP TF UTP::width() const {
    return values.size() - 1;
}

DTP TF UTP::mass() const {
    return primitives.back();
}

DTP void UTP::get_inv_cdf( VF &inv_cdf_values, TF &mul_coeff, PI nb_bins ) const {
    using namespace std;

    mul_coeff = nb_bins / mass();
    inv_cdf_values.resize( nb_bins + 1 );
    for( PI n = 0; n + 1 < primitives.size(); ++n ) {
        const TF p0 = primitives[ n + 0 ];
        const TF p1 = primitives[ n + 1 ];
        const TF v0 = values[ n + 0 ];
        const TF v1 = values[ n + 1 ];
        
        const TF by = p0 * mul_coeff;
        const TF ey = p1 * mul_coeff;
        for( PI y = PI( ceil( by ) ); y < min( PI( floor( ey ) ), nb_bins ); ++y ) {
            if ( const TF a = v1 - v0 ) {
                const TF d = v0 * v0 - 2 * a * ( p0 - y / mul_coeff );
                inv_cdf_values[ y ] = n + ( a > 0 ? 
                    ( sqrt( d ) + v0 ) / a :
                    ( sqrt( d ) - v0 ) / a
                );
            } else if ( v0 ) {
                inv_cdf_values[ y ] = n + ( y / mul_coeff - p0 ) / v0;
            } else {
                inv_cdf_values[ y ] = n;
            }
        }
    }
    inv_cdf_values.back() = primitives.size() - 1;
}

DTP void UTP::plot( std::ostream &fs ) const {
    fs << "pyplot.plot( ["; 
    for( PI n = 0; n < values.size(); ++n )
        fs << n << ", ";    
    fs << " ], [";
    for( PI n = 0; n < values.size(); ++n )
        fs << values[ n ] << ", ";    
    fs << "] )\n";
}

#undef DTP
#undef UTP


} // namespace usdot
