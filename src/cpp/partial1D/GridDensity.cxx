#pragma once

#include <algorithm>
#include <tl/support/P.h>
#include "GridDensity.h"

namespace usdot {

    
#define DTP template<class TF>
#define UTP GridDensity<TF>

DTP UTP::GridDensity( const Vec &original_values, TF filter, PI mul_x, TF cut_ratio ) {
    using namespace std;

    // 
    Vec der_forward; der_forward.reserve( 2 * mul_x * values.size() );
    Vec forward; forward.reserve( 2 * mul_x * values.size() );
    TF der_prev = 0, prev = 0, m = 0;
    for( TF v : original_values ) {
        for( PI x = 0; x < mul_x; ++x ) {
            const TF n = filter * prev + ( 1 - filter ) * v;
            const TF der_n = filter * der_prev + prev - v;
            der_forward.push_back( der_n );
            forward.push_back( n );

            m = max( m, n );

            der_prev = der_n;
            prev = n;
        }
    }
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
    offset = values.size() - forward.size();
}

DTP TF UTP::value( TF x ) const {
    x += offset;
    if ( x < 0 || x >= values.size() - 1 )
        return 0;

    PI i( x );
    TF f = x - i;
    return values[ i ] * ( 1 - f ) + values[ i + 1 ] * f;
}

DTP TF UTP::der_value( TF x ) const {
    x += offset;
    if ( x < 0 || x >= der_values.size() - 1 )
        return 0;

    PI i( x );
    TF f = x - i;
    return der_values[ i ] * ( 1 - f ) + der_values[ i + 1 ] * f;
}


#undef DTP
#undef UTP


} // namespace usdot
