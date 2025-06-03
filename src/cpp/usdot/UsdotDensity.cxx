#pragma once

#include "UsdotDensity.h"
#include <stdexcept>

namespace usdot {
    
#define DTP template<class TF>
#define UTP UsdotDensity<TF>

DTP UTP::UsdotDensity( TF beg_positions, TF end_positions, const VF &values ) : beg_positions( beg_positions ), end_positions( end_positions ), values( values ) {
    if ( values.empty() )
        throw std::runtime_error( "UsdotDensity: empty values" );
    if ( beg_positions >= end_positions )
        throw std::runtime_error( "UsdotDensity: beg_positions >= end_positions" );

    // normalization + primitives
    _compute_primitives();
}

DTP void UTP::_compute_primitives() {
    x_primitives.resize( values.size() + 1 );
    primitives.resize( values.size() + 1 );
    TF o = 0, x_o = 0;
    for( PI i = 0; i < values.size(); ++i ) {
        const TF x0 = beg_positions + ( end_positions - beg_positions ) * ( i + 0 ) / values.size();
        const TF x1 = beg_positions + ( end_positions - beg_positions ) * ( i + 1 ) / values.size();
        const TF v0 = values[ i + 0 ]; 

        x_primitives[ i ] = x_o;
        primitives[ i ] = o;

        x_o += v0 * ( x1 * x1 - x0 * x0 ) / 2;
        o += ( x1 - x0 ) * v0;
    }

    x_primitives.back() = x_o;
    primitives.back() = o;

    // normalization (actually, o should be equal to 1 by construction)
    for( TF &v : this->x_primitives )
        v /= o;
    for( TF &v : this->primitives )
        v /= o;
    for( TF &v : this->values )
        v /= o;

    //
    primitive_positions = primitives;
}

DTP TF UTP::inv_primitive( TF x ) const {
    return beg_positions + ( end_positions - beg_positions ) * primitive_positions.index_f( x ) / values.size();
}

DTP TF UTP::primitive( TF x ) const {
    auto index = ( x - beg_positions ) / ( end_positions - beg_positions ) * values.size();
    if ( index < 0 )
        return 0;
    if ( index >= values.size() )
        return 1;
    PI index_i = PI( index );
    TF index_f = index - index_i;
    return primitives[ index_i ] + values[ index_i ] * index_f * ( end_positions - beg_positions ) / values.size();
}

DTP TF UTP::x_primitive( TF x ) const {
    auto index = ( x - beg_positions ) / ( end_positions - beg_positions ) * values.size();
    if ( index < 0 )
        return 0;
    if ( index >= values.size() )
        return x_primitives.back();
    PI index_i = PI( index );
    const TF x0 = beg_positions + TF( end_positions - beg_positions ) * index_i / values.size();
    return x_primitives[ index_i ] + values[ index_i ] * ( x * x - x0 * x0 ) / 2;
}

DTP TF UTP::value( TF x ) const {
    auto index = ( x - beg_positions ) / ( end_positions - beg_positions ) * values.size();
    if ( index < 0 || index >= values.size() )
        return 0;
    return values[ PI( index ) ];
}

DTP TF UTP::x_integral( TF x0, TF x1 ) const {
    return x_primitive( x1 ) - x_primitive( x0 );
}

DTP TF UTP::integral( TF x0, TF x1 ) const {
    return primitive( x1 ) - primitive( x0 );
}

DTP TF UTP::min_x() const {
    return beg_positions;
}

DTP TF UTP::max_x() const {
    return end_positions;
}

DTP TF UTP::ptp_x() const {
    return max_x() - min_x();
}

DTP void UTP::plot( std::ostream &fs, std::string linestyle, double linewidth ) const {
    fs << "pyplot.plot( ["; 
    for( PI i = 0; i < values.size(); ++i ) {
        const TF x0 = beg_positions + ( end_positions - beg_positions ) * ( i + 0 ) / values.size();
        const TF x1 = beg_positions + ( end_positions - beg_positions ) * ( i + 1 ) / values.size();
        fs << x0 << ", " << x1 << ", ";
    }
    fs << " ], [";
    for( PI n = 0; n < values.size(); ++n ) 
        fs << values[ n ] << ", " << values[ n ] << ", ";
    fs << "], '" << linestyle << "', linewidth = " << linewidth << " )\n";

    // fs << "pyplot.plot( ["; 
    // for( PI i = 0; i < primitives.size(); ++i ) {
    //     const TF x0 = beg_positions + ( end_positions - beg_positions ) * ( i + 0 ) / values.size();
    //     fs << x0 << ", ";
    // }
    // fs << " ], [";
    // for( PI n = 0; n < primitives.size(); ++n ) 
    //     fs << primitives[ n ] << ", ";
    // fs << "] )\n";
}

#undef DTP
#undef UTP

} // namespace usdot
