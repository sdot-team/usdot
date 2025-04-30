#pragma once

#include "utility/linspace.h"
#include "GridDensity.h"
#include <algorithm>
#include <stdexcept>
#include <limits>

namespace usdot {

    
#define DTP template<class TF,int regular_positions>
#define UTP GridDensity<TF,regular_positions>

DTP UTP::GridDensity( TF beg_positions, TF end_positions, const VF &values ) : GridDensity( linspace<TF>( beg_positions, end_positions, values.size() ), values ) {
}

DTP UTP::GridDensity( const VF &values ) : GridDensity( linspace<TF>( 0, values.size() - 1, values.size() ), values ) {
}

DTP UTP::GridDensity( const VF &positions, const VF &values ) : positions( positions ), values( values ) {
    using namespace std;

    // primitives
    x_primitives.resize( this->values.size() );
    primitives.resize( this->values.size() );
    TF o = 0, x_o = 0;
    for( PI i = 0; i + 1 < this->values.size(); ++i ) {
        const TF x0 = this->positions[ i + 0 ]; 
        const TF x1 = this->positions[ i + 1 ]; 
        const TF v0 = this->values[ i + 0 ]; 
        const TF v1 = this->values[ i + 1 ];
        const TF x6 = ( x0 + x1 ) / 6;

        x_primitives[ i ] = x_o;
        primitives[ i ] = o;

        // ( x1 - x0 ) * ( f * v0 * x0 + (f^3 * (v0 - v1) * (x0 - x1) ) / 3 + ( f^2 * (-2 * v0 * x0 + v1 * x0 + v0 * x1 ) ) / 2 )
        x_o += v0 * ( x1 * x6 - x0 * x0 / 3 ) + v1 * ( x1 * x1 / 3 - x0 * x6 );
        o += ( x1 - x0 ) * ( v0 + v1 ) / 2;
    }

    x_primitives.back() = x_o;
    primitives.back() = o;

    // normalization
    for( TF &v : this->x_primitives )
        v /= o;
    for( TF &v : this->primitives )
        v /= o;
    for( TF &v : this->values )
        v /= o;

    // opt_pos
    const PI opt_pos_size = this->values.size() - 1;
    opt_pos_begs.resize( opt_pos_size, numeric_limits<PI>::max() );
    opt_pos_beg = this->positions.front();
    opt_pos_end = this->positions.back();
    opt_pos_coeff = opt_pos_size / ( opt_pos_end - opt_pos_beg );
    for( PI i = 0; i + 1 < this->values.size(); ++i ) {
        const PI x0 = floor( ( this->positions[ i + 0 ] - opt_pos_beg ) * opt_pos_coeff );
        const PI x1 = ceil( ( this->positions[ i + 1 ] - opt_pos_beg ) * opt_pos_coeff );
        for( PI n = x0; n < min( x1, opt_pos_size ); ++n )
            opt_pos_begs[ n ] = min( opt_pos_begs[ n ], i );
    }

    // opt_cdf
    const PI mul_opt_cdf = 2;
    const PI opt_cdf_size = mul_opt_cdf * this->values.size();
    opt_cdf_begs.resize( opt_cdf_size, numeric_limits<PI>::max() );
    for( PI i = 0; i + 1 < this->primitives.size(); ++i ) {
        const PI x0 = floor( this->primitives[ i + 0 ] * opt_cdf_size );
        const PI x1 = ceil( this->primitives[ i + 1 ] * opt_cdf_size );
        for( PI n = x0; n < min( x1, opt_cdf_size ); ++n )
            opt_cdf_begs[ n ] = min( opt_cdf_begs[ n ], i );
    }
}

DTP TF UTP::x_primitive( TF x ) const {
    if ( regular_positions ) {
        if ( regular_positions == 1 )
            x = ( x - opt_pos_beg ) * opt_pos_coeff;
        if ( x < 0 )
            return 0;

        const PI i( x );
        if ( i >= x_primitives.size() - 1 )
            return x_primitives.back();

        TF p0 = i + 0;
        TF p1 = i + 1;
        if ( regular_positions == 1 ) {
            p0 = opt_pos_beg + p0 / opt_pos_coeff;
            p1 = opt_pos_beg + p1 / opt_pos_coeff;
        }
        
        const TF v0 = values[ i + 0 ];
        const TF v1 = values[ i + 1 ];
        const TF dv = v1 - v0;
        const TF dp = p1 - p0;
        const TF f = x - i;
        const TF f2 = f * f;
        const TF f3 = f2 * f;
        TF res = ( f * v0 * p0 + f3 * dv * dp / 3 + ( f2 * ( v1 * p0 + v0 * p1 - 2 * v0 * p0 ) ) / 2 );
        if ( regular_positions == 1 )
            res /= opt_pos_coeff;
        return x_primitives[ i ] + res;
    }

    if ( x < opt_pos_beg )
        return 0;
    if ( x > opt_pos_end )
        return x_primitives.back();

    const PI ox = std::min( PI( ( x - opt_pos_beg ) * opt_pos_coeff ), opt_pos_begs.size() - 1 );
    for( PI i = opt_pos_begs[ ox ]; ; ++i ) {
        if ( positions[ i + 1 ] >= x ) {
            const TF p0 = positions[ i + 0 ];
            const TF p1 = positions[ i + 1 ];
            if ( const TF dp = p1 - p0 ) {
                const TF v0 = values[ i + 0 ];
                const TF v1 = values[ i + 1 ];
                const TF dv = v1 - v0;
                const TF f = ( x - p0 ) / dp;
                const TF f2 = f * f;
                const TF f3 = f2 * f;
        
                return x_primitives[ i ] + dp * ( f * v0 * p0 + f3 * dv * dp / 3 + ( f2 * ( v1 * p0 + v0 * p1 - 2 * v0 * p0 ) ) / 2 );
            }
            return x_primitives[ i ];
        }
    }
}

DTP TF UTP::primitive( TF x ) const {
    if ( regular_positions ) {
        if ( regular_positions == 1 )
            x = ( x - opt_pos_beg ) * opt_pos_coeff;
        if ( x < 0 )
            return 0;

        const PI i( x );
        if ( i >= primitives.size() - 1 )
            return primitives.back();
        
        TF f = x - i;
        const TF v0 = values[ i + 0 ];
        const TF v1 = values[ i + 1 ];
        TF res = f * v0 + f * f * ( v1 - v0 ) / 2;
        if ( regular_positions == 1 )
            res /= opt_pos_coeff;
        return primitives[ i ] + res;
    }

    if ( x < opt_pos_beg )
        return 0;
    if ( x > opt_pos_end )
        return primitives.back();

    const PI ox = std::min( PI( ( x - opt_pos_beg ) * opt_pos_coeff ), opt_pos_begs.size() - 1 );
    for( PI i = opt_pos_begs[ ox ]; ; ++i ) {
        if ( positions[ i + 1 ] >= x ) {
            const TF x0 = positions[ i + 0 ];
            const TF x1 = positions[ i + 1 ];
            if ( const TF dx = x1 - x0 ) {
                const TF v0 = values[ i + 0 ];
                const TF v1 = values[ i + 1 ];
                const TF xd = x - x0;
              
                TF res = xd * v0 + xd * xd / dx * ( v1 - v0 ) / 2;
                return primitives[ i ] + res;
            }
            return primitives[ i ];
        }
    }
}

DTP TF UTP::value( TF x ) const {
    if ( regular_positions ) {
        if ( regular_positions == 1 )
            x = ( x - opt_pos_beg ) * opt_pos_coeff;
        if ( x < 0 || x >= values.size() - 1 )
            return 0;
    
        PI i( x );
        TF f = x - i;
        return values[ i ] * ( 1 - f ) + values[ i + 1 ] * f;
    }

    if ( x < opt_pos_beg || x > opt_pos_end )
        return 0;

    const PI ox = std::min( PI( ( x - opt_pos_beg ) * opt_pos_coeff ), opt_pos_begs.size() - 1 );
    for( PI i = opt_pos_begs[ ox ]; ; ++i ) {
        if ( positions[ i + 1 ] >= x ) {
            if ( const TF de = positions[ i + 1 ] - positions[ i + 0 ] ) {
                TF f = ( x - positions[ i + 0 ] ) / de;
                return values[ i + 0 ] * ( 1 - f ) + values[ i + 1 ] * f;
            }
            return values[ i ];
        }
    }
}

DTP TF UTP::position( PI i ) const {
    if ( regular_positions == 2 )
        return i;
    if ( regular_positions == 1 )
        return opt_pos_beg + i / opt_pos_coeff;
    return positions[ i ];
}

DTP TF UTP::inv_cdf_no_check( TF u ) const {
    using namespace std;

    const PI os = opt_cdf_begs.size() - 1;
    const PI ox = std::min( PI( u * os ), os );
    for( PI i = opt_cdf_begs[ ox ]; ; ++i ) {
        if ( primitives[ i + 1 ] >= u ) {
            const TF p0 = primitives[ i + 0 ];
            const TF p1 = primitives[ i + 1 ];
            if ( const TF dp = p1 - p0 ) {
                const TF x0 = position( i + 0 );
                const TF x1 = position( i + 1 );
                const TF v0 = values[ i + 0 ];
                const TF v1 = values[ i + 1 ];
                const TF dx = x1 - x0;
                const TF dv = v1 - v0;
                
                if ( x0 == x1 || v0 == v1 ) {
                    const TF f = ( u - p0 ) / dp;
                    return x0 * ( 1 - f ) + x1 * f;
                }

                const TF c = ( primitives[ i ] - u );
                const TF a = dv / dx;

                const TF d = max( TF( 0 ), v0 * v0 - 2 * a * c );
                const TF x = ( sqrt( d ) - v0 ) / a; // 
                return x0 + x;

            }
            return positions[ i ];
        }
    }
}

DTP TF UTP::der_inv_cdf_no_check( TF u ) const {
    using namespace std;

    const PI os = opt_cdf_begs.size() - 1;
    const PI ox = std::min( PI( u * os ), os );
    for( PI i = opt_cdf_begs[ ox ]; ; ++i ) {
        if ( primitives[ i + 1 ] >= u ) {
            const TF p0 = primitives[ i + 0 ];
            const TF p1 = primitives[ i + 1 ];
            if ( const TF dp = p1 - p0 ) {
                const TF x0 = position( i + 0 );
                const TF x1 = position( i + 1 );
                const TF v0 = values[ i + 0 ];
                const TF v1 = values[ i + 1 ];
                const TF dx = x1 - x0;
                const TF dv = v1 - v0;
                
                if ( x0 == x1 || v0 == v1 )
                    return ( x1 - x0 ) / dp;

                const TF c = ( primitives[ i ] - u );
                const TF a = dv / dx;

                const TF d = max( TF( 0 ), v0 * v0 - 2 * a * c );
                return d ? pow( d, - 0.5 ) : sqrt( numeric_limits<TF>::max() );
            }
            return sqrt( numeric_limits<TF>::max() );
        }
    }
}

DTP TF UTP::inv_cdf( TF u ) const {
    if ( u > 1 )
        return 1;
    if ( u < 0 )
        return 0;
    return inv_cdf_no_check( u );
}

DTP TF UTP::inv_cdf( TF u, TF ext_density ) const {
    if ( u > 1 )
        return opt_pos_end + ( u - 1 ) / ext_density;
    if ( u < 0 )
        return opt_pos_beg + u / ext_density;
    return inv_cdf_no_check( u );
}

DTP TF UTP::der_inv_cdf( TF u, TF ext_density ) const {
    if ( u > 1 )
        return 1 / ext_density;
    if ( u < 0 )
        return 1 / ext_density;
    return der_inv_cdf_no_check( u );
}

DTP TF UTP::cdf( TF x, TF ext_density ) const {
    if ( x < opt_pos_beg )
        return ( x - opt_pos_beg ) * ext_density;
    if ( x > opt_pos_end )
        return 1 + ( x - opt_pos_end ) * ext_density;
    return cdf( x );
}

DTP TF UTP::x_integral( TF x0, TF x1 ) const {
    return x_primitive( x1 ) - x_primitive( x0 );
}

DTP TF UTP::integral( TF x0, TF x1 ) const {
    return cdf( x1 ) - cdf( x0 );
}

DTP TF UTP::min_x() const {
    return opt_pos_beg;
}

DTP TF UTP::max_x() const {
    return opt_pos_end;
}

DTP TF UTP::ptp_x() const {
    return max_x() - min_x();
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
