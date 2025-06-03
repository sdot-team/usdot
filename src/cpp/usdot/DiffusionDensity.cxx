#pragma once

#include "utility/linspace.h"
#include "DiffusionDensity.h"
#include <algorithm>
#include <stdexcept>
// #include <limits>

namespace usdot {
    
#define DTP template<class TF>
#define UTP DiffusionDensity<TF>

DTP UTP::DiffusionDensity( TF beg_original_positions, TF end_original_positions, const VF &original_values, TF start_flattening_ratio ) : DiffusionDensity( linspace( beg_original_positions, end_original_positions, original_values.size() ), original_values, start_flattening_ratio ) {
}

DTP UTP::DiffusionDensity( const VF &original_values, TF start_flattening_ratio ) : DiffusionDensity( 0, 1, original_values, start_flattening_ratio ) {
}

DTP UTP::DiffusionDensity( const VF &positions, const VF &original_values, TF start_flattening_ratio ) : indexed_positions( positions ), original_values( original_values ) {
    this->current_flattening_ratio = -1;
    this->coeff_flattening_ratio = 3;

    // extended positions
    const PI n = positions.size();
    const PI e = n - 1;
    this->extended_positions.resize( n + 2 * e );

    const TF beg_x = positions.front();
    const TF end_x = positions.back();
    for( PI i = 0; i < e; ++i ) {
        TF f = TF( i ) / e;
        this->extended_positions[ i ] = ( 2 - f ) * beg_x + ( f - 1 ) * end_x;
    }
    for( PI i = 0; i < n; ++i )
        this->extended_positions[ e + i ] = positions[ i ];
    for( PI i = 0; i < e; ++i ) {
        TF f = TF( i + 1 ) / e;
        this->extended_positions[ e + n + i ] = ( 1 - f ) * beg_x + ( 1 + f ) * end_x;
    }

    // normalization
    TF o = 0;
    for( PI i = 0; i + 1 < this->original_values.size(); ++i )
        o += ( positions[ i + 1 ] - positions[ i + 0 ] ) * ( original_values[ i + 0 ] + original_values[ i + 1 ] ) / 2;
    for( auto &v : this->original_values )
        v /= o;

    // start with no flattening
    set_flattening_ratio( start_flattening_ratio );
}

DTP void UTP::set_flattening_ratio( TF flattening_ratio ) {
    if ( current_flattening_ratio == flattening_ratio )
        return;

    current_flattening_ratio = flattening_ratio;
    _compute_values_for( flattening_ratio );
    _compute_primitives();
}

DTP void UTP::_compute_primitives() {
    x_primitives.resize( values.size() );
    primitives.resize( values.size() );
    TF o = 0, x_o = 0;
    for( PI i = 0; i + 1 < values.size(); ++i ) {
        const TF x0 = indexed_positions.positions[ i + 0 ]; 
        const TF x1 = indexed_positions.positions[ i + 1 ]; 
        const TF v0 = values[ i + 0 ]; 
        const TF v1 = values[ i + 1 ];
        const TF x6 = ( x0 + x1 ) / 6;

        x_primitives[ i ] = x_o;
        primitives[ i ] = o;

        // ( x1 - x0 ) * ( f * v0 * x0 + (f^3 * dv * (x0 - x1) ) / 3 + ( f^2 * (-2 * v0 * x0 + v1 * x0 + v0 * x1 ) ) / 2 )
        x_o += v0 * ( x1 * x6 - x0 * x0 / 3 ) + v1 * ( x1 * x1 / 3 - x0 * x6 );
        o += ( x1 - x0 ) * ( v0 + v1 ) / 2;
    }

    x_primitives.back() = x_o;
    primitives.back() = o;

    // normalization (actually, o should be equal to 1 by construction)
    if ( normalization ) {
        sys_div = o;
        for( TF &v : this->x_primitives )
            v /= sys_div;
        for( TF &v : this->primitives )
            v /= sys_div;
        for( TF &v : this->values )
            v /= sys_div;
    } else
        sys_div = 1;

    max_value = 0;
    for( const TF &v : values )
        if ( v > max_value )
            max_value = v;
}


DTP void UTP::compute_derivatives( PI nb_derivatives ) {
    const auto &positions = indexed_positions.positions;
    const PI s = extended_positions.size();
    const TF d = pow( ptp_x(), -2 );
    const PI n = positions.size();
    const PI e = ( s - n ) / 2;

    der_primitives.clear();
    der_values.clear();

    // 0
    if ( current_flattening_ratio == 0 ) {
        throw std::runtime_error( "TODO der 0" );
    }

    VF der_1;
    if ( nb_derivatives >= 1 ) {
        VF rhs( s, 0 );
        const TF dc = coeff_flattening_ratio * pow( current_flattening_ratio, coeff_flattening_ratio - 1 );
        for( PI i = 0; i < s; ++i ) {
            auto ep = [&]() { return pow( extended_positions[ i - 0 ] - extended_positions[ i - 1 ], -2 ); };
            auto en = [&]() { return pow( extended_positions[ i + 1 ] - extended_positions[ i + 0 ], -2 ); };

            TF r = 0;
            if ( i )
                r += ep() * sys_vec[ i - 1 ];
            if ( i + 1 < s )
                r += en() * sys_vec[ i + 1 ];
            
            if ( i >= e && i < n + e )
                r -= ( ep() + en() - d ) * sys_vec[ i ] + d * original_values[ i - e ];
            else if ( i == 0 )
                r -= en() * sys_vec[ i ];
            else if ( i + 1 == s )
                r -= ep() * sys_vec[ i ];
            else
                r -= ( ep() + en() ) * sys_vec[ i ];

            rhs[ i ] = dc * r;
        }

        // munally solve if current_flattening_ratio == 1
        if ( current_flattening_ratio == 1 ) {
            auto ep = [&]( PI i ) { return pow( extended_positions[ i - 0 ] - extended_positions[ i - 1 ], -2 ); };
            auto en = [&]( PI i ) { return pow( extended_positions[ i + 1 ] - extended_positions[ i + 0 ], -2 ); };
            TF v0 = 0, v1 = v0 - rhs[ 0 ] / en( 0 );
            der_1.resize( s );
            der_1[ 0 ] = v0;
            der_1[ 1 ] = v1;
            for( PI i = 1; i + 1 < s; ++i ) {
                const TF v2 = ( 
                    + ( ep( i ) + en( i ) ) * v1
                    - ep( i ) * v0
                    - rhs[ i ]
                ) / ( en( i ) );

                der_1[ i + 1 ] = v2;

                v0 = v1;
                v1 = v2;
            }

            TF sum = 0;
            for( PI i = 0; i < n; ++i )
                sum += rhs[ e + i ] - der_1[ e + i ];        
            sum /= n;
        
            for( TF &v : der_1 )
                v += sum;
        } else
            der_1 = sys_mat.solve_using_ldlt( rhs );

        // normalization
        TF der_mul = 0;
        if ( normalization ) {
            for( PI i = 0; i + 1 < values.size(); ++i ) {
                const TF x0 = positions[ i + 0 ]; 
                const TF x1 = positions[ i + 1 ]; 
                const TF v0 = der_1[ e + i + 0 ]; 
                const TF v1 = der_1[ e + i + 1 ];
                der_mul += ( x1 - x0 ) * ( v0 + v1 ) / 2;
            }
        }

        //
        VF rst( n );
        for( PI i = 0; i < n; ++i )
            rst[ i ] = der_1[ e + i ] / sys_div - sys_vec[ e + i ] * der_mul / pow( sys_div, 2 );
        der_values.push_back( rst );
        der_primitives.push_back( _primitive_of( rst ) );
    }

    VF der_2;
    if ( nb_derivatives >= 2 ) {
        const TF c2 = ( coeff_flattening_ratio - 1 ) * coeff_flattening_ratio * pow( current_flattening_ratio, coeff_flattening_ratio - 2 );
        const TF c1 = coeff_flattening_ratio * pow( current_flattening_ratio, coeff_flattening_ratio - 1 );
        // const TF c0 = pow( current_flattening_ratio, coeff_flattening_ratio );
        VF rhs( s, 0 );
        for( PI i = 0; i < s; ++i ) {
            auto ep = [&]() { return pow( extended_positions[ i - 0 ] - extended_positions[ i - 1 ], -2 ); };
            auto en = [&]() { return pow( extended_positions[ i + 1 ] - extended_positions[ i + 0 ], -2 ); };

            TF r = 0;
            if ( i )
                r += c2 * ep() * sys_vec[ i - 1 ] + 2 * c1 * ep() * der_1[ i - 1 ];
            if ( i + 1 < s )
                r += c2 * en() * sys_vec[ i + 1 ] + 2 * c1 * en() * der_1[ i + 1 ];
            
            if ( i >= e && i < n + e )
                r -= c2 * ( ( ep() + en() - d ) * sys_vec[ i ] + d * original_values[ i - e ] ) + 2 * c1 * ( ep() + en() - d ) * der_1[ i ];
            else if ( i == 0 )
                r -= c2 * en() * sys_vec[ i ] + 2 * c1 * en() * der_1[ i ];
            else if ( i + 1 == s )
                r -= c2 * ep() * sys_vec[ i ] + 2 * c1 * ep() * der_1[ i ];
            else
                r -= c2 * ( ep() + en() ) * sys_vec[ i ] + 2 * c1 * ( ep() + en() ) * der_1[ i ];

            rhs[ i ] = r;
        }

        // manually solve if current_flattening_ratio == 1
        if ( current_flattening_ratio == 1 ) {
            throw std::runtime_error( "TODO der 1" );

            // auto ep = [&]( PI i ) { return pow( extended_positions[ i - 0 ] - extended_positions[ i - 1 ], -2 ); };
            // auto en = [&]( PI i ) { return pow( extended_positions[ i + 1 ] - extended_positions[ i + 0 ], -2 ); };
            // TF v0 = 0, v1 = v0 - rhs[ 0 ] / en( 0 );
            // der.resize( s );
            // der[ 0 ] = v0;
            // der[ 1 ] = v1;
            // for( PI i = 1; i + 1 < s; ++i ) {
            //     const TF v2 = ( 
            //         + ( ep( i ) + en( i ) ) * v1
            //         - ep( i ) * v0
            //         - rhs[ i ]
            //     ) / ( en( i ) );

            //     der[ i + 1 ] = v2;

            //     v0 = v1;
            //     v1 = v2;
            // }

            // TF sum = 0;
            // for( PI i = 0; i < n; ++i )
            //     sum += rhs[ e + i ] - der[ e + i ];        
            // sum /= n;
        
            // for( TF &v : der )
            //     v += sum;
        } else
            der_2 = sys_mat.solve_using_ldlt( rhs );

        // normalization
        TF dm1 = 0;
        TF dm2 = 0;
        if ( normalization ) {
            for( PI i = 0; i + 1 < values.size(); ++i ) {
                const TF x0 = positions[ i + 0 ]; 
                const TF x1 = positions[ i + 1 ]; 
                dm1 += ( x1 - x0 ) * ( der_1[ e + i + 0 ] + der_1[ e + i + 1 ] ) / 2;
                dm2 += ( x1 - x0 ) * ( der_2[ e + i + 0 ] + der_2[ e + i + 1 ] ) / 2;
            }
        }

        //
        VF rst( n );
        // d'' = V'' * sys_mul
        //     - 2 * V' * sum( x' ) / sum( x )^2
        //     - V * sum( x'' ) / sum( x )^2
        //     + 2 * V * sum( x' )^2 / sum( x )^3
        for( PI i = 0; i < n; ++i )
            rst[ i ] = der_2[ e + i ] / sys_div
                     - 2 * der_1[ e + i ] * dm1 / pow( sys_div, 2 )
                     + 2 * sys_vec[ e + i ] * pow( dm1, 2 ) / pow( sys_div, 3 )
                     - sys_vec[ e + i ] * dm2 / pow( sys_div, 2 )
                     ;
        der_values.push_back( rst );
        der_primitives.push_back( _primitive_of( rst ) );
    }
}

DTP TF UTP::x2_integral( TF b, TF e, TF di ) const {
    // TODO: optimize
    const auto &positions = indexed_positions.positions;
    using namespace std;
    TF res = 0;
    for( PI i = 1; i < positions.size(); ++i ) {
        const TF x0 = positions[ i - 1 ];
        const TF x1 = positions[ i - 0 ];
        const TF i0 = max( x0, b );
        const TF i1 = min( x1, e );
        if ( i1 <= i0 )
            continue;
        const TF v0 = values[ i - 1 ];
        const TF v1 = values[ i - 0 ];
        // int( ( x - di ) ^ 2 * ( v0 + ( x - x0 ) / ( x1 - x0 ) * ( v1 - v0 ) ) )
        // res += (
        //     i1 * (
        //         + v0 * ( -12 * di * di * x1 - 4 * i1 * i1 * ( 2 * di + x1 ) + 6 * i1 * di * ( di + 2 * x1 ) + 3 * i1 * i1 * i1 )
        //         + v1 * ( +12 * di * di * x0 + 4 * i1 * i1 * ( 2 * di + x0 ) - 6 * i1 * di * ( di + 2 * x0 ) - 3 * i1 * i1 * i1 )
        //     ) - 
        //     i0 * (
        //         + v0 * ( -12 * di * di * x1 - 4 * i0 * i0 * ( 2 * di + x1 ) + 6 * i0 * di * ( di + 2 * x1 ) + 3 * i0 * i0 * i0 )
        //         + v1 * ( +12 * di * di * x0 + 4 * i0 * i0 * ( 2 * di + x0 ) - 6 * i0 * di * ( di + 2 * x0 ) - 3 * i0 * i0 * i0 )
        //     )
        // ) / ( 12 * ( x0 - x1 ) );
        res += (
            pow( i1 - di, 3 ) * (
                - v0 * ( 3 * i1 + di - 4 * x1 )
                + v1 * ( 3 * i1 + di - 4 * x0 )
            ) - 
            pow( i0 - di, 3 ) * (
                - v0 * ( 3 * i0 + di - 4 * x1 )
                + v1 * ( 3 * i0 + di - 4 * x0 )
            )
        ) / ( 12 * ( x1 - x0 ) );
    }
    return res;
}

DTP TF UTP::x_primitive( TF x ) const {
    auto index = indexed_positions.index( x );
    if ( index.f > 1 )
        return x_primitives.back();
    if ( index.f < 0 )
        return 0;
    
    const TF dx = index.x1 - index.x0;
    const PI i0 = index.i - 1;
    const PI i1 = index.i - 0;

    const TF v0 = values[ i0 ];
    const TF v1 = values[ i1 ];
    const TF dv = v1 - v0;

    const TF f1 = index.f;
    const TF f2 = f1 * f1;
    const TF f3 = f2 * f1;

    return x_primitives[ i0 ] + dx * ( f1 * v0 * index.x0 + f3 * dv * dx / 3 + ( f2 * ( v1 * index.x0 + v0 * index.x1 - 2 * v0 * index.x0 ) ) / 2 );
}

DTP void UTP::disk_integral( TF x0, TF x1, TF c, TF r, TF &mass, TF &dmass ) const {
    using namespace std;

    const auto &positions = indexed_positions.positions;
    auto i0 = indexed_positions.bounded_index( x0 );
    auto i1 = indexed_positions.bounded_index( x1 );
    if ( i0.n == i1.n && i0.f == i1.f )
        return;

    /* sage:
        r, x, c, v0, v1, p0, p1 = var('r,x,c,v0,v1,p0,p1')
        assume( r > 0 )
        assume( x > 0 )
        assume( c > 0 )
        assume( x - c > 0 )
        assume( p1 - p0 > 0 )
        assume( r^2 - ( x - c )^2 > 0 )
        integral( 2 * ( r^2 - ( x - c )^2 )^(1/2) * ( v0 + ( v1 - v0 ) * ( x - p0 ) / ( p1 - p0 ) ), x )
    */
    auto p = [&]( auto x, const TF p0, const TF p1, const TF v0, const TF v1 ) -> TF {
        const TF dp = p0 - p1;
        if ( dp == 0 )
            return 0;
        const TF as = asin( clamp( ( c - x ) / r, TF( -1 ), TF( 1 ) ) );
        const TF x2 = pow( x, 2 );
        const TF r2 = pow( r, 2 );
        const TF c2 = pow( c, 2 );
        const TF c3 = pow( c, 3 );
        const TF s2 = max( TF( 0 ), r2 - pow( c - x, 2 ) );
        const TF s3 = pow( s2, TF( 3 ) / 2 );
        const TF s1 = sqrt( s2 );
        const TF dv = v0 - v1;
        const TF d0 = p0 * dv / dp - v0;
        return 
            + (
                + s1 * ( c * x - c2 )
                - as * r2 * c
                - 2 * s3 / 3
            ) * dv / dp 
             
            + (
                + s1 * ( c - x )
                + r2 * as
            ) * d0
            ;
    };

    auto contrib = [&]( const TF p0, const TF p1, const TF v0, const TF v1 ) {
        if ( p0 == p1 )
            return;

        mass += p( p1, p0, p1, v0, v1 ) - p( p0, p0, p1, v0, v1 );

        const TF a0 = acos( clamp( TF( ( p0 - c ) / r ), TF( -1 ), TF( +1 ) ) );
        const TF a1 = acos( clamp( TF( ( p1 - c ) / r ), TF( -1 ), TF( +1 ) ) );
        const TF dc = cos( a1 ) - cos( a0 );
        const TF ds = sin( a1 ) - sin( a0 );
        dmass -= (
            ds * ( v1 - v0 ) +
            ( a1 - a0 ) * ( cos( a1 ) * v0 - cos( a0 ) * v1 )
        ) / dc;
    };

    if ( i0.n < i1.n ) {
        TF p0 = positions[ i0.n - 1 ] * ( 1 - i0.f ) + positions[ i0.n - 0 ] * i0.f;
        TF v0 = values[ i0.n - 1 ] * ( 1 - i0.f ) + values[ i0.n - 0 ] * i0.f;
        TF acc = 0;
        for( PI i = i0.n; i < i1.n; ++i ) {
            const TF pi = positions[ i ];
            const TF vi = values[ i ];

            contrib( p0, pi, v0, vi );

            v0 = vi;
            p0 = pi;
        }

        const TF p1 = positions[ i1.n - 1 ] * ( 1 - i1.f ) + positions[ i1.n - 0 ] * i1.f;
        const TF v1 = values[ i1.n - 1 ] * ( 1 - i1.f ) + values[ i1.n - 0 ] * i1.f;
        contrib( p0, p1, v0, v1 );
    } else {
        const TF p0 = positions[ i0.n - 1 ] * ( 1 - i0.f ) + positions[ i0.n - 0 ] * i0.f;
        const TF p1 = positions[ i0.n - 1 ] * ( 1 - i1.f ) + positions[ i0.n - 0 ] * i1.f;
        const TF v0 = values[ i0.n - 1 ] * ( 1 - i0.f ) + values[ i0.n - 0 ] * i0.f;
        const TF v1 = values[ i0.n - 1 ] * ( 1 - i1.f ) + values[ i0.n - 0 ] * i1.f;
        contrib( p0, p1, v0, v1 );
    }

}

DTP TF UTP::disk_integral( TF x0, TF x1, TF c, TF r ) const {
    using namespace std;

    const auto &positions = indexed_positions.positions;
    auto i0 = indexed_positions.bounded_index( x0 );
    auto i1 = indexed_positions.bounded_index( x1 );
    if ( i0.n == i1.n && i0.f == i1.f )
        return 0;

    /* sage:
        r, x, c, v0, v1, p0, p1 = var('r,x,c,v0,v1,p0,p1')
        assume( r > 0 )
        assume( x > 0 )
        assume( c > 0 )
        assume( x - c > 0 )
        assume( p1 - p0 > 0 )
        assume( r^2 - ( x - c )^2 > 0 )
        integral( 2 * ( r^2 - ( x - c )^2 )^(1/2) * ( v0 + ( v1 - v0 ) * ( x - p0 ) / ( p1 - p0 ) ), x )
    */
    auto p = [&]( auto x, const TF p0, const TF p1, const TF v0, const TF v1 ) -> TF {
        const TF dp = p0 - p1;
        if ( dp == 0 )
            return 0;
        const TF as = asin( clamp( ( c - x ) / r, TF( -1 ), TF( 1 ) ) );
        const TF x2 = pow( x, 2 );
        const TF r2 = pow( r, 2 );
        const TF c2 = pow( c, 2 );
        const TF c3 = pow( c, 3 );
        const TF s2 = max( TF( 0 ), r2 - pow( c - x, 2 ) );
        const TF s3 = pow( s2, TF( 3 ) / 2 );
        const TF s1 = sqrt( s2 );
        const TF dv = v0 - v1;
        const TF d0 = p0 * dv / dp - v0;
        return 
            + (
                + s1 * ( c * x - c2 )
                - as * r2 * c
                - 2 * s3 / 3
            ) * dv / dp 
             
            + (
                + s1 * ( c - x )
                + r2 * as
            ) * d0
            ;
    };

    auto pa = [&]( const TF p0, const TF p1, const TF v0, const TF v1 ) -> TF {
        return p( p1, p0, p1, v0, v1 ) - p( p0, p0, p1, v0, v1 );
    };

    if ( i0.n < i1.n ) {
        TF p0 = positions[ i0.n - 1 ] * ( 1 - i0.f ) + positions[ i0.n - 0 ] * i0.f;
        TF v0 = values[ i0.n - 1 ] * ( 1 - i0.f ) + values[ i0.n - 0 ] * i0.f;
        TF acc = 0;
        for( PI i = i0.n; i < i1.n; ++i ) {
            const TF pi = positions[ i ];
            const TF vi = values[ i ];

            acc += pa( p0, pi, v0, vi );

            v0 = vi;
            p0 = pi;
        }

        const TF p1 = positions[ i1.n - 1 ] * ( 1 - i1.f ) + positions[ i1.n - 0 ] * i1.f;
        const TF v1 = values[ i1.n - 1 ] * ( 1 - i1.f ) + values[ i1.n - 0 ] * i1.f;
        acc += pa( p0, p1, v0, v1 );

        return acc;
    }

    const TF p0 = positions[ i0.n - 1 ] * ( 1 - i0.f ) + positions[ i0.n - 0 ] * i0.f;
    const TF p1 = positions[ i0.n - 1 ] * ( 1 - i1.f ) + positions[ i0.n - 0 ] * i1.f;
    const TF v0 = values[ i0.n - 1 ] * ( 1 - i0.f ) + values[ i0.n - 0 ] * i0.f;
    const TF v1 = values[ i0.n - 1 ] * ( 1 - i1.f ) + values[ i0.n - 0 ] * i1.f;
    return pa( p0, p1, v0, v1 );
}

DTP TF UTP::primitive( TF x ) const {
    auto index = indexed_positions.index( x );
    if ( index.f > 1 )
        return primitives.back();
    if ( index.f < 0 )
        return 0;

    const TF v0 = values[ index.n - 1 ];
    const TF v1 = values[ index.n - 0 ];
    const TF xd = x - index.x0;
    
    if ( const TF dx = index.x1 - index.x0 )
        return primitives[ index.n - 1 ] + xd * v0 + xd * xd / dx * ( v1 - v0 ) / 2;
    return primitives[ index.n - 1 ];
}

DTP TF UTP::value( TF x ) const {
    auto index = indexed_positions.index( x );
    if ( index.f > 1 || index.f < 0 )
        return 0;

    const TF v0 = values[ index.n - 1 ];
    const TF v1 = values[ index.n - 0 ];
    return v0 * ( 1 - index.f ) + v1 * index.f;
}

DTP TF UTP::value( TF x, PI num_der ) const {
    auto index = indexed_positions.index( x );
    if ( index.f > 1 || index.f < 0 )
        return 0;

    const VF &values = num_der ? this->der_values[ num_der - 1 ] : this->values;
    const TF v0 = values[ index.n - 1 ];
    const TF v1 = values[ index.n - 0 ];
    return v0 * ( 1 - index.f ) + v1 * index.f;
}

DTP TF UTP::primitive( TF x, PI num_der ) const {
    const VF &primitives = num_der ? der_primitives[ num_der - 1 ] : this->primitives;
    const VF &values = num_der ? der_values[ num_der - 1 ] : this->values;

    auto index = indexed_positions.index( x );
    if ( index.f > 1 )
        return primitives.back();
    if ( index.f < 0 )
        return 0;

    if ( const TF dx = index.x1 - index.x0 ) {
        const TF v0 = values[ index.n - 1 ];
        const TF v1 = values[ index.n - 0 ];
        const TF xd = x - index.x0;
        
        TF res = xd * v0 + xd * xd / dx * ( v1 - v0 ) / 2;
        return primitives[ index.n - 1 ] + res;
    }
    return primitives[ index.n - 1 ];
}

DTP typename UTP::VF UTP::_primitive_of( const VF &values ) const {
    const auto &positions = indexed_positions.positions;
    const PI n = values.size();
    VF res( n );

    TF o = 0;
    for( PI i = 0; i + 1 < n; ++i ) {
        const TF x0 = positions[ i + 0 ]; 
        const TF x1 = positions[ i + 1 ]; 
        const TF v0 = values[ i + 0 ]; 
        const TF v1 = values[ i + 1 ];

        res[ i ] = o;

        o += ( x1 - x0 ) * ( v0 + v1 ) / 2;
    }
    res.back() = o;

    return res;
}

DTP void UTP::_compute_values_for( TF flattening_ratio ) {
    using namespace std;

    const auto &positions = indexed_positions.positions;
    const PI s = extended_positions.size();
    const PI n = positions.size();
    const PI e = ( s - n ) / 2;

    //
    der_values.clear();

    // specific values
    if ( flattening_ratio == 0 ) {
        values = original_values;
        return;
    }

    if ( flattening_ratio == 1 ) {
        const TF v = 1 / ( ptp_x() );

        values.resize( n );
        for( PI i = 0; i < n; ++i )
            values[ i ] = v;

        sys_vec.resize( s );
        for( PI i = 0; i < s; ++i )
            sys_vec[ i ] = v;

        return;
    }

    //
    auto ep = [&]( PI i ) { return pow( extended_positions[ i - 0 ] - extended_positions[ i - 1 ], -2 ); };
    auto en = [&]( PI i ) { return pow( extended_positions[ i + 1 ] - extended_positions[ i + 0 ], -2 ); };
    const TF c = pow( flattening_ratio, coeff_flattening_ratio );
    const TF d = pow( ptp_x(), -2 );
    sys_mat = { n + 2 * e };
    for( PI i = 0; i < n + 2 * e; ++i ) {

        if ( i )
            sys_mat( i, i - 1 ) = - c * ep( i );
        
        if ( i == 0 )
            sys_mat( i, i ) = c * en( i );
        else if ( i + 1 == n + 2 * e )
            sys_mat( i, i ) = c * ep( i );
        else if ( i >= e && i < n + e )
            sys_mat( i, i ) = c * ( ep( i ) + en( i ) ) + ( 1 - c ) * d;
        else
            sys_mat( i, i ) = c * ( ep( i ) + en( i ) );
    }

    VF E( n + 2 * e, 0 );
    for( PI i = 0; i < n; ++i )
        E[ e + i ] = ( 1 - c ) * d * original_values[ i ];

    sys_mat.inplace_ldlt_decomposition();
    sys_vec = sys_mat.solve_using_ldlt( E );

    values.resize( n );
    for( PI i = 0; i < n; ++i )
        values[ i ] = sys_vec[ e + i ];
}

DTP TF UTP::integral( TF x0, TF x1, PI num_der_lag_ratio ) const {
    return primitive( x1, num_der_lag_ratio ) - primitive( x0, num_der_lag_ratio );
}


DTP TF UTP::x_integral( TF x0, TF x1 ) const {
    return x_primitive( x1 ) - x_primitive( x0 );
}

DTP TF UTP::integral( TF x0, TF x1 ) const {
    return primitive( x1 ) - primitive( x0 );
}

DTP TF UTP::min_x() const {
    return indexed_positions.min_x();
}

DTP TF UTP::max_x() const {
    return indexed_positions.max_x();
}

DTP TF UTP::ptp_x() const {
    return max_x() - min_x();
}

DTP void UTP::plot( std::ostream &fs, std::string linestyle, double linewidth ) const {
    const auto &positions = indexed_positions.positions;
    fs << "pyplot.plot( ["; 
    for( PI n = 0; n < positions.size(); ++n )
        fs << positions[ n ] << ", ";    
    fs << " ], [";
    for( PI n = 0; n < values.size(); ++n )
        fs << values[ n ] << ", ";    
    fs << "], '" << linestyle << "', linewidth = " << linewidth << " )\n";
}

#undef DTP
#undef UTP

} // namespace usdot
