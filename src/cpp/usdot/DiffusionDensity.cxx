#pragma once

#include "utility/linspace.h"
#include "DiffusionDensity.h"
#include <algorithm>
#include <stdexcept>
#include <limits>

namespace usdot {
    
#define DTP template<class TF>
#define UTP DiffusionDensity<TF>

DTP UTP::DiffusionDensity( TF beg_original_positions, TF end_original_positions, const VF &original_values ) : DiffusionDensity( linspace( beg_original_positions, end_original_positions, original_values.size() ), original_values ) {
}

DTP UTP::DiffusionDensity( const VF &original_values ) : DiffusionDensity( 0, 1, original_values ) {
}

DTP UTP::DiffusionDensity( const VF &positions, const VF &original_values ) {
    this->current_flattening_ratio = -1;
    
    this->original_values = original_values;
    this->positions = positions;

    // normalization
    TF o = 0;
    for( PI i = 0; i + 1 < this->original_values.size(); ++i )
        o += ( positions[ i + 1 ] - positions[ i + 0 ] ) * ( original_values[ i + 0 ] + original_values[ i + 1 ] ) / 2;
    for( auto &v : this->original_values )
        v /= o;

    // opt_pos
    const TF dx = positions.back() - positions.front();
    const PI opt_size = positions.size();
    opt_pos_beg_inds.resize( opt_size );
    opt_pos_beg_x = positions.front();
    opt_pos_mul_x = opt_size / dx;

    for( PI &ind : opt_pos_beg_inds )
        ind = opt_size;
    for( PI i = 1; i < positions.size(); ++i ) {
        const PI x0 = PI( floor( ( positions[ i - 1 ] - opt_pos_beg_x ) * opt_pos_mul_x ) );
        const PI x1 = PI( ceil( ( positions[ i - 0 ] - opt_pos_beg_x ) * opt_pos_mul_x ) );
        for( PI n = x0; n < std::min( x1, opt_size ); ++n )
            opt_pos_beg_inds[ n ] = std::min( opt_pos_beg_inds[ n ], i );
    }
}

DTP void UTP::set_flattening_ratio( TF flattening_ratio ) {
    using namespace std;

    if ( current_flattening_ratio == flattening_ratio )
        return;
    current_flattening_ratio = flattening_ratio;

    //
    _compute_values_for( flattening_ratio );

    // primitives
    x_primitives.resize( this->values.size() );
    primitives.resize( this->values.size() );
    TF o = 0, x_o = 0;
    for( PI i = 0; i + 1 < this->values.size(); ++i ) {
        const TF v0 = this->values[ i + 0 ]; 
        const TF v1 = this->values[ i + 1 ];
        const TF x0 = positions[ i + 0 ]; 
        const TF x1 = positions[ i + 1 ]; 
        const TF x6 = ( x0 + x1 ) / 6;

        x_primitives[ i ] = x_o;
        primitives[ i ] = o;

        // ( x1 - x0 ) * ( f * v0 * x0 + (f^3 * (v0 - v1) * (x0 - x1) ) / 3 + ( f^2 * (-2 * v0 * x0 + v1 * x0 + v0 * x1 ) ) / 2 )
        x_o += v0 * ( x1 * x6 - x0 * x0 / 3 ) + v1 * ( x1 * x1 / 3 - x0 * x6 );
        o += ( x1 - x0 ) * ( v0 + v1 ) / 2;
    }

    x_primitives.back() = x_o;
    primitives.back() = o;

    // normalization (actually, o should be equal to 1 by construction)
    sys_div = o;
    for( TF &v : this->x_primitives )
        v /= sys_div;
    for( TF &v : this->primitives )
        v /= sys_div;
    for( TF &v : this->values )
        v /= sys_div;
}

DTP void UTP::compute_derivatives( PI nb_derivatives ) {
    const PI n = positions.size();
    const PI s = sys_vec.size();
    const PI e = ( s - n ) / 2;

    der_primitives.clear();
    der_values.clear();

    // 0
    if ( current_flattening_ratio == 0 ) {
        throw std::runtime_error( "TODO" );
    }

    // const TF c = flattening_ratio * pow( n, 2 );
    // sys_mat = { n + 2 * e };
    // for( PI i = 0; i < n + 2 * e; ++i ) {
    //     if ( i )
    //         sys_mat( i, i - 1 ) = - c;
        
    //     if ( i == 0 || i + 1 == n + 2 * e )
    //         sys_mat( i, i ) = c;
    //     else if ( i >= e && i < n + e )
    //         sys_mat( i, i ) = 2 * c + 1 - flattening_ratio;
    //     else
    //         sys_mat( i, i ) = 2 * c;
    // }

    // VF E( n + 2 * e, 0 );
    // for( PI i = 0; i < n; ++i )
    //     E[ e + i ] = ( 1 - flattening_ratio ) * original_values[ i ];

    // 
    VF rhs( s, 0 );
    const TF c = pow( n, 2 );
    for( PI i = 0; i < s; ++i ) {
        TF r = 0;
        if ( i )
            r += c * sys_vec[ i - 1 ];
        if ( i + 1 < s )
            r += c * sys_vec[ i + 1 ];
        
        if ( i >= e && i < n + e )
            r -= ( 2 * c - 1 ) * sys_vec[ i ] + original_values[ i - e ];
        else if ( i && i + 1 < s )
            r -= 2 * c * sys_vec[ i ];
        else
            r -= c * sys_vec[ i ];

        rhs[ i ] = r;
    }

    // en 1
    VF der;
    if ( current_flattening_ratio == 1 ) {
        TF v0 = 0, v1 = v0;
        der.resize( s );
        der[ 0 ] = v0;
        der[ 1 ] = v1;
        for( PI i = 2; i < s; ++i ) {
            const TF v2 = 2 * v1 - v0 - rhs[ i - 1 ] / c;

            der[ i ] = v2;

            v0 = v1;
            v1 = v2;
        }

        TF sum = 0;
        for( PI i = 0; i < n; ++i )
            sum += rhs[ e + i ] - der[ e + i ];        
        sum /= n;
    
        for( TF &v : der )
            v += sum;
    } else
        der = sys_mat.solve_using_ldlt( rhs );

    // normalization
    TF der_mul = 0;
    for( PI i = 0; i + 1 < values.size(); ++i ) {
        const TF x0 = positions[ i + 0 ]; 
        const TF x1 = positions[ i + 1 ]; 
        const TF v0 = der[ e + i + 0 ]; 
        const TF v1 = der[ e + i + 1 ];
        der_mul += ( x1 - x0 ) * ( v0 + v1 ) / 2;
    }

    //
    VF rst( n );
    for( PI i = 0; i < n; ++i )
        rst[ i ] = der[ e + i ] / sys_div - sys_vec[ e + i ] * der_mul / pow( sys_div, 2 );
    der_values.push_back( rst );
    der_primitives.push_back( _primitive_of( rst ) );
}

DTP TF UTP::x_primitive( TF x ) const {
    // if ( x < opt_pos_beg )
    //     return 0;
    // if ( x > opt_pos_end )
    //     return x_primitives.back();

    // const PI ox = std::min( PI( ( x - opt_pos_beg ) * opt_pos_coeff ), opt_pos_begs.size() - 1 );
    // for( PI i = opt_pos_begs[ ox ]; ; ++i ) {
    //     if ( positions[ i + 1 ] >= x ) {
    //         const TF p0 = positions[ i + 0 ];
    //         const TF p1 = positions[ i + 1 ];
    //         if ( const TF dp = p1 - p0 ) {
    //             const TF v0 = values[ i + 0 ];
    //             const TF v1 = values[ i + 1 ];
    //             const TF dv = v1 - v0;
    //             const TF f = ( x - p0 ) / dp;
    //             const TF f2 = f * f;
    //             const TF f3 = f2 * f;
        
    //             return x_primitives[ i ] + dp * ( f * v0 * p0 + f3 * dv * dp / 3 + ( f2 * ( v1 * p0 + v0 * p1 - 2 * v0 * p0 ) ) / 2 );
    //         }
    //         return x_primitives[ i ];
    //     }
    // }
    throw std::runtime_error( "TODO" );

    return 0;
}

DTP TF UTP::primitive( TF x ) const {
    if ( x < opt_pos_beg_x )
        return 0;

    const PI ox = PI( ( x - opt_pos_beg_x ) * opt_pos_mul_x );
    if ( ox >= opt_pos_beg_inds.size() )
        return primitives.back();

    for( PI i = opt_pos_beg_inds[ ox ]; i < positions.size(); ++i ) {
        if ( positions[ i ] >= x ) {
            const TF x0 = positions[ i - 1 ];
            const TF x1 = positions[ i - 0 ];
            if ( const TF dx = x1 - x0 ) {
                const TF v0 = values[ i - 1 ];
                const TF v1 = values[ i - 0 ];
                const TF xd = x - x0;
              
                TF res = xd * v0 + xd * xd / dx * ( v1 - v0 ) / 2;
                return primitives[ i - 1 ] + res;
            }
            return primitives[ i - 1 ];
        }
    }

    return primitives.back();
}

DTP TF UTP::value( TF x ) const {
    if ( x < opt_pos_beg_x )
        return 0;

    const PI ox = PI( ( x - opt_pos_beg_x ) * opt_pos_mul_x );
    if ( ox >= opt_pos_beg_inds.size() )
        return 0;

    for( PI i = opt_pos_beg_inds[ ox ]; i < positions.size(); ++i ) {
        if ( positions[ i ] >= x ) {
            if ( const TF de = positions[ i - 0 ] - positions[ i - 1 ] ) {
                TF f = ( x - positions[ i - 1 ] ) / de;
                return values[ i - 1 ] * ( 1 - f ) + values[ i - 0 ] * f;
            }
            return values[ i ];
        }
    }

    return 0;
}

DTP TF UTP::value( TF x, PI num_der ) const {
    // 
    if ( num_der == 0 )
        return value( x );
    const VF &values = this->der_values[ num_der - 1 ];

    //
    if ( x < opt_pos_beg_x )
        return 0;

    const PI ox = PI( ( x - opt_pos_beg_x ) * opt_pos_mul_x );
    if ( ox >= opt_pos_beg_inds.size() )
        return 0;

    for( PI i = opt_pos_beg_inds[ ox ]; i < positions.size(); ++i ) {
        if ( positions[ i ] >= x ) {
            if ( const TF de = positions[ i - 0 ] - positions[ i - 1 ] ) {
                TF f = ( x - positions[ i - 1 ] ) / de;
                return values[ i - 1 ] * ( 1 - f ) + values[ i - 0 ] * f;
            }
            return values[ i ];
        }
    }

    return 0;
}

DTP UTP::VF UTP::_primitive_of( const VF &values ) const {
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

    const PI n = original_values.size();
    const PI e = n;
    const PI s = n + 2 * e;

    //
    der_values.clear();

    // specific values
    if ( flattening_ratio == 0 ) {
        values = original_values;
        return;
    }

    if ( flattening_ratio == 1 ) {
        const TF v = ptp_x();

        values.resize( n );
        for( PI i = 0; i < n; ++i )
            values[ i ] = v;

        sys_vec.resize( s );
        for( PI i = 0; i < s; ++i )
            sys_vec[ i ] = v;

        return;
    }

    //
    const TF c = flattening_ratio * pow( n, 2 );
    sys_mat = { n + 2 * e };
    for( PI i = 0; i < n + 2 * e; ++i ) {
        if ( i )
            sys_mat( i, i - 1 ) = - c;
        
        if ( i == 0 || i + 1 == n + 2 * e )
            sys_mat( i, i ) = c;
        else if ( i >= e && i < n + e )
            sys_mat( i, i ) = 2 * c + 1 - flattening_ratio;
        else
            sys_mat( i, i ) = 2 * c;
    }

    VF E( n + 2 * e, 0 );
    for( PI i = 0; i < n; ++i )
        E[ e + i ] = ( 1 - flattening_ratio ) * original_values[ i ];

    sys_mat.inplace_ldlt_decomposition();
    sys_vec = sys_mat.solve_using_ldlt( E );

    values.resize( n );
    for( PI i = 0; i < n; ++i )
        values[ i ] = sys_vec[ e + i ];
}

// DTP TF UTP::value( TF x, PI num_der_lag_ratio ) {
//     // not a derivative ? 
//     if ( num_der_lag_ratio == 0 )
//         return value( x );

//     // precomputation(s)
//     while ( der_values.size() < num_der_lag_ratio )
//         _append_der_value();

//     //
//     const VF &values = der_values[ num_der_lag_ratio - 1 ];
//     if ( regular_positions ) {
//         if ( regular_positions == 1 )
//             x = ( x - opt_pos_beg ) * opt_pos_coeff;
//         if ( x < 0 || x >= values.size() - 1 )
//             return 0;
    
//         PI i( x );
//         TF f = x - i;
//         return values[ i ] * ( 1 - f ) + values[ i + 1 ] * f;
//     }

//     if ( x < opt_pos_beg || x > opt_pos_end )
//         return 0;

//     const PI ox = std::min( PI( ( x - opt_pos_beg ) * opt_pos_coeff ), opt_pos_begs.size() - 1 );
//     for( PI i = opt_pos_begs[ ox ]; ; ++i ) {
//         if ( positions[ i + 1 ] >= x ) {
//             if ( const TF de = positions[ i + 1 ] - positions[ i + 0 ] ) {
//                 TF f = ( x - positions[ i + 0 ] ) / de;
//                 return values[ i + 0 ] * ( 1 - f ) + values[ i + 1 ] * f;
//             }
//             return values[ i ];
//         }
//     }
// }

// DTP TF UTP::primitive( TF x, PI num_der_lag_ratio ) {
//     if ( num_der_lag_ratio == 0 )
//         return primitive( x );

//     // precomputation(s)
//     while ( der_values.size() < num_der_lag_ratio )
//         _append_der_value();

//     const auto &primitives = der_primitives[ num_der_lag_ratio - 1 ];
//     const auto &values = der_values[ num_der_lag_ratio - 1 ];
//     if ( regular_positions ) {
//         if ( regular_positions == 1 )
//             x = ( x - opt_pos_beg ) * opt_pos_coeff;
//         if ( x < 0 )
//             return 0;

//         const PI i( x );
//         if ( i >= primitives.size() - 1 )
//             return primitives.back();
        
//         TF f = x - i;
//         const TF v0 = values[ i + 0 ];
//         const TF v1 = values[ i + 1 ];
//         TF res = f * v0 + f * f * ( v1 - v0 ) / 2;
//         if ( regular_positions == 1 )
//             res /= opt_pos_coeff;
//         return primitives[ i ] + res;
//     }

//     if ( x < opt_pos_beg )
//         return 0;
//     if ( x > opt_pos_end )
//         return primitives.back();

//     const PI ox = std::min( PI( ( x - opt_pos_beg ) * opt_pos_coeff ), opt_pos_begs.size() - 1 );
//     for( PI i = opt_pos_begs[ ox ]; ; ++i ) {
//         if ( positions[ i + 1 ] >= x ) {
//             const TF x0 = positions[ i + 0 ];
//             const TF x1 = positions[ i + 1 ];
//             if ( const TF dx = x1 - x0 ) {
//                 const TF v0 = values[ i + 0 ];
//                 const TF v1 = values[ i + 1 ];
//                 const TF xd = x - x0;
                
//                 TF res = xd * v0 + xd * xd / dx * ( v1 - v0 ) / 2;
//                 return primitives[ i ] + res;
//             }
//             return primitives[ i ];
//         }
//     }
// }

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
    return positions.front();
}

DTP TF UTP::max_x() const {
    return positions.back();
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
