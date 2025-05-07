#pragma once

#include "utility/linspace.h"
#include "StretchedDensity.h"
#include <algorithm>
#include <stdexcept>
#include <limits>

namespace usdot {
    
#define DTP template<class TF>
#define UTP StretchedDensity<TF>

DTP UTP::StretchedDensity( TF beg_original_positions, TF end_original_positions, const VF &original_values ) : StretchedDensity( linspace( beg_original_positions, end_original_positions, original_values.size() ), original_values ) {
}

DTP UTP::StretchedDensity( const VF &original_values ) : StretchedDensity( 0, 1, original_values ) {
}

DTP UTP::StretchedDensity( const VF &original_positions_, const VF &original_values_ ) {
    original_positions = original_positions_;
    original_values = original_values_;

    // normalize
    TF itg = 0;
    for( PI i = 1; i < original_positions.size(); ++i ) {
        const TF x0 = original_positions[ i - 1 ];
        const TF x1 = original_positions[ i - 0 ];
        const TF y0 = original_values[ i - 1 ];
        const TF y1 = original_values[ i - 0 ];
        itg += ( x1 - x0 ) * ( y0 + y1 ) / 2;
    }
    for( TF &v : this->original_values )
        v /= itg;

    // flat_value
    const TF dx = original_positions.back() - original_positions.front();
    flat_value = 1 / dx;
    
    // flat_positions
    // flat_positions
    TF op = original_positions.front();
    flat_positions.push_back( op );
    for( PI i = 1; i < original_positions.size(); ++i ) {
        const TF x0 = original_positions[ i - 1 ];
        const TF x1 = original_positions[ i - 0 ];
        const TF y0 = original_values[ i - 1 ];
        const TF y1 = original_values[ i - 0 ];

        op += ( x1 - x0 ) * ( y0 + y1 ) / ( 2 * flat_value );
        flat_positions.push_back( op );
    }

    // opt_pos_beg_x
    const PI opt_size = original_positions.size();
    opt_pos_beg_x = original_positions.front();
    opt_pos_beg_inds.resize( opt_size );
    opt_pos_mul_x = opt_size / dx;

    // x_primitives
}

DTP void UTP::set_flattening_ratio( TF flattening_ratio ) {
    this->flattening_ratio = flattening_ratio;

    const PI n = original_positions.size();

    // new positions and values
    positions.clear();
    values.clear();

    TF v0 = ( 1 - flattening_ratio ) * original_values.front() + flattening_ratio * flat_value;
    TF op = original_positions.front();
    positions.push_back( op );
    values.push_back( v0 );
    for( PI i = 1; i + 1 < original_positions.size(); ++i ) {
        const TF ox0 = original_positions[ i - 1 ];
        const TF ox1 = original_positions[ i + 0 ];
        const TF ox2 = original_positions[ i + 1 ];
 
        const TF oy = original_values[ i ];

        const TF nx0 = ( 1 - flattening_ratio ) * ox0 + flattening_ratio * flat_positions[ i - 1 ];
        const TF nx1 = ( 1 - flattening_ratio ) * ox1 + flattening_ratio * flat_positions[ i - 0 ];
        const TF nx2 = ( 1 - flattening_ratio ) * ox2 + flattening_ratio * flat_positions[ i + 1 ];

        const TF nv = ( ox2 - ox0 ) / ( nx2 - nx0 ) * oy;

        positions.push_back( nx1 );
        values.push_back( nv );
    }

    // for( PI i = 0; i < n; ++i ) {
    //     const TF p = ( 1 - flattening_ratio ) * original_positions[ i ] + flattening_ratio * flat_positions[ i ];
    //     const TF v = ( 1 - flattening_ratio ) * original_values[ i ] + flattening_ratio * flat_value;
    //     positions[ i ] = p;
    //     values[ i ] = v;
    // }

    //
    primitives = _primitive_of( positions, values );

    //
    for( PI &ind : opt_pos_beg_inds )
        ind = n;
    for( PI i = 1; i  < this->values.size(); ++i ) {
        const PI x0 = floor( ( this->positions[ i - 1 ] - opt_pos_beg_x ) * opt_pos_mul_x );
        const PI x1 = ceil( ( this->positions[ i - 0 ] - opt_pos_beg_x ) * opt_pos_mul_x );
        for( PI n = x0; n < std::min( x1, opt_pos_beg_inds.size() ); ++n )
            opt_pos_beg_inds[ n ] = std::min( opt_pos_beg_inds[ n ], i );
    }
}

// fs << "pyplot.plot( [ ";
// for( auto x : xs )
//     fs << x << ", ";
// fs << " ], [ ";
// for( auto x : xs )
//     fs << func( x ) << ", ";
// fs << " ] )\n";

// DTP TF UTP::x_primitive( TF x ) const {
//     if ( regular_positions ) {
//         if ( regular_positions == 1 )
//             x = ( x - opt_pos_beg ) * opt_pos_coeff;
//         if ( x < 0 )
//             return 0;

//         const PI i( x );
//         if ( i >= x_primitives.size() - 1 )
//             return x_primitives.back();

//         TF p0 = i + 0;
//         TF p1 = i + 1;
//         if ( regular_positions == 1 ) {
//             p0 = opt_pos_beg + p0 / opt_pos_coeff;
//             p1 = opt_pos_beg + p1 / opt_pos_coeff;
//         }
        
//         const TF v0 = values[ i + 0 ];
//         const TF v1 = values[ i + 1 ];
//         const TF dv = v1 - v0;
//         const TF dp = p1 - p0;
//         const TF f = x - i;
//         const TF f2 = f * f;
//         const TF f3 = f2 * f;
//         TF res = ( f * v0 * p0 + f3 * dv * dp / 3 + ( f2 * ( v1 * p0 + v0 * p1 - 2 * v0 * p0 ) ) / 2 );
//         if ( regular_positions == 1 )
//             res /= opt_pos_coeff;
//         return x_primitives[ i ] + res;
//     }

//     if ( x < opt_pos_beg )
//         return 0;
//     if ( x > opt_pos_end )
//         return x_primitives.back();

//     const PI ox = std::min( PI( ( x - opt_pos_beg ) * opt_pos_coeff ), opt_pos_begs.size() - 1 );
//     for( PI i = opt_pos_begs[ ox ]; ; ++i ) {
//         if ( positions[ i + 1 ] >= x ) {
//             const TF p0 = positions[ i + 0 ];
//             const TF p1 = positions[ i + 1 ];
//             if ( const TF dp = p1 - p0 ) {
//                 const TF v0 = values[ i + 0 ];
//                 const TF v1 = values[ i + 1 ];
//                 const TF dv = v1 - v0;
//                 const TF f = ( x - p0 ) / dp;
//                 const TF f2 = f * f;
//                 const TF f3 = f2 * f;
        
//                 return x_primitives[ i ] + dp * ( f * v0 * p0 + f3 * dv * dp / 3 + ( f2 * ( v1 * p0 + v0 * p1 - 2 * v0 * p0 ) ) / 2 );
//             }
//             return x_primitives[ i ];
//         }
//     }
// }

DTP TF UTP::primitive( TF x ) const {
    if ( x < opt_pos_beg_x )
        return 0;

    const PI ox = ( x - opt_pos_beg_x ) * opt_pos_mul_x;
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

    const PI ox = ( x - opt_pos_beg_x ) * opt_pos_mul_x;
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

DTP UTP::VF UTP::_primitive_of( const VF &positions, const VF &values ) const {
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

// DTP void UTP::_set_values( TF t ) {
//     const PI n = original_values.size();
//     const PI p = n - 1;

//     //
//     positions = original_positions;
//     der_values.clear();

//     // values
//     if ( t == 0 ) {
//         values = original_values;
//         return;
//     }

//     if ( t == 1 ) {
//         values.resize( n );
//         for( PI i = 0; i < n; ++i )
//             values[ i ] = 1 / TF( n - 1 );
//         return;
//     }

//     //
//     system = { n };
//     VF E( n );
//     for( PI i = 1; i < p; ++i ) {
//         system.tridiag_value( i, i - 1 ) = - t;
//         system.tridiag_value( i, i ) = 1 + t;
//         system.line_value( i ) = 1;
//     }
//     system.tridiag_value( p, p - 1 ) = - t;
//     system.tridiag_value( p, p ) = 1;
//     system.tridiag_value( 0, 0 ) = 1;
//     system.line_value( 0 ) = 0.5;
//     system.line_value( p ) = 0.5;
//     for( PI i = 0; i < n; ++i )
//         E[ i ] = ( 1 - t ) * original_values[ i ];
//     system.line_value( n ) = 0;

//     system.inplace_ldlt_decomposition();
//     values = system.solve_using_ldlt( E, 1 );
// }

// DTP void UTP::_append_der_value() {
//     const PI n = values.size();
//     const PI p = n - 1;

//     if ( regular_positions != 2 )
//         throw std::runtime_error( "TODO" );

//     // first derivative ?
//     VF E( n );
//     if ( der_values.empty() ) {
//         if ( current_lag_ratio == 1 ) {
//             // [   1,  -1,   0,   0,   0,   0,   0,   0,   0,   0, 0.5 ]
//             // [  -1,   2,  -1,   0,   0,   0,   0,   0,   0,   0,   1 ]
//             // [   0,  -1,   2,  -1,   0,   0,   0,   0,   0,   0,   1 ]
//             // [   0,   0,  -1,   2,  -1,   0,   0,   0,   0,   0,   1 ]
//             // [   0,   0,   0,  -1,   2,  -1,   0,   0,   0,   0,   1 ]
//             // [   0,   0,   0,   0,  -1,   2,  -1,   0,   0,   0,   1 ]
//             // [   0,   0,   0,   0,   0,  -1,   2,  -1,   0,   0,   1 ]
//             // [   0,   0,   0,   0,   0,   0,  -1,   2,  -1,   0,   1 ]
//             // [   0,   0,   0,   0,   0,   0,   0,  -1,   2,  -1,   1 ]
//             // [   0,   0,   0,   0,   0,   0,   0,   0,  -1,   1, 0.5 ]
//             // [ 0.5,   1,   1,   1,   1,   1,   1,   1,   1, 0.5,   0 ]
//             TF ay = 0, by = 0;
//             for( PI i = 0; i < n; ++i ) {
//                 const TF dv = values[ i ] - original_values[ i ];
//                 const TF an = TF( i ) * TF( i ) / 2;
//                 const TF bn = an + 1;
//                 ay += an * dv;
//                 by += bn * dv;
//             }

//             const TF Z_l = TF( p - 1 ) * p * ( 2 * p - 1 ) / 12 + TF( p ) * TF( p ) / 4;
//             const TF lam = ( by - ay ) / p;

//             E[ p ] = ( ay - Z_l * lam ) / p;
//             E[ p - 1 ] = E[ p ] + lam / 2 + original_values[ p ] - values[ p ];
//             for( PI i = p - 1; i--; )
//                 E[ i ] = 2 * E[ i + 1 ] - E[ i + 2 ] + lam + original_values[ i + 1 ] - values[ i + 1 ];

//             der_values.push_back( std::move( E ) );
//             der_primitives.push_back( _primitive_of( der_values.back() ) );
//             return;
//         }

//         for( PI i = 0; i < n; ++i ) {
//             TF e = - original_values[ i ];
//             if ( i && i + 1 < n )
//                 e -= values[ i ];
//             if ( i )
//                 e += values[ i - 1 ];
//             if ( i + 1 < n )
//                 e += values[ i + 1 ];
//             E[ i ] = e;
//         }
//         der_values.push_back( system.solve_using_ldlt( E, 0 ) );
//         der_primitives.push_back( _primitive_of( der_values.back() ) );
//         return;
//     }

//     // 
//     const auto &prev_der_values = der_values.back();
//     const PI nder = der_values.size() + 1;
//     if ( current_lag_ratio == 1 ) {
//         auto diu = [&]( PI i ) {
//             TF dv = 0;
//             if ( i && i + 1 < n )
//                 dv -= prev_der_values[ i ];
//             if ( i )
//                 dv += prev_der_values[ i - 1 ];
//             if ( i + 1 < n )
//                 dv += prev_der_values[ i + 1 ];
//             return nder * dv;
//         };

//         TF ay = 0, by = 0;
//         for( PI i = 0; i < n; ++i ) {
//             const TF an = TF( i ) * TF( i ) / 2;
//             const TF bn = an + 1;
//             const TF dv = diu( i );
//             ay += an * dv;
//             by += bn * dv;
//         }

//         const TF Z_l = TF( p - 1 ) * p * ( 2 * p - 1 ) / 12 + TF( p ) * TF( p ) / 4;
//         const TF lam = ( by - ay ) / p;

//         E[ p ] = ( ay - Z_l * lam ) / p;
//         E[ p - 1 ] = E[ p ] + lam / 2 - diu( p );
//         for( PI i = p - 1; i--; )
//             E[ i ] = 2 * E[ i + 1 ] - E[ i + 2 ] + lam - diu( i + 1 );

//         der_values.push_back( std::move( E ) );
//         der_primitives.push_back( _primitive_of( der_values.back() ) );
//         return;
//     }

//     for( PI i = 0; i < n; ++i ) {
//         TF e = 0;
//         if ( i && i + 1 < n )
//             e -= prev_der_values[ i ];
//         if ( i )
//             e += prev_der_values[ i - 1 ];
//         if ( i + 1 < n )
//             e += prev_der_values[ i + 1 ];
//         E[ i ] = nder * e;
//     }
//     der_values.push_back( system.solve_using_ldlt( E, 0 ) );
//     der_primitives.push_back( _primitive_of( der_values.back() ) );
// }

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

// DTP TF UTP::integral( TF x0, TF x1, PI num_der_lag_ratio ) {
//     return primitive( x1, num_der_lag_ratio ) - primitive( x0, num_der_lag_ratio );
// }


// DTP TF UTP::x_integral( TF x0, TF x1 ) const {
//     return x_primitive( x1 ) - x_primitive( x0 );
// }

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
    for( TF p : positions )
        fs << p << ", ";    
    fs << " ], [";
    for( TF v : values )
        fs << v << ", ";    
    fs << "] )\n";
}

#undef DTP
#undef UTP

} // namespace usdot
