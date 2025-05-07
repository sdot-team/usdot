#pragma once

#include "utility/linspace.h"
#include "ConvolutedDensity.h"
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace usdot {
    
#define DTP template<class TF,int regular_positions>
#define UTP ConvolutedDensity<TF,regular_positions>

DTP UTP::ConvolutedDensity( TF beg_original_positions, TF end_original_positions, const VF &original_values ) {
    this->beg_original_positions = beg_original_positions;
    this->end_original_positions = end_original_positions;
    this->original_values = original_values;
    current_lag_ratio = -1;

    if ( regular_positions == 0 )
        this->original_positions = linspace<TF>( beg_original_positions, end_original_positions, original_values.size() );

    // normalization
    TF o = 0;
    for( PI i = 0; i + 1 < this->original_values.size(); ++i )
        o += ( position( i + 1 ) - position( i + 0 ) ) * ( this->original_values[ i + 0 ] + this->original_values[ i + 1 ] ) / 2;
    for( auto &v : this->original_values )
        v /= o;
}

DTP UTP::ConvolutedDensity( const VF &original_values ) : ConvolutedDensity( 0, original_values.size() - 1, original_values ) {
}

DTP UTP::ConvolutedDensity( const VF &original_positions, const VF &original_values ) {
    this->beg_original_positions = original_positions.front();
    this->end_original_positions = original_positions.back();
    this->original_positions = original_positions;
    this->original_values = original_values;
    current_lag_ratio = -1;

    // normalization
    TF o = 0;
    for( PI i = 0; i + 1 < this->original_values.size(); ++i )
        o += ( position( i + 1 ) - position( i + 0 ) ) * ( this->original_values[ i + 0 ] + this->original_values[ i + 1 ] ) / 2;
    for( auto &v : this->original_values )
        v /= o;
}

DTP void UTP::set_lag_ratio( TF t ) {
    using namespace std;

    if ( current_lag_ratio == t )
        return;
    current_lag_ratio = t;

    //
    _set_values( t );

    // primitives
    x_primitives.resize( this->values.size() );
    primitives.resize( this->values.size() );
    TF o = 0, x_o = 0;
    for( PI i = 0; i + 1 < this->values.size(); ++i ) {
        const TF v0 = this->values[ i + 0 ]; 
        const TF v1 = this->values[ i + 1 ];
        const TF x0 = position( i + 0 ); 
        const TF x1 = position( i + 1 ); 
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
    for( TF &v : this->x_primitives )
        v /= o;
    for( TF &v : this->primitives )
        v /= o;
    for( TF &v : this->values )
        v /= o;

    // opt_pos
    if ( regular_positions == 0 ) {
        const PI opt_pos_size = this->values.size() - 1;
        opt_pos_begs.resize( opt_pos_size, numeric_limits<PI>::max() );
        opt_pos_beg = beg_original_positions;
        opt_pos_end = end_original_positions;
        opt_pos_coeff = opt_pos_size / ( opt_pos_end - opt_pos_beg );
        for( PI i = 0; i + 1 < this->values.size(); ++i ) {
            const PI x0 = floor( ( this->positions[ i + 0 ] - opt_pos_beg ) * opt_pos_coeff );
            const PI x1 = ceil( ( this->positions[ i + 1 ] - opt_pos_beg ) * opt_pos_coeff );
            for( PI n = x0; n < min( x1, opt_pos_size ); ++n )
                opt_pos_begs[ n ] = min( opt_pos_begs[ n ], i );
        }
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

DTP TF UTP::derivative( TF x ) const {
    if ( regular_positions ) {
        if ( regular_positions == 1 )
            x = ( x - opt_pos_beg ) * opt_pos_coeff;
        if ( x < 0 || x >= values.size() - 1 )
            return 0;

        if ( regular_positions == 1 )
            throw std::runtime_error( "TODO" );

        PI i( x );
        return values[ i + 1 ] - values[ i ];
    }

    throw std::runtime_error( "TODO" );
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

DTP TF UTP::derivative( TF x, PI num_der_lag_ratio ) {
    throw std::runtime_error( "TODO" );
}

DTP UTP::VF UTP::_primitive_of( const VF &values ) const {
    const PI n = values.size();
    VF res( n );

    TF o = 0;
    for( PI i = 0; i + 1 < n; ++i ) {
        const TF x0 = position( i + 0 ); 
        const TF x1 = position( i + 1 ); 
        const TF v0 = values[ i + 0 ]; 
        const TF v1 = values[ i + 1 ];

        res[ i ] = o;

        o += ( x1 - x0 ) * ( v0 + v1 ) / 2;
    }
    res.back() = o;

    return res;
}

DTP void UTP::_set_values( TF t ) {
    const PI n = original_values.size();
    const PI p = n - 1;

    //
    positions = original_positions;
    der_values.clear();

    // values
    if ( t == 0 ) {
        values = original_values;
        return;
    }

    if ( t == 1 ) {
        values.resize( n );
        for( PI i = 0; i < n; ++i )
            values[ i ] = 1 / TF( n - 1 );
        return;
    }

    //
    system = { n };
    VF E( n );
    for( PI i = 1; i < p; ++i ) {
        system.tridiag_value( i, i - 1 ) = - t;
        system.tridiag_value( i, i ) = 1 + t;
        system.line_value( i ) = 1;
    }
    system.tridiag_value( p, p - 1 ) = - t;
    system.tridiag_value( p, p ) = 1;
    system.tridiag_value( 0, 0 ) = 1;
    system.line_value( 0 ) = 0.5;
    system.line_value( p ) = 0.5;
    for( PI i = 0; i < n; ++i )
        E[ i ] = ( 1 - t ) * original_values[ i ];
    system.line_value( n ) = 0;

    system.inplace_ldlt_decomposition();
    values = system.solve_using_ldlt( E, 1 );
}

DTP void UTP::_append_der_value() {
    const PI n = values.size();
    const PI p = n - 1;

    if ( regular_positions != 2 )
        throw std::runtime_error( "TODO" );

    // first derivative ?
    VF E( n );
    if ( der_values.empty() ) {
        if ( current_lag_ratio == 1 ) {
            // [   1,  -1,   0,   0,   0,   0,   0,   0,   0,   0, 0.5 ]
            // [  -1,   2,  -1,   0,   0,   0,   0,   0,   0,   0,   1 ]
            // [   0,  -1,   2,  -1,   0,   0,   0,   0,   0,   0,   1 ]
            // [   0,   0,  -1,   2,  -1,   0,   0,   0,   0,   0,   1 ]
            // [   0,   0,   0,  -1,   2,  -1,   0,   0,   0,   0,   1 ]
            // [   0,   0,   0,   0,  -1,   2,  -1,   0,   0,   0,   1 ]
            // [   0,   0,   0,   0,   0,  -1,   2,  -1,   0,   0,   1 ]
            // [   0,   0,   0,   0,   0,   0,  -1,   2,  -1,   0,   1 ]
            // [   0,   0,   0,   0,   0,   0,   0,  -1,   2,  -1,   1 ]
            // [   0,   0,   0,   0,   0,   0,   0,   0,  -1,   1, 0.5 ]
            // [ 0.5,   1,   1,   1,   1,   1,   1,   1,   1, 0.5,   0 ]
            TF ay = 0, by = 0;
            for( PI i = 0; i < n; ++i ) {
                const TF dv = values[ i ] - original_values[ i ];
                const TF an = TF( i ) * TF( i ) / 2;
                const TF bn = an + 1;
                ay += an * dv;
                by += bn * dv;
            }

            const TF Z_l = TF( p - 1 ) * p * ( 2 * p - 1 ) / 12 + TF( p ) * TF( p ) / 4;
            const TF lam = ( by - ay ) / p;

            E[ p ] = ( ay - Z_l * lam ) / p;
            E[ p - 1 ] = E[ p ] + lam / 2 + original_values[ p ] - values[ p ];
            for( PI i = p - 1; i--; )
                E[ i ] = 2 * E[ i + 1 ] - E[ i + 2 ] + lam + original_values[ i + 1 ] - values[ i + 1 ];

            der_values.push_back( std::move( E ) );
            der_primitives.push_back( _primitive_of( der_values.back() ) );
            return;
        }

        for( PI i = 0; i < n; ++i ) {
            TF e = - original_values[ i ];
            if ( i && i + 1 < n )
                e -= values[ i ];
            if ( i )
                e += values[ i - 1 ];
            if ( i + 1 < n )
                e += values[ i + 1 ];
            E[ i ] = e;
        }
        der_values.push_back( system.solve_using_ldlt( E, 0 ) );
        der_primitives.push_back( _primitive_of( der_values.back() ) );
        return;
    }

    // 
    const auto &prev_der_values = der_values.back();
    const PI nder = der_values.size() + 1;
    if ( current_lag_ratio == 1 ) {
        auto diu = [&]( PI i ) {
            TF dv = 0;
            if ( i && i + 1 < n )
                dv -= prev_der_values[ i ];
            if ( i )
                dv += prev_der_values[ i - 1 ];
            if ( i + 1 < n )
                dv += prev_der_values[ i + 1 ];
            return nder * dv;
        };

        TF ay = 0, by = 0;
        for( PI i = 0; i < n; ++i ) {
            const TF an = TF( i ) * TF( i ) / 2;
            const TF bn = an + 1;
            const TF dv = diu( i );
            ay += an * dv;
            by += bn * dv;
        }

        const TF Z_l = TF( p - 1 ) * p * ( 2 * p - 1 ) / 12 + TF( p ) * TF( p ) / 4;
        const TF lam = ( by - ay ) / p;

        E[ p ] = ( ay - Z_l * lam ) / p;
        E[ p - 1 ] = E[ p ] + lam / 2 - diu( p );
        for( PI i = p - 1; i--; )
            E[ i ] = 2 * E[ i + 1 ] - E[ i + 2 ] + lam - diu( i + 1 );

        der_values.push_back( std::move( E ) );
        der_primitives.push_back( _primitive_of( der_values.back() ) );
        return;
    }

    for( PI i = 0; i < n; ++i ) {
        TF e = 0;
        if ( i && i + 1 < n )
            e -= prev_der_values[ i ];
        if ( i )
            e += prev_der_values[ i - 1 ];
        if ( i + 1 < n )
            e += prev_der_values[ i + 1 ];
        E[ i ] = nder * e;
    }
    der_values.push_back( system.solve_using_ldlt( E, 0 ) );
    der_primitives.push_back( _primitive_of( der_values.back() ) );
}

DTP TF UTP::value( TF x, PI num_der_lag_ratio ) {
    // not a derivative ? 
    if ( num_der_lag_ratio == 0 )
        return value( x );

    // precomputation(s)
    while ( der_values.size() < num_der_lag_ratio )
        _append_der_value();

    //
    const VF &values = der_values[ num_der_lag_ratio - 1 ];
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

DTP TF UTP::primitive( TF x, PI num_der_lag_ratio ) {
    if ( num_der_lag_ratio == 0 )
        return primitive( x );

    // precomputation(s)
    while ( der_values.size() < num_der_lag_ratio )
        _append_der_value();

    const auto &primitives = der_primitives[ num_der_lag_ratio - 1 ];
    const auto &values = der_values[ num_der_lag_ratio - 1 ];
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

DTP TF UTP::integral( TF x0, TF x1, PI num_der_lag_ratio ) {
    return primitive( x1, num_der_lag_ratio ) - primitive( x0, num_der_lag_ratio );
}


DTP TF UTP::position( PI i ) const {
    if ( regular_positions == 2 )
        return i;
    if ( regular_positions == 1 )
        return opt_pos_beg + i / opt_pos_coeff;
    return positions[ i ];
}

DTP TF UTP::x_integral( TF x0, TF x1 ) const {
    return x_primitive( x1 ) - x_primitive( x0 );
}

DTP TF UTP::integral( TF x0, TF x1 ) const {
    return primitive( x1 ) - primitive( x0 );
}

DTP TF UTP::min_x() const {
    return beg_original_positions;
}

DTP TF UTP::max_x() const {
    return end_original_positions;
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
