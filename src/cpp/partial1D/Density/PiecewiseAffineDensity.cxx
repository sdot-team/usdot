#pragma once

#include "PiecewiseAffineDensity.h"
#include <tl/support/ASSERT.h>

namespace usdot {

template<class TF>
class PiecewiseAffineDensityIterator : public DensityIterator<TF> {
public:
    virtual void barycenter( TF &pint, TF &area, TF l0, TF l1 ) const override {
        using namespace std;
        // v = y0 - y1 * x0 / ( x1 - x0 ) * x 
        //   + y1 / ( x1 - x0 ) * x^2
        // i = ( y0 - y1 * x0 / ( x1 - x0 ) ) * ( l1^2 - l0^2 ) / 2
        //   + y1 / ( x1 - x0 ) * ( l1^3 - l0^3 ) / 3
        // i = ( y0 - dy * x0 ) * ( l1^2 - l0^2 ) / 2
        //   + dy * ( l1^3 - l0^3 ) / 3
        const TF x0 = cl->xs[ i + 0 ];
        const TF x1 = cl->xs[ i + 1 ];
        const TF y0 = cl->ys[ i ][ 0 ];
        const TF y1 = cl->ys[ i ][ 1 ];
        const TF dy = y1 / ( x1 - x0 );
        if ( x1 == x0 )
            return;
        area += ( l1 - l0 ) * ( y0 + dy * ( ( l1 + l0 ) / 2 - x0 ) );
        pint += ( y0 - dy * x0 ) * ( pow( l1, 2 ) - pow( l0, 2 ) ) / 2
              + dy * ( pow( l1, 3 ) - pow( l0, 3 ) ) / 3;
    }
    
    virtual TF integral( TF l0, TF l1 ) const override {
        // v = y0 + y1 * ( x - x0 ) / ( x1 - x0 )
        // v = y0 - y1 * x0 / ( x1 - x0 ) + x * y1 / ( x1 - x0 )
        // p = ( y0 - y1 * x0 / ( x1 - x0 ) ) * x + y1 / ( x1 - x0 ) * x^2 / 2
        // i = ( y0 - y1 * x0 / ( x1 - x0 ) ) * ( l1 - l0 ) + y1 / ( x1 - x0 ) * ( l1^2 - l0^2 ) / 2
        // i = ( y0 - y1 * x0 / ( x1 - x0 ) ) * ( l1 - l0 ) + y1 / ( x1 - x0 ) * ( l1 - l0 ) * ( l1 + l0 ) / 2
        // i = ( l1 - l0 ) * ( y0 + y1 / ( x1 - x0 ) * ( ( l1 + l0 ) / 2 - x0 ) )
        const TF x0 = cl->xs[ i + 0 ];
        const TF x1 = cl->xs[ i + 1 ];
        const TF y0 = cl->ys[ i ][ 0 ];
        const TF y1 = cl->ys[ i ][ 1 ];
        return x1 == x0 ? 0 : ( l1 - l0 ) * ( y0 + y1 * ( ( l1 + l0 ) / 2 - x0 ) / ( x1 - x0 ) );
    }

    virtual void display( Displayer &ds ) const override {
        DS_OBJECT( PiecewiseAffineDensityIterator, cl, this->x0, this->x1, i );
    }
    
    virtual TF value( TF pos ) const override {
        return cl->ys[ i ][ 0 ] + cl->ys[ i ][ 1 ] * ( pos - cl->xs[ i + 0 ] ) / ( cl->xs[ i + 1 ] - cl->xs[ i + 0 ] );
    }

    virtual bool move_backward() override {
        if ( i == 0 )
            return false;
        --i;

        this->x0 = cl->xs[ i + 0 ];
        this->x1 = cl->xs[ i + 1 ];
        return true;
    }

    virtual bool move_forward() override {
        if ( i + 1 >= cl->ys.size() )
            return false;
        ++i;

        this->x0 = cl->xs[ i + 0 ];
        this->x1 = cl->xs[ i + 1 ];
        return true;
    }

    const PiecewiseAffineDensity<TF> *cl;
    PI i;
};

#define DTP template<class TF>
#define UTP PiecewiseAffineDensity<TF>

DTP UTP::PiecewiseAffineDensity( const Vec<TF> &xs, const Vec<Vec<TF,2>> &ys ) : xs( xs ), ys( ys ) {
    ASSERT( xs.size() == ys.size() + 1 );
}

DTP CdfApproximation<TF> UTP::cdf_approximation( TF epsilon ) const {
    Vec<TF> rxs( FromReservationSize(), ys.size() );
    Vec<TF> rys( FromReservationSize(), ys.size() );
    Vec<TF> rzs( FromReservationSize(), ys.size() );

    rxs << xs[ 0 ];
    rys << 0;

    TF y = 0;
    for( PI i = 0; i < ys.size(); ++i ) {
        const TF y0 = ys[ i ][ 0 ];
        const TF y1 = ys[ i ][ 1 ];
        const TF x0 = xs[ i + 0 ];
        const TF x1 = xs[ i + 1 ];
        const TF m = y1 / 3;

        y += ( x1 - x0 ) * ( y0 + y1 / 2 );

        rxs << xs[ i + 1 ];
        rys << y;
        rzs << m;
    }

    P( rxs, rys, rzs );

    return { rxs, rys, rzs };
}

DTP RcPtr<Density<TF>> UTP::convoluted( RcPtr<Convolution<TF>> convolution ) const {
    //return new PiecewiseAffineDensity<TF>{ x0, x1, h, this->std + std };
    TODO;
    return const_cast<PiecewiseAffineDensity *>( this );
}

DTP RcPtr<DensityIterator<TF>> UTP::iterator() const {
    auto *res = new PiecewiseAffineDensityIterator<TF>;
    ASSERT( xs.size() >= 2 );
    res->x0 = xs[ 0 ];
    res->x1 = xs[ 1 ];
    res->cl = this;
    res->i = 0;
    return res;
}

DTP void UTP::display( Displayer &ds ) const {
    DS_OBJECT( PiecewiseAffineDensity, xs, ys );
}


DTP TF UTP::min_x( TF eps ) const {
    return xs.front();
}

DTP TF UTP::max_x( TF eps ) const {
    return xs.back();
}

DTP TF UTP::mass() const {
    TF res = 0;
    for( PI i = 0; i < ys.size(); ++i )
        res += ( xs[ i + 1 ] - xs[ i + 0 ] ) * ( ys[ i ][ 0 ] + ys[ i ][ 1 ] / 2 );
    return res;
}

} // namespace usdot

#undef DTP
#undef UTP
