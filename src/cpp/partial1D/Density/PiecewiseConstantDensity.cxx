#pragma once

#include "ConvolutedPiecewiseConstantDensity.h"
#include "PiecewiseConstantDensity.h"
#include <tl/support/ASSERT.h>

namespace usdot {

template<class TF>
class PiecewiseConstantDensityIterator : public DensityIterator<TF> {
public:
    virtual void barycenter( TF &pint, TF &area, TF l0, TF l1 ) const override {
        TF a = ( l1 - l0 ) * cl->ys[ i ];
        pint += a * ( l1 + l0 ) / 2;
        area += a;
    }
    
    virtual TF integral( TF l0, TF l1 ) const override {
        return ( l1 - l0 ) * cl->ys[ i ];
    }

    virtual void display( Displayer &ds ) const override {
        DS_OBJECT( PiecewiseConstantDensityIterator, cl, this->x0, this->x1, i );
    }
    
    virtual TF value( TF pos ) const override {
        return cl->ys[ i ];
    }

    virtual TF inv_cdf( TF dy ) const override {
        return this->x0 + dy / cl->ys[ i ];
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

    const PiecewiseConstantDensity<TF> *cl;
    PI i;
};

#define DTP template<class TF>
#define UTP PiecewiseConstantDensity<TF>

DTP UTP::PiecewiseConstantDensity( const Vec<TF> &xs, const Vec<TF> &ys ) : xs( xs ), ys( ys ) {
    ASSERT( xs.size() == ys.size() + 1 );
}

DTP CdfApproximation<TF> UTP::cdf_approximation( TF epsilon ) const {
    Vec<TF> rxs( FromReservationSize(), ys.size() );
    Vec<TF> rys( FromReservationSize(), ys.size() );
    Vec<TF> rzs( FromReservationSize(), ys.size() );
    auto app = [&]( TF x, TF y ) {
        if ( rxs.back() != x || rys.back() != y ) {
            rxs << x;
            rys << y;
            rzs << 0;
        }
    };

    rxs << xs[ 0 ];
    rys << 0;

    TF y = 0;
    for( PI i = 0; i < ys.size(); ++i ) {
        y += ( xs[ i + 1 ] - xs[ i + 0 ] ) * ys[ i ];
        app( xs[ i + 1 ], y );
    }


    return { rxs, rys, rzs };
}

DTP RcPtr<Density<TF>> UTP::convoluted( RcPtr<Convolution<TF>> convolution ) const {
    return new ConvolutedPiecewiseConstantDensity<TF>( xs, ys, convolution );
}

DTP RcPtr<DensityIterator<TF>> UTP::iterator() const {
    auto *res = new PiecewiseConstantDensityIterator<TF>;
    ASSERT( xs.size() >= 2 );
    res->x0 = xs[ 0 ];
    res->x1 = xs[ 1 ];
    res->cl = this;
    res->i = 0;
    return res;
}

DTP void UTP::display( Displayer &ds ) const {
    DS_OBJECT( PiecewiseConstantDensity, xs, ys );
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
        res += ( xs[ i + 1 ] - xs[ i + 0 ] ) * ys[ i ];
    return res;
}

#undef DTP
#undef UTP

} // namespace usdot
