#pragma once

#include "ConvolutedPiecewiseConstantDensity.h"
#include <tl/support/ASSERT.h>
#include <tl/support/P.h>

namespace usdot {

template<class TF>
class ConvolutedPiecewiseConstantDensityIterator : public DensityIterator<TF> {
public:
    virtual void display( Displayer &ds ) const override {
        DS_OBJECT( ConvolutedPiecewiseConstantDensityIterator, cl );
    }

    virtual void barycenter( TF &pint, TF &area, TF l0, TF l1 ) const override {
        // TF a = ( l1 - l0 ) * cl->ys[ i ];
        // pint += a * ( l1 + l0 ) / 2;
        // area += a;
        TODO;
    }
    
    virtual TF integral( TF l0, TF l1 ) const override {
        TF res = 0;
        for( PI i = 0; i < cl->ys.size(); ++i )
            res += cl->ys[ i ] * cl->convolution->step_integral( l0, l1, cl->xs[ i + 0 ], cl->xs[ i + 1 ] );
        return res;
    }

    virtual TF value( TF pos ) const override {
        TF res = 0;
        for( PI i = 0; i < cl->ys.size(); ++i )
            res += cl->ys[ i ] * cl->convolution->step_value( pos, cl->xs[ i + 0 ], cl->xs[ i + 1 ] );
        return res;
    }

    virtual bool move_backward() override {
        return false;
    }

    virtual bool move_forward() override {
        return false;
    }

    const ConvolutedPiecewiseConstantDensity<TF> *cl;
};

#define DTP template<class TF>
#define UTP ConvolutedPiecewiseConstantDensity<TF>

DTP UTP::ConvolutedPiecewiseConstantDensity( const Vec<TF> &xs, const Vec<TF> &ys, RcPtr<Convo> convolution ) : convolution( convolution ), xs( xs ), ys( ys ) {
    ASSERT( xs.size() == ys.size() + 1 );
    ASSERT( xs.size() >= 2 );
}

DTP RcPtr<Density<TF>> UTP::convoluted( RcPtr<Convolution<TF>> convolution ) const {
    return new ConvolutedPiecewiseConstantDensity( xs, ys, this->convolution->convoluted( convolution ) );
}

DTP RcPtr<DensityIterator<TF>> UTP::iterator() const {
    auto *res = new ConvolutedPiecewiseConstantDensityIterator<TF>;
    res->x0 = min_x();
    res->x1 = max_x();
    res->cl = this;
    return res;
}

DTP void UTP::display( Displayer &ds ) const {
    DS_OBJECT( ConvolutedPiecewiseConstantDensity, convolution, xs, ys );
}

DTP TF UTP::min_x( TF eps ) const {
    return xs.front() - convolution->numerical_width( eps );
}

DTP TF UTP::max_x( TF eps ) const {
    return xs.back() + convolution->numerical_width( eps );
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
