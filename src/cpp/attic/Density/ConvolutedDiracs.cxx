#pragma once

#include <tl/support/operators/max.h>
#include <tl/support/operators/min.h>
#include <tl/support/operators/sum.h>
#include "ConvolutedDiracs.h"
// #include <tl/support/P.h>

namespace usdot {

template<class TF,class Cl>
class ConvolutedDiracsDensityIterator : public DensityIterator<TF> {
public:
    virtual void barycenter( TF &pint, TF &area, TF l0, TF l1 ) const override {
        TODO;
    }
    
    virtual TF integral( TF l0, TF l1 ) const override {
        TF res = 0;
        for( PI i = 0; i < cl->positions.size(); ++i )
            res += cl->weights[ i ] * cl->convolution->dirac_integral( l0, l1, cl->positions[ i ] );
        return res;
    }

    virtual void display( Displayer &ds ) const override {
        DS_OBJECT( ConvolutedDiracsDensityIterator, cl );
    }
    
    virtual TF value( TF pos ) const override {
        TF res = 0;
        for( PI i = 0; i < cl->positions.size(); ++i )
            res += cl->weights[ i ] * cl->convolution->dirac_value( pos, cl->positions[ i ] );
        return res;
    }

    virtual bool move_backward() override {
        return false;
    }

    virtual bool move_forward() override {
        return false;
    }

    const Cl *cl;
};

#define DTP template<class TF>
#define UTP ConvolutedDiracs<TF>

DTP UTP::ConvolutedDiracs( Vec<TF> positions, Vec<TF> weights, RcPtr<Convo> convolution ) : convolution( convolution ), positions( positions ), weights( weights ) {
}

DTP RcPtr<Density<TF>> UTP::convoluted( RcPtr<Convolution<TF>> convolution ) const {
    return new ConvolutedDiracs( positions, weights, this->convolution->convoluted( convolution ) );
}

DTP RcPtr<DensityIterator<TF>> UTP::iterator() const {
    auto *res = new ConvolutedDiracsDensityIterator<TF,UTP>;
    const TF d = convolution->numerical_width();
    res->x1 = max( positions ) + d;
    res->x0 = min( positions ) - d;
    res->cl = this;
    return res;
}

DTP void UTP::display( Displayer &ds ) const {
    DS_OBJECT( ConvolutedDiracs, convolution, positions, weights );
}

DTP TF UTP::min_x( TF eps ) const {
    return min( positions ) - convolution->numerical_width( eps );
}

DTP TF UTP::max_x( TF eps ) const {
    return max( positions ) + convolution->numerical_width( eps );
}

DTP TF UTP::mass() const {
    return sum( weights );
}

#undef DTP
#undef UTP

} // namespace usdot
