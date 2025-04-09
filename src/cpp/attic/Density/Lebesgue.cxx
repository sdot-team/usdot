#pragma once

#include "ConvolutedLebesgue.h"
#include "Lebesgue.h"

namespace usdot {


template<class TF>
class LebesgueDensityIterator : public DensityIterator<TF> {
public:
    virtual void barycenter( TF &pint, TF &area, TF l0, TF l1 ) const override {
        const TF a = ( l1 - l0 ) * h;
        pint += a * ( l1 + l0 ) / 2;
        area += a;    
    }
    
    virtual TF integral( TF l0, TF l1 ) const override {
        return ( l1 - l0 ) * h;        
    }

    virtual void display( Displayer &ds ) const override {
        DS_OBJECT( LebesgueDensityIterator, this->x0, this->x1, h );
    }
    
    virtual TF value( TF pos ) const override {
        return h;
    }

    virtual bool move_backward() override {
        return false;
    }

    virtual bool move_forward() override {
        return false;
    }

    TF h; ///<
};

#define DTP template<class TF>
#define UTP Lebesgue<TF>

DTP UTP::Lebesgue( TF x0, TF x1, TF h ) : x0( x0 ), x1( x1 ), h( h ) {
}

DTP CdfApproximation<TF> UTP::cdf_approximation( TF epsilon ) const {
    Vec<TF> rys{ 0, ( x1 - x0 ) * h };
    Vec<TF> rxs{ x0, x1 };
    Vec<TF> rzs{ 0, 0 };

    return { rxs, rys, rzs };
}

DTP RcPtr<Density<TF>> UTP::convoluted( RcPtr<Convolution<TF>> convolution ) const {
    return new ConvolutedLebesgue{ x0, x1, h, convolution };
}

DTP RcPtr<DensityIterator<TF>> UTP::iterator() const {
    if ( x0 == x1 )
        return {};

    auto *res = new LebesgueDensityIterator<TF>;
    if ( x0 > x1 ) {
        res->x0 = x1;
        res->x1 = x0;
        res->h = -h;
    } else {
        res->x0 = x0;
        res->x1 = x1;
        res->h = h;
    }
    return res;
}

DTP TF UTP::min_x( TF eps ) const {
    return x0;
}

DTP TF UTP::max_x( TF eps ) const {
    return x1;
}

DTP void UTP::display( Displayer &ds ) const {
    DS_OBJECT( Lebesgue, x0, x1, h );
}

DTP TF UTP::mass() const {
    return ( x1 - x0 ) * h;
}

#undef DTP
#undef UTP


} // namespace usdot
