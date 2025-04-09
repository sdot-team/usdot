#pragma once

#include <tl/support/operators/sum.h>
#include "ConvolutedImage.h"
#include <tl/support/P.h>

template<class TF,class Cl>
class ConvolutedImageDensityIterator : public DensityIterator<TF> {
public:
    virtual void barycenter( TF &pint, TF &area, TF l0, TF l1 ) const override {
        TODO;
    }
    
    virtual TF integral( TF l0, TF l1 ) const override {
        TF res = 0;
        for( PI x = 0; x < cl->h.size(); ++x ) {
            const TF vh = cl->h[ x ];
            const TF x0 = x + 0;
            const TF x1 = x + 1;

            res += vh * cl->int_conv.primitive( l0 - x1 );
            res -= vh * cl->int_conv.primitive( l0 - x0 );
            res += vh * cl->int_conv.primitive( l1 - x0 );
            res -= vh * cl->int_conv.primitive( l1 - x1 );
        }
        return res;
    }

    virtual void display( Displayer &ds ) const override {
        DS_OBJECT( ConvolutedImageDensityIterator, cl );
    }
    
    virtual TF value( TF pos ) const override {
        TF res = 0;
        for( PI x = 0; x < cl->h.size(); ++x ) {
            const TF vh = cl->h[ x ];
            const TF x0 = x + 0;
            const TF x1 = x + 1;

            res += vh * ( 
                + cl->int_conv.value( pos - x0 )
                - cl->int_conv.value( pos - x1 ) 
            );
        }
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

#define DTP template<class TF,class IntConv>
#define UTP ConvolutedImage<TF,IntConv>

DTP UTP::ConvolutedImage( Vec<TF> h, IntConv int_conv ) : int_conv( int_conv ), h( h ) {
}

DTP RcPtr<Density<TF>> UTP::convoluted( TF std ) const {
    //return new ConvolutedImage<TF>{ x0, x1, h, this->std + std };
    TODO;
    return const_cast<ConvolutedImage *>( this );
}

DTP RcPtr<DensityIterator<TF>> UTP::iterator() const {
    auto *res = new ConvolutedImageDensityIterator<TF,UTP>;
    const TF d = int_conv.numerical_width();
    res->x1 = h.size() + d;
    res->x0 = - d;
    res->cl = this;
    return res;
}

DTP void UTP::display( Displayer &ds ) const {
    DS_OBJECT( ConvolutedImage, int_conv, h );
}

DTP TF UTP::mass() const {
    return sum( h );
}

#undef DTP
#undef UTP
